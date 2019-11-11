subroutine pointimplicit(dtmp,dprs,aYi,delt)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!     i(subroutine)-m(multi)-t(time)-s(scale)-s(stanford)
!
!     this subroutine was created by kazushige maeda, 2014/05/30
!                         written for f90 by hiroaki watanabe, 2015/10/03
!
!     MTS method including 1st-order euler integrator
!     Reference: Gou, X., Sun, W., Chen, Z., Ju, Y.,
!                "A dynamic multi-timescale method for combustion modeling
!                 with detailed and reduced chemical kinetic mechanisms,"
!                Combustion and Flame, 157 (2010) 1111-1121
!
!     Stanford model
!     Reference: Hong, Z., Davidson, D. F., Hanson, R. K.,
!                "An improved H2/O2 mechanism based on recent shocktube/
!                 laser absorption measurements,"
!                Combustion and Flame, 158 (2011) 633-644
!
!     Units    : second, mole, cubic centimeter, calory, Kelvin
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use mod_globals
  use mod_parameter, only : imax, ldim, lr, lu, lv, le, l1, ls1, lp1, l2, ls2,     &
       &                         lp2, l3, ls3, lp3, l4, ls4, lp4, l5, ls5, lp5, l6,&
       &                         ls6, lp6, l7, ls7, lp7, l8, ls8, lp8, l9, ls9,    &
       &                         lp9, lfl, lpr, lpe, lpp, a05, a12, a13, a14,      &
       &                         a15, dru, dplt
  use mod_function
  implicit none
  real(8),parameter :: alim = 1.d0 - 1.d-13
  real(8), intent(inout) :: dtmp, dprs, aYi(9)
  real(8), intent(in) :: delt
  real(8), allocatable :: adns(:), afm(:), ab(:), dml(:), dmlr(:), dmlru(:), dplh(:,:,:), atmw(:), adnsr(:), arhs(:)

  ! number of reaction constant
  integer, parameter :: lsp = 9
  integer, parameter :: lnr = 32
  integer, parameter :: len = 25
  integer, parameter :: ler = 28
  integer, parameter :: lnn = 31

  real(8) :: at, afmt, aeng, afdd, afd, af, addd, adt, druo, phti
  real(8) :: ar, au, av, awk1, awk2, awk3, totaldens, totalh
  real(8) :: atmp, arho, art, atmpr, atmpl, aruo, atime
  real(8) :: ae(lnr), akf(lnr), ax(lsp), ah(lsp), as(lsp), ag(len), ahh(len), &
  &          akb(lnr), akfc(lnn), akbc(lnn), ars(len), adf(ler), akft(len),   &
  &          adb(ler), akbt(len), adtd(lsp), admd(lsp), akc(lsp), aw(lsp),    &
  &          a1diag(8), akfr4(lsp), akfr5(lsp), akbr4(lsp), akbr5(lsp),       &
  &          askrwp(lsp), atrpr4(lsp), atrpr5(lsp), add02, add03, add04,      &
  &          add05, atrb2, atrb3, atrb4, adfg2, adfg3, adfg4, am04, am05,     &
  &          am09, am21, am23, am24, am25, atrp2, atrp3, atrp4, atrp5,        &
  &          atrpl2, atrpl3, atrpl4, adlp2, adlp3, adlp4, atra2, atra3,       &
  &          atra4, atrf2, atrf3, atrf4, adpi2, adpi3, adpi4, adpi5, arr,     &
  &          arf09, arf10, arf21, arf22, arf23, arf24, arf25,                 &
  &          arb09, arb10, arb21, arb22, arb23, arb24, arb25,                 &
  &          asktwc, asktwp, atmaxi, acon, apt, aptr, acn, ajcb,              &
  &          akfr2, akfr3, akbr2, akbr3, askrwc, atrr2, atrr3, atrr4, atrp4r, &
  &          atrp5r, dppd(lsp), dcpo, ahti(lsp)

  integer :: i, j, k, itm
! constant parameter
  real(8), parameter :: arc = 1.9872d0, arp = 82.06d0
  real(8), parameter :: arcr = 1.d0 / arc
  real(8), parameter :: arur = 1.d0 / dru
! limiter
  real(8), parameter :: aemn = -130.d0, aemx = 130.d0
  real(8), parameter :: akmn = 1.d-35, akmx = 1.d+35
  real(8), parameter :: aeps = 1.d-35
! reaction frozen temperature
  real(8), parameter :: a400 = 400.d0
  real(8), parameter :: aerr = 1.d-12
  integer, parameter :: itr = 10 ! number of iteration in Newton method
! constant parameter for Troe (for reduction)
  real(8), parameter :: atrc8 =  -0.4d0 - 0.67d0*log10(0.8d0)
  real(8), parameter :: atrn8 = 0.75d0 - 1.27d0*log10(0.8d0)
  real(8), parameter :: atrc7 = -0.4d0 - 0.67d0*log10(0.7d0)
  real(8), parameter :: atrn7 = 0.75d0 - 1.27d0*log10(0.7d0)
  real(8), parameter :: atrf7 = log10(0.7d0)
  real(8), parameter :: atrf8 = log10(0.8d0)

 !reaction parameters
  !modified Arrhenius equation, k = A * T^n * exp(-E/RT) : acf=A, ant=n, aea=E
  real(8), parameter :: acf(1:lnr) =                                          &
  & (/ 1.04000d+14, 5.59000d+13, 3.70000d+19, 5.59000d+13, 5.69000d+18,       &
  &    5.59000d+13, 2.65000d+19, 8.59000d+14, 9.55000d+15, 1.74000d+12,       &
  &    7.59000d+13, 2.89000d+13, 1.30000d+11, 4.20000d+14, 6.06000d+27,       &
  &    1.00000d+26, 3.57000d+04, 3.82000d+12, 8.79000d+14, 2.17000d+08,       &
  &    7.08000d+13, 1.45000d+12, 3.66000d+06, 1.63000d+13, 1.21000d+07,       &
  &    1.02000d+13, 8.43000d+11, 5.84000d+18, 9.03000d+14, 4.58000d+19,       &
  &    6.16000d+15, 4.71000d+18 /)
  real(8), parameter :: ant(1:lnr) =                                          &
  & (/ 0.00000d+00, 0.20000d+00,-1.00000d+00, 0.20000d+00,-1.10000d+00,       &
  &    0.20000d+00,-1.30000d+00, 0.00000d+00, 0.00000d+00, 0.00000d+00,       &
  &    0.00000d+00, 0.00000d+00, 0.00000d+00, 0.00000d+00,-3.31000d+00,       &
  &   -2.44000d+00, 2.40000d+00, 0.00000d+00, 0.00000d+00, 1.52000d+00,       &
  &    0.00000d+00, 0.00000d+00, 2.08700d+00, 0.00000d+00, 2.00000d+00,       &
  &    0.00000d+00, 0.00000d+00,-1.10000d+00, 0.00000d+00,-1.40000d+00,       &
  &   -0.50000d+00,-1.00000d+00 /)
  real(8), parameter :: aea(1:lnr) =                                          &
  & (/  15286.d+00,      0.d+00,      0.d+00,      0.d+00,      0.d+00,       &
  &         0.d+00,      0.d+00,  48560.d+00,  42203.d+00,    318.d+00,       &
  &      7269.d+00,   -500.d+00,  -1603.d+00,  11980.d+00, 120770.d+00,       &
  &    120160.d+00,  -2111.d+00,   7948.d+00,  19170.d+00,   3457.d+00,       &
  &       300.d+00,      0.d+00,  -1450.d+00,   -445.d+00,   5200.d+00,       &
  &      3577.d+00,   3970.d+00, 104380.d+00,  96070.d+00, 104380.d+00,       &
  &         0.d+00,      0.d+00 /)

  !difference of mole
  integer, parameter :: app(1:len) =                                          &
  & (/ 0, 1, 1, 1,-1, 0, 0, 0,-1,-1,                                          &
  &    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                          &
  &   -1,-1,-1, 1, 1 /)

  !adjust unit
  real(8), parameter :: autf(1:ler) =                                         &
  & (/ 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d+0, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,  &
  &    1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,  &
  &    1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-6, 1.d-6 /)

  real(8), parameter :: autb(1:ler) =                                         &
  & (/ 1.d-3, 1.d+0, 1.d+0, 1.d+0, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,  &
  &    1.d-6, 1.d-6, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3, 1.d-3,  &
  &    1.d-3, 1.d-3, 1.d-3, 1.d-6, 1.d-6, 1.d-6, 1.d-3, 1.d-3 /)

!///////////////////////////////////////////////////////////////////////
!//////////////////////  stanford reaction model   /////////////////////
!///////////////////////////////////////////////////////////////////////
!
!        1 :  h    + o2           = oh   + o
!        2 :  h    + o2   + (h2o) = ho2  + (h2o)        !pressure-dependent
!        3 :  h    + o2   + (o2)  = ho2  + (o2)         !pressure-dependent
!        4 :  h    + o2   + (M)   = ho2  + (M)          !pressure-dependent
!        5 :  h2o2        + (M)   = 2oh  + (M)          !pressure-dependent
!        6 :  oh   + h2o2         = h2o  + ho2          !non-Arrhenius
!        7 :  oh   + ho2          = h2o  + o2
!        8 :  ho2  + ho2          = h2o2 + o2           !non-Arrhenius
!        9 :  h2o         + [M]   = h    + oh  + [M]
!        10:  h2o         + [h2o] = oh   + h   + [h2o]
!        11:  oh   + oh           = h2o  + o
!        12:  o    + h2           = h    + oh           !non-Arrhenius
!        13:  h2   + oh           = h2o  + h
!        14:  h    + ho2          = oh   + oh
!        15:  h    + ho2          = h2o  + o
!        16:  h    + ho2          = h2   + o2
!        17:  o    + ho2          = oh   + o2
!        18:  h2o2 + h            = ho2  + h2
!        19:  h2o2 + h            = h2o  + oh
!        20:  h2o2 + o            = oh   + ho2
!        21:  h2          + [M]   = h    + h   + [M]
!        22:  h2          + [h2]  = h    + h   + [h2]
!        23:  h2          + [M]   = h    + h   + [M]    !M = N2 and O2 case
!        24:  o    + o    + [M]   = o2         + [M]
!        25:  o    + h    + [M]   = oh         + [M]
!
!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////

!----------------------------------------------------------------------
  if (bdbg(1:2) == 'on') write(lwrt,*) 'in ipdys'
!-----------------------------------------------------------------------

! allocate
  allocate(adns(lsp),adnsr(lsp),afm(lsp),ab(6),arhs(lsp))
  allocate(dml(lsp), dmlr(lsp), dmlru(lsp), dplh(7,lsp,3), atmw(lsp))

! chemical_database

  dml(1:lsp) = (/2.016d0,  32.000d0, 1.008d0, 16.000d0, 17.008d0, 18.016d0, &
       &                33.008d0, 34.016d0, 28.016d0 /)

! reciprocal molecular weight
  dmlr(1:lsp)  = 1.d0/dml(1:lsp)
  dmlru(1:lsp) = dmlr(1:lsp)*dru

  dplh(1:7,1:lsp,1:3) = reshape((/                      &
! --  300-1000K --
! H2
    &  2.34433112d+00,  7.98052075d-03, -1.94781510d-05,  2.01572094d-08,       &
    & -7.37611761d-12, -9.17935173d+02,  6.83010238d-01,                        &
! O2
    &  3.78245636d+00, -2.99673416d-03,  9.84730201d-06, -9.68129509d-09,       &
    &  3.24372837d-12, -1.06394356d+03,  3.65767573d+00,                        &
! H
    &  2.50000000d+00,  7.05332819d-13, -1.99591964d-15,  2.30081632d-18,       &
    & -9.27732332d-22,  2.54736599d+04, -4.46682853d-01,                        &
! O
    &  3.16826710d+00, -3.27931884d-03,  6.64306396d-06, -6.12806624d-09,       &
    &  2.11265971d-12,  2.91222592d+04,  2.05193346d+00,                        &
! OH
    &  3.99201543d+00, -2.40131752d-03,  4.61793841d-06, -3.88113333d-09,       &
    &  1.36411470d-12,  3.37227356d+03, -1.03925458d-01,                        &
! H2O
    &  4.19864056d+00, -2.03643410d-03,  6.52040211d-06, -5.48797062d-09,       &
    &  1.77197817d-12, -3.02937267d+04, -8.49032208d-01,                        &
! HO2
    &  4.30179807d+00, -4.74912097d-03,  2.11582905d-05, -2.42763914d-08,       &
    &  9.29225225d-12,  2.64018485d+02,  3.71666220d+00,                        &
! H2O2
    &  4.27611269d+00, -5.42822417d-04,  1.67335701d-05, -2.15770813d-08,       &
    &  8.62454363d-12, -1.77025821d+04,  3.43505074d+00,                        &
! N2
    &  3.298677d+00,    1.4082404d-03,  -3.963222d-06,    5.641515d-09,         &
    & -2.444854d-12,   -1.0208999d+03,   3.950372d+00,                          &

! -- 1000-5000K --
! H2
    &  3.33727920d+00, -4.94024731d-05,  4.99456778d-07, -1.79566394d-10,       &
    &  2.00255376d-14, -9.50158922d+02, -3.20502331d+00,                        &
! O2
    &  3.28253784d+00,  1.48308754d-03, -7.57966669d-07,  2.09470555d-10,       &
    & -2.16717794d-14, -1.08845772d+03,  5.45323129d+00,                        &
! H
    &  2.50000001d+00, -2.30842973d-11,  1.61561948d-14, -4.73515235d-18,       &
    &  4.98197357d-22,  2.54736599d+04, -4.46682914d-01,                        &
! O
    &  2.56942078d+00, -8.59741137d-05,  4.19484589d-08, -1.00177799d-11,       &
    &  1.22833691d-15,  2.92175791d+04,  4.78433864d+00,                        &
! OH
    &  3.09288767d+00,  5.48429716d-04,  1.26505228d-07, -8.79461556d-11,       &
    &  1.17412376d-14,  3.61585000d+03,  4.47669610d+00,                        &
! H2O
    &  3.03399249d+00,  2.17691804d-03, -1.64072518d-07, -9.70419870d-11,       &
    &  1.68200992d-14, -3.00042971d+04,  4.96677010d+00,                        &
! HO2
    &  4.17228741d+00,  1.88117627d-03, -3.46277286d-07,  1.94657549d-11,       &
    &  1.76256905d-16,  3.10206839d+01,  2.95767672d+00,                        &
! H2O2
    &  4.16500285d+00,  4.90831694d-03, -1.90139225d-06,  3.71185986d-10,       &
    & -2.87908305d-14, -1.78617877d+04,  2.91615662d+00,                        &
! N2
    &  2.926640d+00,    1.4879768d-03,  -5.684760d-07,    1.0097038d-10,        &
    & -6.753351d-15,   -9.227977d+02,    5.980528d+00,                          &

! -- perfectgas --
! H2
    &  0.35d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! O2
    &  0.35d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! H
    &  0.25d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! O
    &  0.25d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! OH
    &  0.35d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! H2O
    &  0.40d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! HO2
    &  0.40d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! H2O2
    &  0.35d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                            &
! N2
    &  0.35d+01, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0                             &
    &  /), (/7,lsp,3/))


  atmp = dtmp
  atmpr = 1/dtmp

  atime   = 0.d0
  totalh  = 0.d0
  atmw(:) = aYi(:)*dmlr(:)      ! amlf(:)/atw
  druo    = dru*(sum(atmw(:)))
  itm = 1 + int(a05 + sign(a05,(atmp - dplt)))

  !dprs   = totaldens*druo*atmp
  totaldens = dprs*atmpr/druo

  adns(:) = aYi(:)*totaldens
  arho    = sum(adns(1:lsp-1))
  adnsr   = 1.d0 / arho

  do i = 1, lsp
    call calc_phti(i,itm,dtmp,phti,dmlr,dplh)
    ahti(i) = phti
    totalh   = totalh  + adns(i)*ahti(i)/arho
  end do
  aeng = dprs - totalh*arho
  art = arcr*atmpr

!     [E/R]
!     activation energy

      ae(1:lnr)   = min(aemx,max(aemn,-aea(1:lnr)*art))

!     [kf]
!     forward reaction rate constant
!     Arrhenius-s form : kf=AT**n exp(-Ea/RT)

      akf(1:lnr)  = acf(1:lnr)*atmp**ant(1:lnr)*exp(ae(1:lnr))

!     [Xi]
!     concentration of species Xi
!     mole density

      ax(1:lsp)   = adns(1:lsp) * dmlr(1:lsp)

!     [M]
!     Third Body

      am04  = ax(1)*1.5d0  +ax(2)*0.d0  +ax(3)  +ax(4)  +ax(5)  +ax(6)*0.d0  +ax(7)         &
  &         + ax(8)      +ax(9)
      am05  = ax(1)        +ax(2)      +ax(3)  +ax(4)  +ax(5)  +ax(6)*9.0d0  +ax(7)         &
  &         + ax(8)      +ax(9)*1.5d0
      am09  = ax(1)*3.0d0  +ax(2)*1.5d0  +ax(3)  +ax(4)  +ax(5)  +ax(6)*0.0d0              &
  &         + ax(7)  + ax(8)    +ax(9)*2.0d0
      am21  = ax(1)*0.0d0  +ax(2)*0.0d0  +ax(3)  +ax(4)  +ax(5)  +ax(6)*14.4d0             &
  &         + ax(7) + ax(8)    +ax(9)*0.0d0
      am23  = ax(2)        +ax(9)
      am24  = ax(1)*2.5d0  +ax(2)      +ax(3)  +ax(4)  +ax(5)  +ax(6)*12d0   +ax(7)         &
  &         + ax(8)      +ax(9)
      am25  = ax(1)*2.5d0  +ax(2)      +ax(3)  +ax(4)  +ax(5)  +ax(6)*12d0   +ax(7)         &
  &         + ax(8)      +ax(9)

!-----------------------------------------------------------
!  the reaction rate constant of #02

      atrp2  = akf(3)*ax(6)*1.0d-3/akf(2)
      atrpl2 = log10(atrp2+1.0d-30)
      adlp2  = 1.d0/(atrn8-0.14d0*(atrpl2+atrc8))
      atra2  = (atrpl2+atrc8)*adlp2
      atrb2  = 1.d0/(1.d0 + atra2*atra2)
      atrf2  = 0.8d0**atrb2
      adpi2  = 1.d0/(1.d0 + atrp2)

      akf(2) = akf(2)*atrp2*adpi2*atrf2
!----------------------------------------------------------
!  the reaction rate constant of #03

      atrp3  = akf(5)*ax(2)*1.0d-3/akf(4)
      atrpl3 = log10(atrp3+1.0d-30)
      adlp3  = 1.d0/(atrn7-0.14d0*(atrpl3+atrc7))
      atra3  = (atrpl3+atrc7)*adlp3
      atrb3  = 1.d0/(1.d0 + atra3*atra3)
      atrf3  = 0.7d0**atrb3
      adpi3  = 1.d0/(1.d0 + atrp3)

      akf(3) = akf(4)*atrp3*adpi3*atrf3

!----------------------------------------------------------
!  the reaction rate constant of #04

      atrp4  = akf(7)*am04*1.0d-3/akf(6)
      atrp4r = 1.0d0/atrp4
      atrpl4 = log10(atrp4+1.0d-30)
      adlp4  = 1.d0/(atrn7-0.14d0*(atrpl4+atrc7))
      atra4  = (atrpl4+atrc7)*adlp4
      atrb4  = 1.d0/(1.d0 + atra4*atra4)
      atrf4  = 0.7d0**atrb4
      adpi4  = 1.d0/(1.d0 + atrp4)
      atrpr4(1  ) = akf(7)*1.5d-3*dmlr(1  )/akf(6)
      atrpr4(2  ) = 0.d0
      atrpr4(3:5) = akf(7)*1.0d-3*dmlr(3:5)/akf(6)
      atrpr4(6  ) = 0.d0
      atrpr4(7:9) = akf(7)*1.0d-3*dmlr(7:9)/akf(6)

      akf(4) = akf(6)*atrp4*adpi4*atrf4

!----------------------------------------------------------
!  the reaction rate constant of #05

      atrp5  = akf(9)*am05*1.0d-3/akf(8)
      atrp5r = 1.0d0/atrp5
      adpi5  = 1.d0/(1.d0 + atrp5)
      atrpr5(1:5) = akf(9)*1.0d-3*dmlr(1:5)/akf(8)
      atrpr5(6  ) = akf(9)*9.0d-3*dmlr(6  )/akf(8)
      atrpr5(7:8) = akf(9)*1.0d-3*dmlr(7:8)/akf(8)
      atrpr5(9  ) = akf(9)*1.5d-3*dmlr(9  )/akf(8)

      akf(5) = akf(8)*atrp5*adpi5

!-----------------------------------------------------------
!  the reaction rate constant of #06 & #08 & #12

!      akf06 = akf06a + akf06b
!      akf08 = akf08a + akf08b
!      akf12 = akf12a + akf12b
!------------------------------------------------------------

!     re-contaning akf
      akf(6:ler) = akf(10:lnr)


! [kb]
! backward reaction rate constant

! From the analysis in Intel Adviser,
! using user function in the loop takes longer time to calculate.
! In order to improve calculation performance,
! I refrain from using user function.
!
! Demerit is that the code is not beautiful without using user function


      itm = 1 + int(a05 + sign(a05,(atmp - dplt)))

      atmpl = log(atmp)

      ! make s/R[-], h/RT[-]
      as(1:8) =       atmpl* dplh(1,1:8,itm)                                 &
      &             + atmp *(dplh(2,1:8,itm)                                 &
      &             + atmp *(dplh(3,1:8,itm)*a12                             &
      &             + atmp *(dplh(4,1:8,itm)*a13                             &
      &             + atmp *(dplh(5,1:8,itm)*a14))))                         &
      &             +        dplh(7,1:8,itm)

      acon = atmpr
!     acon = arur*atmpr
!     acon = 1.0d0/(dru*atmp)

      ah(1:8) =    (atmp*(dplh(1,1:8,itm)                                    &
      &           + atmp*(dplh(2,1:8,itm)*a12                                &
      &           + atmp*(dplh(3,1:8,itm)*a13                                &
      &           + atmp*(dplh(4,1:8,itm)*a14                                &
      &           + atmp*(dplh(5,1:8,itm)*a15)))))                           &
      &           +       dplh(6,1:8,itm))                                   &
      &           * acon

      dcpo = sum(         (dplh(1,1:8,itm)                                 &
      &            + atmp*(dplh(2,1:8,itm)                               &
      &            + atmp*(dplh(3,1:8,itm)                               &
      &            + atmp*(dplh(4,1:8,itm)                               &
      &            + atmp*(dplh(5,1:8,itm))))))                          &
      &            * dmlru(1:8) * adns(1:8) ) * ar

      apt = arp*atmp
      aptr = 1.d0 / apt

      ahh(1 ) = ah(5) + ah(4) - ah(3) - ah(2)
      ahh(2 ) = ah(7) - ah(3) - ah(2)
      ahh(3 ) = ah(7) - ah(3) - ah(2)
      ahh(4 ) = ah(7) - ah(3) - ah(2)
      ahh(5 ) = ah(5) + ah(5) - ah(8)
      ahh(6 ) = ah(6) + ah(7) - ah(5) - ah(8)
      ahh(7 ) = ah(6) + ah(2) - ah(5) - ah(7)
      ahh(8 ) = ah(8) + ah(2) - ah(7) - ah(7)
      ahh(9 ) = ah(3) + ah(5) - ah(6)
      ahh(10) = ah(5) + ah(3) - ah(6)
      ahh(11) = ah(6) + ah(4) - ah(5) - ah(5)
      ahh(12) = ah(3) + ah(5) - ah(4) - ah(1)
      ahh(13) = ah(6) + ah(3) - ah(1) - ah(5)
      ahh(14) = ah(5) + ah(5) - ah(3) - ah(7)
      ahh(15) = ah(6) + ah(4) - ah(3) - ah(7)
      ahh(16) = ah(1) + ah(2) - ah(3) - ah(7)
      ahh(17) = ah(5) + ah(2) - ah(4) - ah(7)
      ahh(18) = ah(7) + ah(1) - ah(8) - ah(3)
      ahh(19) = ah(6) + ah(5) - ah(8) - ah(3)
      ahh(20) = ah(5) + ah(7) - ah(8) - ah(4)
      ahh(21) = ah(3) + ah(3) - ah(1)
      ahh(22) = ah(3) + ah(3) - ah(1)
      ahh(23) = ah(3) + ah(3) - ah(1)
      ahh(24) = ah(2) - ah(4) - ah(4)
      ahh(25) = ah(5) - ah(4) - ah(3)

      ag(1 ) = as(5)+as(4)-as(3)-as(2) - ahh(1 )
      ag(2 ) = as(7)-as(3)-as(2)       - ahh(2 )
      ag(3 ) = as(7)-as(3)-as(2)       - ahh(3 )
      ag(4 ) = as(7)-as(3)-as(2)       - ahh(4 )
      ag(5 ) = as(5)+as(5)-as(8)       - ahh(5 )
      ag(6 ) = as(6)+as(7)-as(5)-as(8) - ahh(6 )
      ag(7 ) = as(6)+as(2)-as(5)-as(7) - ahh(7 )
      ag(8 ) = as(8)+as(2)-as(7)-as(7) - ahh(8 )
      ag(9 ) = as(3)+as(5)-as(6)       - ahh(9 )
      ag(10) = as(5)+as(3)-as(6)       - ahh(10)
      ag(11) = as(6)+as(4)-as(5)-as(5) - ahh(11)
      ag(12) = as(3)+as(5)-as(4)-as(1) - ahh(12)
      ag(13) = as(6)+as(3)-as(1)-as(5) - ahh(13)
      ag(14) = as(5)+as(5)-as(3)-as(7) - ahh(14)
      ag(15) = as(6)+as(4)-as(3)-as(7) - ahh(15)
      ag(16) = as(1)+as(2)-as(3)-as(7) - ahh(16)
      ag(17) = as(5)+as(2)-as(4)-as(7) - ahh(17)
      ag(18) = as(7)+as(1)-as(8)-as(3) - ahh(18)
      ag(19) = as(6)+as(5)-as(8)-as(3) - ahh(19)
      ag(20) = as(5)+as(7)-as(8)-as(4) - ahh(20)
      ag(21) = as(3)+as(3)-as(1)       - ahh(21)
      ag(22) = as(3)+as(3)-as(1)       - ahh(22)
      ag(23) = as(3)+as(3)-as(1)       - ahh(23)
      ag(24) = as(2)-as(4)-as(4)       - ahh(24)
      ag(25) = as(5)-as(4)-as(3)       - ahh(25)

      akb(1    ) = akf(1    )/min(akmx,max(akmn,exp(max(aemn,min(ag(1    ),aemx)))     ))  ! 01
      akb(2 :4 ) = akf(2 :4 )/min(akmx,max(akmn,exp(max(aemn,min(ag(2 :4 ),aemx)))*apt ))  ! 02~04
      akb(5    ) = akf(5    )/min(akmx,max(akmn,exp(max(aemn,min(ag(5    ),aemx)))*aptr))  ! 05
      akb(6 :7 ) = akf(6 :7 )/min(akmx,max(akmn,exp(max(aemn,min(ag(6    ),aemx)))     ))  ! 06a,b
      akb(8    ) = akf(8    )/min(akmx,max(akmn,exp(max(aemn,min(ag(7    ),aemx)))     ))  ! 07
      akb(9 :10) = akf(9:10 )/min(akmx,max(akmn,exp(max(aemn,min(ag(8    ),aemx)))     ))  ! 08a,b
      akb(11:12) = akf(11:12)/min(akmx,max(akmn,exp(max(aemn,min(ag(9 :10),aemx)))*aptr))  ! 09~10
      akb(13   ) = akf(13   )/min(akmx,max(akmn,exp(max(aemn,min(ag(11   ),aemx)))     ))  ! 11
      akb(14:15) = akf(14:15)/min(akmx,max(akmn,exp(max(aemn,min(ag(12   ),aemx)))     ))  ! 12a,b
      akb(16:23) = akf(16:23)/min(akmx,max(akmn,exp(max(aemn,min(ag(13:20),aemx)))     ))  ! 13~20
      akb(24:26) = akf(24:26)/min(akmx,max(akmn,exp(max(aemn,min(ag(21:23),aemx)))*aptr))  ! 21~23
      akb(27:28) = akf(27:28)/min(akmx,max(akmn,exp(max(aemn,min(ag(24:25),aemx)))*apt ))  ! 24~25

!-----------------------------------------------------------------------
! adjust unit

      akf(1:ler) = akf(1:ler) * autf(1:ler)
      akb(1:ler) = akb(1:ler) * autb(1:ler)

!-----------------------------------------------------------------------
!   make kf[ ][ ] : reaction rate constant * mole density

      akfc(1 ) = akf(1 ) * ax(3) * ax(2)
      akfc(2 ) = akf(2 ) * ax(3) * ax(2)
      akfc(3 ) = akf(3 ) * ax(3) * ax(2)
      akfc(4 ) = akf(4 ) * ax(3) * ax(2)
      akfc(5 ) = akf(5 ) * ax(8)
      akfc(6 ) = akf(6 ) * ax(5) * ax(8)
      akfc(7 ) = akf(7 ) * ax(5) * ax(8)
      akfc(8 ) = akfc(6) + akfc(7)                 !06a+06b
      akfc(9 ) = akf(8 ) * ax(5) * ax(7)
      akfc(10) = akf(9 ) * ax(7) * ax(7)
      akfc(11) = akf(10) * ax(7) * ax(7)
      akfc(12) = akfc(10) + akfc(11)               !08a+08b
      akfc(13) = akf(11) * ax(6)          * am09
      akfc(14) = akf(12) * ax(6)          * ax(6)
      akfc(15) = akf(13) * ax(5) * ax(5)
      akfc(16) = akf(14) * ax(4) * ax(1)
      akfc(17) = akf(15) * ax(4) * ax(1)
      akfc(18) = akfc(16) + akfc(17)               !12a+12b
      akfc(19) = akf(16) * ax(1) * ax(5)
      akfc(20) = akf(17) * ax(3) * ax(7)
      akfc(21) = akf(18) * ax(3) * ax(7)
      akfc(22) = akf(19) * ax(3) * ax(7)
      akfc(23) = akf(20) * ax(4) * ax(7)
      akfc(24) = akf(21) * ax(8) * ax(3)
      akfc(25) = akf(22) * ax(8) * ax(3)
      akfc(26) = akf(23) * ax(8) * ax(4)
      akfc(27) = akf(24) * ax(1)         * am21
      akfc(28) = akf(25) * ax(1)         * ax(1)
      akfc(29) = akf(26) * ax(1)         * am23
      akfc(30) = akf(27) * ax(4) * ax(4) * am24
      akfc(31) = akf(28) * ax(4) * ax(3) * am25


!     make dkf/dT*[ ][ ]
      ! adf :: 1/kf*dkf/dT
      adf(1    ) = (ant(1     ) - ae(1     )) * atmpr
      adf(2    ) = (ant(2     ) - ae(2     )) * atmpr
      adf(3    ) = (ant(4     ) - ae(4     )) * atmpr
      adf(4    ) = (ant(6     ) - ae(6     )) * atmpr
      adf(5    ) = (ant(8     ) - ae(8     )) * atmpr
      adf(6:ler) = (ant(10:lnr) - ae(10:lnr)) * atmpr

      add02  = (ant(3) - ant(2) +  ae(3)- ae(2)) * atmpr
      add03  = (ant(5) - ant(4) +  ae(5)- ae(4)) * atmpr
      add04  = (ant(7) - ant(6) +  ae(7)- ae(6)) * atmpr
      add05  = (ae(9) - ae(8)) * atmpr

      ! about pressure dependent reaction
      adfg2  = 1.d0/(1.d0 + atrb2)/(1.d0 + atrb2)
      adfg3  = 1.d0/(1.d0 + atrb3)/(1.d0 + atrb3)
      adfg4  = 1.d0/(1.d0 + atrb4)/(1.d0 + atrb4)

      atrr2  = -2.d0*atrf8*atra2*adfg2*atrn8*adlp2*adlp2
      atrr3  = -2.d0*atrf7*atra3*adfg3*atrn7*adlp3*adlp3
      atrr4  = -2.d0*atrf7*atra4*adfg4*atrn7*adlp4*adlp4

      adf(2)  = adf(2) + (adpi2 + atrr2)*add02
      adf(3)  = adf(3) + (adpi3 + atrr3)*add03
      adf(4)  = adf(4) + (adpi4 + atrr4)*add04
      adf(5)  = adf(5) + adpi5*add05

      ! akft :: adf*kf[][] = dkf/dT*[][]
      akft(1 :5 ) = adf(1 :5 ) * akfc(1 :5 )
      akft(6    ) = adf(6    ) * akfc(6    ) + adf(7 ) * akfc(7 )
      akft(7    ) = adf(8    ) * akfc(9    )
      akft(8    ) = adf(9    ) * akfc(10   ) + adf(10) * akfc(11)
      akft(9:11 ) = adf(11:13) * akfc(13:15)
      akft(12   ) = adf(14   ) * akfc(16   ) + adf(15) * akfc(17)
      akft(13:25) = adf(16:28) * akfc(19:31)


!     make dkf/dri*[][]

      akfr2      = (adpi2 + atrr2/(atrp2  + aemn)) * (adnsr(6 ) + aemn) * akbc(2)
      akfr3      = (adpi3 + atrr3/(atrp3  + aemn)) * (adnsr(2 ) + aemn) * akbc(3)
      akfr4(1  ) = (adpi4 + atrr4*(atrp4r + aemn)) * (atrp4r    + aemn) * akbc(4) * atrpr4(1  )
      akfr4(2  ) = 0.d0
      akfr4(3:5) = (adpi4 + atrr4*(atrp4r + aemn)) * (atrp4r    + aemn) * akbc(4) * atrpr4(3:5)
      akfr4(6  ) = 0.d0
      akfr4(7:9) = (adpi4 + atrr4*(atrp4r + aemn)) * (atrp4r    + aemn) * akbc(4) * atrpr4(7:9)
      akfr5(1:9) = adpi5 * (atrp5r + aemn) * akbc(5) * atrpr5(1:9)

!     make kb[ ][ ]

      akbc(1 ) = akb(1 ) * ax(5) * ax(4)
      akbc(2 ) = akb(2 ) * ax(7)
      akbc(3 ) = akb(3 ) * ax(7)
      akbc(4 ) = akb(4 ) * ax(7)
      akbc(5 ) = akb(5 ) * ax(5) * ax(5)
      akbc(6 ) = akb(6 ) * ax(6) * ax(7)
      akbc(7 ) = akb(7 ) * ax(6) * ax(7)
      akbc(8 ) = akbc(6) + akbc(7)
      akbc(9 ) = akb(8 ) * ax(6) * ax(2)
      akbc(10) = akb(9 ) * ax(2) * ax(8)
      akbc(11) = akb(10) * ax(2) * ax(8)
      akbc(12) = akbc(10) + akbc(11)
      akbc(13) = akb(11) * ax(3) * ax(5) * am09
      akbc(14) = akb(12) * ax(5) * ax(3) * ax(6)
      akbc(15) = akb(13) * ax(6) * ax(4)
      akbc(16) = akb(14) * ax(3) * ax(5)
      akbc(17) = akb(15) * ax(3) * ax(5)
      akbc(18) = akbc(16) + akbc(17)
      akbc(19) = akb(16) * ax(6) * ax(3)
      akbc(20) = akb(17) * ax(5) * ax(5)
      akbc(21) = akb(18) * ax(6) * ax(4)
      akbc(22) = akb(19) * ax(1) * ax(2)
      akbc(23) = akb(20) * ax(5) * ax(2)
      akbc(24) = akb(21) * ax(7) * ax(1)
      akbc(25) = akb(22) * ax(6) * ax(5)
      akbc(26) = akb(23) * ax(5) * ax(7)
      akbc(27) = akb(24) * ax(3) * ax(3) * am21
      akbc(28) = akb(25) * ax(3) * ax(3) * ax(1)
      akbc(29) = akb(26) * ax(3) * ax(3) * am23
      akbc(30) = akb(27) * ax(2)         * am24
      akbc(31) = akb(28) * ax(5)         * am25


!     make dkb/dT*[ ][ ]

      adb(1 ) = adf(1 ) - (ahh(1 )        )*atmpr
      adb(2 ) = adf(2 ) - (ahh(2 ) + 1.0d0)*atmpr
      adb(3 ) = adf(3 ) - (ahh(3 ) + 1.0d0)*atmpr
      adb(4 ) = adf(4 ) - (ahh(4 ) + 1.0d0)*atmpr
      adb(5 ) = adf(5 ) - (ahh(5 ) - 1.0d0)*atmpr
      adb(6 ) = adf(6 ) - (ahh(6 )        )*atmpr   !06a
      adb(7 ) = adf(7 ) - (ahh(6 )        )*atmpr   !06b
      adb(8 ) = adf(8 ) - (ahh(7 )        )*atmpr
      adb(9 ) = adf(9 ) - (ahh(8 )        )*atmpr   !08a
      adb(10) = adf(10) - (ahh(8 )        )*atmpr   !08b
      adb(11) = adf(11) - (ahh(9 ) - 1.0d0)*atmpr
      adb(12) = adf(12) - (ahh(10) - 1.0d0)*atmpr
      adb(13) = adf(13) - (ahh(11)        )*atmpr
      adb(14) = adf(14) - (ahh(12)        )*atmpr   !12a
      adb(15) = adf(15) - (ahh(12)        )*atmpr   !12b
      adb(16) = adf(16) - (ahh(13)        )*atmpr
      adb(17) = adf(17) - (ahh(14)        )*atmpr
      adb(18) = adf(18) - (ahh(15)        )*atmpr
      adb(19) = adf(19) - (ahh(16)        )*atmpr
      adb(20) = adf(20) - (ahh(17)        )*atmpr
      adb(21) = adf(21) - (ahh(18)        )*atmpr
      adb(22) = adf(22) - (ahh(19)        )*atmpr
      adb(23) = adf(23) - (ahh(20)        )*atmpr
      adb(24) = adf(24) - (ahh(21) - 1.0d0)*atmpr
      adb(25) = adf(25) - (ahh(22) - 1.0d0)*atmpr
      adb(26) = adf(26) - (ahh(23) - 1.0d0)*atmpr
      adb(27) = adf(27) - (ahh(24) + 1.0d0)*atmpr
      adb(28) = adf(28) - (ahh(25) + 1.0d0)*atmpr

      akbt(1 :5 ) = adb(1 :5 ) * akbc(1 :5 )
      akbt(6    ) = adb(6    ) * akbc(6    ) + adb(7 ) * akbc(7 )
      akbt(7    ) = adb(8    ) * akbc(9    )
      akbt(8    ) = adb(9    ) * akbc(10   ) + adb(10) * akbc(11)
      akbt(9:11 ) = adb(11:13) * akbc(13:15)
      akbt(12   ) = adb(14   ) * akbc(16   ) + adb(15) * akbc(17)
      akbt(13:25) = adb(16:28) * akbc(19:31)


!     make dkb/dri*[][]

      akbr2      = (adpi2 + atrr2/(atrp2  + aemn)) * (adnsr(6 ) + aemn) * akbc(2)
      akbr3      = (adpi3 + atrr3/(atrp3  + aemn)) * (adnsr(2 ) + aemn) * akbc(3)
      akbr4(1  ) = (adpi4 + atrr4*(atrp4r + aemn)) * (atrp4r    + aemn) * akbc(4) * atrpr4(1  )
      akbr4(2  ) = 0.d0
      akbr4(3:5) = (adpi4 + atrr4*(atrp4r + aemn)) * (atrp4r    + aemn) * akbc(4) * atrpr4(3:5)
      akbr4(6  ) = 0.d0
      akbr4(7:9) = (adpi4 + atrr4*(atrp4r + aemn)) * (atrp4r    + aemn) * akbc(4) * atrpr4(7:9)
      akbr5(1:9) = adpi5 * (atrp5r + aemn) * akbc(5) * atrpr5(1:9)


!     re-contaning akfc and akbc

      akfc(6 :7 ) = akfc(8 :9 )
      akfc(8 :11) = akfc(12:15)
      akfc(12:25) = akfc(18:31)

      akbc(6 :7 ) = akbc(8 :9 )
      akbc(8 :11) = akbc(12:15)
      akbc(12:25) = akbc(18:31)


!     reaction rate of each elementary reaction
!     kf[ }[ ] - kb[ ][ ]

      ars(1:len) = akfc(1:len) - akbc(1:len)


!    reaction rate of each species

      aw(1) = dml(1) * (-ars(12) -ars(13) +ars(16) +ars(18) -ars(21)          &
      &                 -ars(22) -ars(23))
      aw(2) = dml(2) * (-ars(1 ) -ars(2 ) -ars(3 ) -ars(4 ) +ars(7 )          &
      &                 +ars(8 ) +ars(16) +ars(17) +ars(24))
      aw(3) = dml(3) * (-ars(1 ) -ars(2 ) -ars(3 ) -ars(4 ) +ars(9 )          &
      &                 +ars(10) +ars(12) +ars(13) -ars(14) -ars(15)          &
      &                 -ars(16) -ars(18) +ars(19) +ars(21) +ars(21)          &
      &                 +ars(22) +ars(22) +ars(23) +ars(23) -ars(25))
      aw(4) = dml(4) * (+ars(1 ) +ars(11) -ars(12) +ars(15) -ars(17)          &
      &                 -ars(20) -ars(24) -ars(24) -ars(25))
      aw(5) = dml(5) * (+ars(1 ) +ars(5 ) +ars(5 ) -ars(6 ) -ars(7 )          &
      &                 +ars(9 ) +ars(10) -ars(11) -ars(11) +ars(12)          &
      &                 -ars(13) +ars(14) +ars(14) +ars(17) +ars(19)          &
      &                 +ars(20) +ars(25))
      aw(6) = dml(6) * (+ars(6 ) +ars(7 ) -ars(9 ) -ars(10) +ars(11)          &
      &                 +ars(13) +ars(15) +ars(19))
      aw(7) = dml(7) * (+ars(2 ) +ars(3 ) +ars(4 ) +ars(6 ) -ars(7 )          &
      &                 -ars(8 ) -ars(8 ) -ars(14) -ars(15) -ars(16)          &
      &                 -ars(17) +ars(18) +ars(20))
      aw(8) = dml(8) * (-ars(5 ) -ars(6 ) +ars(8 ) -ars(18) -ars(19)          &
      &                 -ars(20))


!     source term of H2

!      drst(j,k) = aw(1)

      arhs(1:8) = delt*aw(1:8)

!-----------------------------------------------------------------------
! make jacobian
!-----------------------------------------------------------------------

! dT/drj

    druo = sum(adns(1:lsp)*dmlru(1:lsp))*ar
    aruo = 1/druo
    awk1 = dcpo/druo
    awk2 = 1.d0 / (1.d0 - awk1)
    awk3 = awk1 * dru * atmp
    dppd(1:8) = (ahti(1:8) - awk3*dmlr(1:8)) * awk2
    arr  = aruo*arho
    adtd(1:8) = (dppd(1:8)-dru*atmp*dmlr(1:8))*arr

! elementary reaction concerning with [M], kf[][][M]/[M]

    arf09 = akf(11) * ax(6)
    arf10 = akf(12) * ax(6)
    arf21 = akf(24) * ax(1)
    arf22 = akf(25) * ax(1)
    arf23 = akf(26) * ax(1)
    arf24 = akf(27) * ax(4) * ax(4)
    arf25 = akf(28) * ax(3) * ax(4)

    arb09 = akb(11) * ax(3) * ax(5)
    arb10 = akb(12) * ax(3) * ax(5)
    arb21 = akb(24) * ax(3) * ax(3)
    arb22 = akb(25) * ax(3) * ax(3)
    arb23 = akb(26) * ax(3) * ax(3)
    arb24 = akb(27) * ax(2)
    arb25 = akb(28) * ax(5)

    acn = 0.5d0*delt


!-----------
!   dw1   - H2
!-----------

    ! H2 coms or prod rate of third-body reactions in the reactions including H2 (only No.7 reaction)
    admd(1) =                arf22
    admd(2) =                        arb23
    admd(3) = arb21
    admd(4) = arb21
    admd(5) = arb21
    admd(6) = arb21*14.4d0
    admd(7) = arb21
    admd(8) = arb21
    admd(9) =                        arb23

    ! each comsumption rate in the reactions including H2
    akc(1) = akfc(12) + akfc(13) + akbc(16) + akbc(18) + akfc(21) + akfc(22)  &
    &      + akfc(23)
    akc(2) = 0.d0
    akc(3) = akbc(12) + akbc(13) + akfc(16) + akfc(18) + akbc(21) + akbc(21)  &
    &      + akbc(22) + akbc(22) + akbc(23) + akbc(23)
    akc(4) = 0.d0
    akc(5) = akbc(12)
    akc(6) = akbc(13)
    akc(7) = akfc(16)
    akc(8) = akfc(18)
    akc(9) = 0.d0

    ! H2 comsumption and production (asktw"c" and "p")
    asktwc = akft(12) + akft(13) + akbt(16) + akbt(18) + akft(21)             &
    &      + akft(22) + akft(23)
    asktwp = akbt(12) + akbt(13) + akft(16) + akft(18) + akbt(21)             &
    &      + akbt(22) + akbt(23)

    ! tau^-1 = max(dwj/drj,dwj/drk)  -inverse of characteristic time
    aw(1    ) = admd(1    )*dmlr(1    ) + akc(1    )*adnsr(1    ) + adtd(1    )*asktwc
    aw(2:lsp) = admd(2:lsp)*dmlr(2:lsp) + akc(2:lsp)*adnsr(2:lsp) + adtd(2:lsp)*asktwp

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(1) = 1.0d0 + acn*dml(1)*atmaxi

!-----------
!   dw2   - O2
!-----------

    admd(1) = arf24*2.5d0  ! [M] = H2
    admd(2) = arb24        ! [M] = O2
    admd(3) = arf24        ! [M] = H
    admd(4) = arf24        ! [M] = O
    admd(5) = arf24        ! [M] = OH
    admd(6) = arf24*12.0d0 ! [M] = H2O
    admd(7) = arf24        ! [M] = HO2
    admd(8) = arf24        ! [M] = H2O2
    admd(9) = arf24        ! [M] = N2

    akc(1) = 0.d0
    akc(2) = akfc(1 ) + akfc(2 ) + akfc(3 ) + akfc(4 ) + akbc(7 )           &
    &      + akbc(8 ) + akbc(16) + akbc(17) + akbc(24)
    akc(3) = akfc(16)
    akc(4) = akbc(1 ) + akfc(17) + akfc(24) + akfc(24)
    akc(5) = akbc(1 ) + akfc(7 )
    akc(6) = 0.d0
    akc(7) = akbc(2 ) + akbc(3 ) + akbc(4 ) + akfc(7 ) + akfc(8 )           &
    &      + akfc(8 ) + akfc(16) + akfc(17)
    akc(8) = 0.d0
    akc(9) = 0.d0

    asktwc = akft(1 ) + akft(2 ) + akft(3 ) + akft(4 ) + akbt(7 )           &
    &      + akbt(8 ) + akbt(16) + akbt(17) + akbt(24)
    asktwp = akbt(1 ) + akbt(2 ) + akbt(3 ) + akbt(4 ) + akft(7 )           &
    &      + akft(8 ) + akft(16) + akft(17) + akft(24)

    askrwc    = akfr3
    askrwp(1) = akbr4(1)
    askrwp(2) = 0.d0
    askrwp(3) = akbr4(3)
    askrwp(4) = akbr4(4)
    askrwp(5) = akbr4(5)
    askrwp(6) = akbr2
    askrwp(7) = akbr4(7)
    askrwp(8) = akbr4(8)
    askrwp(9) = akbr4(9)

    aw(1    ) = admd(1    )*dmlr(1    ) + akc(1    )*adnsr(1    ) + adtd(1    )*asktwp + askrwp(1    )
    aw(2    ) = admd(2    )*dmlr(2    ) + akc(2    )*adnsr(2    ) + adtd(2    )*asktwc + askrwc
    aw(3:lsp) = admd(3:lsp)*dmlr(3:lsp) + akc(3:lsp)*adnsr(3:lsp) + adtd(3:lsp)*asktwp + askrwp(3:lsp)

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(2) = 1.0d0 + acn*dml(2)*atmaxi


!----------
!   dw3  -  H
!-----------

    admd(1) = arf09*3.0d0 + arf22 + arf22 + arb25*2.5d0
    admd(2) = arf09*1.5d0 + arf23 + arf23 + arb25
    admd(3) = arb09 + arb21 + arb21 + arf25
    admd(4) = arf09 + arf21 + arf21 + arb25
    admd(5) = arf09 + arf21 + arf21 + arb25
    admd(6) = arf10 + arf21*14.4d0 + arf21*14.4d0 + arb25*12.d0
    admd(7) = arf09 + arf21 + arf21 + arb25
    admd(8) = arf09 + arf21 + arf21 + arb25
    admd(9) = arf09*2.0d0 + arf23 + arf23 + arb25

    akc(1) = akfc(12) + akfc(13) + akbc(16) + akbc(18) + akfc(21) +  akfc(21) &
    &      + akfc(22) + akfc(22) + akfc(23) + akfc(23)
    akc(2) = akbc(16)
    akc(3) = akfc(1 ) + akfc(2 ) + akfc(3 ) + akfc(4 ) + akbc(9 ) + akbc(10)  &
    &      + akbc(12) + akbc(13) + akfc(14) + akfc(15) + akfc(16) + akfc(18)  &
    &      + akfc(19) + 4.d0 * akbc(21) + 4.d0 * akbc(22) + 4.d0 * akbc(23)   &
    &      + akfc(25)
    akc(4) = akbc(1 ) + akfc(12) + akbc(15)
    akc(5) = akbc(1 ) + akfc(13) + akbc(14) + akbc(14) + akbc(19) + akbc(25)
    akc(6) = akfc(9 ) + akfc(10) + akbc(15) + akbc(19)
    akc(7) = akbc(2 ) + akbc(3 ) + akbc(4 ) + akbc(18)
    akc(8) = 0.d0
    akc(9) = 0.d0

    asktwc = akft(1 ) + akft(2 ) + akft(3 ) + akft(4 ) + akbt(9 ) + akbt(10)  &
    &      + akbt(12) + akbt(13) + akft(14) + akft(15) + akft(16) + akft(18)  &
    &      + akft(19) + akbt(21) + akbt(21) + akbt(22) + akbt(22) + akbt(23)  &
    &      + akft(25)
    asktwp = akbt(1 ) + akbt(2 ) + akbt(3 ) + akbt(4 ) + akft(9 ) + akft(10)  &
    &      + akft(12) + akft(13) + akbt(14) + akbt(15) + akbt(16) + akbt(18)  &
    &      + akbt(19) + akft(21) + akft(21) + akft(22) + akft(22) + akft(23)  &
    &      + akbt(25)

    askrwc    = akfr4(3)
    askrwp(1) = akbr4(1)
    askrwp(2) = akbr3
    askrwp(3) = 0.d0
    askrwp(4) = akbr4(4)
    askrwp(5) = akbr4(5)
    askrwp(6) = akbr2
    askrwp(7) = akbr4(7)
    askrwp(8) = akbr4(8)
    askrwp(9) = akbr4(9)

    aw(1:2  ) = admd(1:2  )*dmlr(1:2  ) + akc(1:2  )*adnsr(1:2  ) + adtd(1:2  )*asktwp + askrwp(1:2  )
    aw(3    ) = admd(3    )*dmlr(3    ) + akc(3    )*adnsr(3    ) + adtd(3    )*asktwp + askrwc
    aw(4:lsp) = admd(4:lsp)*dmlr(4:lsp) + akc(4:lsp)*adnsr(4:lsp) + adtd(4:lsp)*asktwc + askrwp(4:lsp)

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(3) = 1.0d0 + acn*dml(3)*atmaxi

!-----------
!   dw4  -  O
!-----------

    admd(1) = arb24*2.5d0 + arb24*2.5d0 + arb25*2.5d0
    admd(2) = arb24       + arb24       + arb25
    admd(3) = arb24       + arb24       + arb25
    admd(4) = arf24       + arf24       + arf25
    admd(5) = arb24       + arb24       + arb25
    admd(6) = arb24*12.d0 + arb24*12.d0 + arb25*12.d0
    admd(7) = arb24       + arb24       + arb25
    admd(8) = arb24       + arb24       + arb25
    admd(9) = arb24       + arb24       + arb25

    akc(1) = 0.d0
    akc(2) = akfc(1 ) + akbc(17) + akbc(24) + akbc(24)
    akc(3) = akfc(1 ) + akbc(12) + akfc(15)
    akc(4) = akbc(1 ) + akbc(11) + akfc(12) + akbc(15) + akfc(17) + akfc(20)  &
    &      + 4.d0 * akfc(24) + akfc(25)
    akc(5) = akfc(11) + akfc(11) + akbc(12) + akbc(17) + akbc(20) + akbc(25)
    akc(6) = 0.d0
    akc(7) = akfc(15) + akbc(20)
    akc(8) = 0.d0
    akc(9) = 0.d0

    asktwc = akbt(1 ) + akbt(11) + akft(12) + akbt(15) + akft(17) + akft(20)  &
    &      + akft(24) + akft(24) + akft(25)
    asktwp = akft(1 ) + akft(11) + akbt(12) + akft(15) + akbt(17) + akbt(20)  &
    &      + akbt(24) + akbt(24) + akbt(25)

    aw(1:3  ) = admd(1:3  )*dmlr(1:3  ) + akc(1:3  )*adnsr(1:3  ) + adtd(1:3  )*asktwp
    aw(4    ) = admd(4    )*dmlr(4    ) + akc(4    )*adnsr(4    ) + adtd(4    )*asktwp
    aw(5:lsp) = admd(5:lsp)*dmlr(5:lsp) + akc(5:lsp)*adnsr(5:lsp) + adtd(5:lsp)*asktwp

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(4) = 1.0d0 + acn*dml(4)*atmaxi

!-----------
!   dw5  -  OH
!-----------

    admd(1) = arf09*3.0d0         + arf25*2.5d0
    admd(2) = arf09*1.5d0         + arf25
    admd(3) = arf09               + arf25
    admd(4) = arf09               + arf25
    admd(5) = arb09               + arb25
    admd(6) =               arf10 + arf25*12.0d0
    admd(7) = arf09               + arf25
    admd(8) = arf09               + arf25
    admd(9) = arf09*2.0d0         + arf25

    akc(1) = akfc(12)
    akc(2) = akfc(1 ) + akbc(7 )
    akc(3) = akfc(1 ) + akbc(13) + akfc(14) + akfc(14) + akfc(19) + akfc(25)
    akc(4) = akbc(11) + akbc(11) + akfc(12) + akfc(17) + akfc(20) + akfc(25)
    akc(5) = akbc(1 ) + 4.d0 * akbc(5 )     + akfc(6 ) + akfc(7 ) + akbc(9 )  &
    &      + akbc(10) + 4.d0 * akfc(11)     + akbc(12) + akfc(13)             &
    &      + 4.d0 * akbc(14)     + akbc(17) + akbc(19) + akbc(20) + akbc(25)
    akc(6) = akbc(6 ) + akbc(7 ) + akfc(9 ) + akfc(10) + akbc(11) + akbc(11)  &
    &      + akbc(13)
    akc(7) = akbc(6 ) + akfc(14) + akfc(14) + akfc(17)
    akc(8) = akfc(5 ) + akfc(5 ) + akfc(19) + akfc(20)
    akc(9) = 0.0d0

    asktwc = akbt(1 ) + akbt(5 ) + akbt(5 ) + akft(6 ) + akft(7 ) + akbt(9 )  &
    &      + akbt(10) + akft(11) + akft(11) + akbt(12) + akft(13) + akbt(14)  &
    &      + akbt(14) + akbt(17) + akbt(19) + akbt(20) + akbt(25)
    asktwp = akft(1 ) + akft(5 ) + akbt(5 ) + akbt(6 ) + akbt(7 ) + akft(9 )  &
    &      + akft(10) + akbt(11) + akbt(11) + akft(12) + akbt(13) + akft(14)  &
    &      + akft(14) + akft(17) + akft(19) + akft(20) + akft(25)

    askrwc    = 2.d0 * akbr5(5)
    askrwp(1) = 2.d0 * akfr5(1)
    askrwp(2) = 2.d0 * akfr5(2)
    askrwp(3) = 2.d0 * akfr5(3)
    askrwp(4) = 2.d0 * akfr5(4)
    askrwp(5) = 0.d0
    askrwp(6) = 2.d0 * akfr5(6)
    askrwp(7) = 2.d0 * akfr5(7)
    askrwp(8) = 2.d0 * akfr5(8)
    askrwp(9) = 2.d0 * akfr5(9)

    aw(1:4  ) = admd(1:4  )*dmlr(1:4  ) + akc(1:4  )*adnsr(1:4  ) + adtd(1:4  )*asktwp + askrwp(1:4)
    aw(5    ) = admd(5    )*dmlr(5    ) + akc(5    )*adnsr(5    ) + adtd(5    )*asktwp + askrwc
    aw(6:lsp) = admd(6:lsp)*dmlr(6:lsp) + akc(6:lsp)*adnsr(6:lsp) + adtd(6:lsp)*asktwp + askrwp(6:lsp)

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(5) = 1.0d0 + acn*dml(5)*atmaxi

!-----------
!   dw6  -  H2O
!-----------

    admd(1) = arb09*3.0d0
    admd(2) = arb09*1.5d0
    admd(3) = arb09
    admd(4) = arb09
    admd(5) = arb09
    admd(6) =             + arf10
    admd(7) = arb09
    admd(8) = arb09
    admd(9) = arb09*2.0d0

    akc(1) = akfc(13)
    akc(2) = 0.d0
    akc(3) = akbc(9 ) + akbc(10) + akfc(15) + akfc(19)
    akc(4) = 0.d0
    akc(5) = akfc(6 ) + akfc(7 ) + akbc(9 ) + akbc(10) + akfc(11) + akfc(11)  &
    &      + akfc(13)
    akc(6) = akbc(6 ) + akbc(7 ) + akfc(9 ) + akfc(10) + akbc(11) + akbc(13)  &
    &      + akbc(15) + akbc(19)
    akc(7) = akfc(7 ) + akfc(15)
    akc(8) = akfc(6 ) + akfc(19)
    akc(9) = 0.0d0

    asktwc = akbt(6 ) + akbt(7 ) + akft(9 ) + akft(10) + akbt(11) + akbt(13)  &
    &      + akbt(15) + akbt(19)
    asktwp = akft(6 ) + akft(7 ) + akbt(9 ) + akbt(10) + akft(11) + akft(13)  &
    &      + akft(15) + akft(19)

    aw(1:5  ) = admd(1:5  )*dmlr(1:5  ) + akc(1:5  )*adnsr(1:5  ) + adtd(1:5  )*asktwp
    aw(6    ) = admd(6    )*dmlr(6    ) + akc(6    )*adnsr(6    ) + adtd(6    )*asktwp
    aw(7:lsp) = admd(7:lsp)*dmlr(7:lsp) + akc(7:lsp)*adnsr(7:lsp) + adtd(7:lsp)*asktwp

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(6) = 1.0d0 + acn*dml(6)*atmaxi

!-----------
!   dw7  -  HO2
!-----------

    admd(1) = 0.d0
    admd(2) = 0.d0
    admd(3) = 0.d0
    admd(4) = 0.d0
    admd(5) = 0.d0
    admd(6) = 0.d0
    admd(7) = 0.d0
    admd(8) = 0.d0
    admd(9) = 0.d0

    akc(1) = akbc(16)
    akc(2) = akfc(2 ) + akfc(3 ) + akfc(4 ) + akbc(7 ) + akbc(8 ) + akbc(8 )  &
    &      + akbc(16) + akbc(17)
    akc(3) = akfc(2 ) + akfc(3 ) + akfc(4 ) + akfc(18)
    akc(4) = akbc(15) + akfc(20)
    akc(5) = akfc(6 ) + akbc(14) + akbc(14) + akbc(17)
    akc(6) = akbc(7 ) + akbc(15)
    akc(7) = akbc(2 ) + akbc(3 ) + akbc(4 ) + akbc(6 ) + akfc(7 )             &
    &      + 4.d0 * akfc(8 )     + akfc(14) + akfc(15) + akfc(16) + akfc(17)  &
    &      + akbc(18) + akbc(20)
    akc(8) = akfc(6 ) + akbc(8 ) + akbc(8 ) + akfc(18) + akfc(20)
    akc(9) = 0.0d0

    asktwc = akbt(2 ) + akbt(3 ) + akbt(4 ) + akbt(6 ) + akft(7 ) + akft(8 )  &
    &      + akft(8 ) + akft(14) + akft(15) + akft(16) + akft(17) + akbt(18)  &
    &      + akbt(20)
    asktwp = akft(2 ) + akft(3 ) + akft(4 ) + akft(6 ) + akbt(7 ) + akbt(8 )  &
    &      + akbt(8 ) + akbt(14) + akbt(15) + akbt(16) + akbt(17) + akft(18)  &
    &      + akft(20)

    askrwc    = akbr4(7)
    askrwp(1) = akbr4(1)
    askrwp(2) = akbr2
    askrwp(3) = akbr4(3)
    askrwp(4) = akbr4(4)
    askrwp(5) = akbr4(5)
    askrwp(6) = akbr2
    askrwp(7) = 0.d0
    askrwp(8) = akbr4(8)
    askrwp(9) = akbr4(9)

    aw(1:6  ) = admd(1:6  )*dmlr(1:6  ) + akc(1:6  )*adnsr(1:6  ) + adtd(1:6  )*asktwp + askrwp(1:6  )
    aw(7    ) = admd(7    )*dmlr(7    ) + akc(7    )*adnsr(7    ) + adtd(7    )*asktwp + askrwc
    aw(8:lsp) = admd(8:lsp)*dmlr(8:lsp) + akc(8:lsp)*adnsr(8:lsp) + adtd(8:lsp)*asktwp + askrwp(8:lsp)

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(7) = 1.0d0 + acn*dml(7)*atmaxi

!----------------
!        dw8  -  H2O2
!----------------

    admd(1) = 0.d0
    admd(2) = 0.d0
    admd(3) = 0.d0
    admd(4) = 0.d0
    admd(5) = 0.d0
    admd(6) = 0.d0
    admd(7) = 0.d0
    admd(8) = 0.d0
    admd(9) = 0.d0

    akc(1) = akbc(18)
    akc(2) = 0.d0
    akc(3) = 0.d0
    akc(4) = 0.d0
    akc(5) = akbc(5 ) + akbc(5 ) + akbc(19) + akbc(20)
    akc(6) = akbc(6 ) + akbc(19)
    akc(7) = akbc(6 ) + akfc(8 ) + akfc(8 ) + akbc(18) + akbc(20)
    akc(8) = akfc(5 ) + akfc(6 ) + akbc(8 ) + akfc(18) + akfc(19) + akfc(20)
    akc(9) = 0.0d0

    asktwc = akft(5 ) + akft(6 ) + akbt(8 ) + akft(18) + akft(19) + akft(20)
    asktwp = akbt(5 ) + akbt(6 ) + akft(8 ) + akbt(18) + akbt(19) + akbt(20)

    askrwc    = akfr5(5)
    askrwp(1) = akbr5(1)
    askrwp(2) = akbr5(2)
    askrwp(3) = akbr5(3)
    askrwp(4) = akbr5(4)
    askrwp(5) = akbr5(5)
    askrwp(6) = akbr5(6)
    askrwp(7) = akbr5(7)
    askrwp(8) = 0.d0
    askrwp(9) = akbr5(9)

    aw(1:7  ) = admd(1:7  )*dmlr(1:7  ) + akc(1:7  )*adnsr(1:7  ) + adtd(1:7  )*asktwp + askrwp(1:7)
    aw(8    ) = admd(8    )*dmlr(8    ) + akc(8    )*adnsr(8    ) + adtd(8    )*asktwc + askrwc
    aw(9    ) = admd(9    )*dmlr(9    ) + akc(9    )*adnsr(9    ) + adtd(9    )*asktwp + askrwp(9  )

    atmaxi = max(aw(1),aw(2),aw(3),aw(4),aw(5),aw(6),aw(7),aw(8),aw(9))
    a1diag(8) = 1.0d0 + acn*dml(8)*atmaxi

!--------------------------------------------------
! [I + sita*Diag(delt/taui)]deltQ = RHSn]
! sita = 1/2 for Crank-Nicholson

    arhs(1:8) = arhs(1:8) / a1diag(1:8)

    adns(1:8) = adns(1:8) + arhs(1:8)

!---------------------------------------------------------------------
! update temperature

    afm(1:lsp) = adns(1:lsp)*dmlru(1:lsp)

    afmt = sum(afm(:))

    at = atmp

! --------------------
! newton-rapson method
! --------------------
    do i = 1, itr

      itm = 1 + int(a05 + sign(a05,(at - dplt)))

      ab(1) = sum(afm(1:lsp)*dplh(1,1:lsp,itm)) - afmt
      ab(2) = sum(afm(1:lsp)*dplh(2,1:lsp,itm)) * a12
      ab(3) = sum(afm(1:lsp)*dplh(3,1:lsp,itm)) * a13
      ab(4) = sum(afm(1:lsp)*dplh(4,1:lsp,itm)) * a14
      ab(5) = sum(afm(1:lsp)*dplh(5,1:lsp,itm)) * a15
      ab(6) = sum(afm(1:lsp)*dplh(6,1:lsp,itm)) + aeng

      afdd = 2.d0*ab(2) + at*( 6.d0*ab(3)                                    &
      &                 + at*(12.d0*ab(4) + at*(20.d0*ab(5))))
      afd  = ab(1) + at*(2.d0*ab(2) + at*(3.d0*ab(3)                         &
      &            + at*(4.d0*ab(4) + at*(5.d0*ab(5)))))
      af   = at*(ab(1) + at*(ab(2) + at*(ab(3) &
      &    + at*(ab(4) + at*ab(5)))))  + ab(6)

      addd = afdd*af/afd

      adt  = af/(afd-a05*addd)
!      adt  = af/afd

      at = at - adt
      if (abs(adt) < aerr) exit

    end do

    atmp = at

    aYi(:) = adns(:)/totaldens
    dtmp = atmp
    atmw(:) = aYi(:)*dmlr(:)      ! amlf(:)/atw
    druo    = dru*(sum(atmw(:)))
    itm = 1 + int(a05 + sign(a05,(atmp - dplt)))

    dprs   = totaldens*druo*atmp
    !totaldens = dprs*atmpr/druo
    !totaldens = sum(adns(:))

  deallocate(adns,afm,ab)

!-----------------------------------------------------------------------

  end subroutine pointimplicit

!=======================================================================
