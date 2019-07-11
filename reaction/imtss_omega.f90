subroutine imtss_omega(dtmp,dprs,aYi,delt,omega_ave)
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
  real(8), intent(inout) :: dtmp, dprs, aYi(9), omega_ave(9)
  real(8), intent(in) :: delt
  real(8), allocatable :: adns(:), afm(:), ab(:), dml(:), dmlr(:), dmlru(:),       &
       &                  dplh(:,:,:), atmw(:), adns0(:)

  real(8) :: at, afmt, aeng, afdd, afd, af, addd, adt, druo, phti, ahti
  real(8) :: atime, atmp, arho, aflag, art, atmpr, atmpl, aptr, adnsr
  real(8) :: ar, au, av, aw, awk1, awk2, awk3, totaldens
  real(8) :: ae01, ae02h, ae02l, ae03h, ae03l, ae04h, ae04l, ae05h, ae05l,    &
  &          ae06a, ae06b, ae07, ae08a, ae08b, ae09, ae10, ae11, ae12a, ae12b,&
  &          ae13, ae14, ae15, ae16, ae17, ae18, ae19, ae20, ae21, ae22, ae23,&
  &          ae24, ae25, akf01, akf02h, akf02l, akf03h, akf03l, akf04h,       &
  &          akf04l, akf05h, akf05l, akf06a, akf06b, akf07, akf08a, akf08b,   &
  &          akf09, akf10, akf11, akf12a, akf12b, akf13, akf14, akf15, akf16, &
  &          akf17, akf18, akf19, akf20, akf21, akf22, akf23, akf24, akf25,   &
  &          ah2, ao2, ah, ao, aoh, ah2o, aho2, ah2o2, an2, am04, am05, am09, &
  &          am21, am23, am24, am25, atrp, atrpl, atra, atrb,                 &
  &          atrf, akf02, akf03, akf04, akf05, akf06, akf08, akf12, as1, as2, &
  &          as3, as4, as5, as6, as7, as8, acon, ah1,  ah3, ah4, ah5, ah6,    &
  &          ah7, ah8, apt, ag, ahh01, ahh02, ahh03, ahh04, ahh05, ahh06,     &
  &          ahh07, ahh08, ahh09, ahh10, ahh11, ahh12, ahh13, ahh14, ahh15,   &
  &          ahh16, ahh17, ahh18, ahh19, ahh20, ahh21, ahh22, ahh23, ahh24,   &
  &          ahh25, akb01, akb02, akb03, akb04, akb05, akb06, akb07, akb08,   &
  &          akb09, akb10, akb11, akb12, akb13, akb14, akb15, akb16, akb17,   &
  &          akb18, akb19, akb20, akb21, akb22, akb23, akb24, akb25, akfc01,  &
  &          akfc02, akfc03, akfc04, akfc05, akfc06, akfc07, akfc08, akfc09,  &
  &          akfc10, akfc11, akfc12, akfc13, akfc14, akfc15, akfc16, akfc17,  &
  &          akfc18, akfc19, akfc20, akfc21, akfc22, akfc23, akfc24, akfc25,  &
  &          akbc01, akbc02, akbc03, akbc04, akbc05, akbc06, akbc07, akbc08,  &
  &          akbc09, akbc10, akbc11, akbc12, akbc13, akbc14, akbc15, akbc16,  &
  &          akbc17, akbc18, akbc19, akbc20, akbc21, akbc22, akbc23, akbc24,  &
  &          akbc25, ar01, ar02, ar03, ar04, ar05, ar06, ar07, ar08, ar09,    &
  &          ar10, ar11, ar12, ar13, ar14, ar15, ar16, ar17, ar18, ar19, ar20,&
  &          ar21, ar22, ar23, ar24, ar25, awh2, awo2, awh, awo, awoh, awh2o, &
  &          awho2, awh2o2, ath2, ato2, ath, ato, atoh, ath2o, atho2, ath2o2, &
  &          afg, adlt, adlt1, adlt2, asgm, atmx, aaa, abb, ajcb, aflgm1,     &
  &          aflgm2, aflgm3, aflgm4, aflgm5, aflgm6, aflgm7, aflgm8
  integer :: i, j, k, l, itm, m
! number of chemical species
  integer, parameter :: lsp = 9
! number of reaction constant
  integer, parameter :: lnr = 32
! constant parameter
  real(8), parameter :: arc = 1.9872d0, arp = 82.06d0
  real(8), parameter :: arcr = 1.d0 / arc
  real(8), parameter :: arur = 1.d0 / dru
! limiter
  real(8), parameter :: aemn = -130.d0, aemx = 130.d0
  real(8), parameter :: akmn = 1.d-35, akmx = 1.d+35
! MTS limiter
  real(8), parameter :: aber = 1.d-13, arer = 1.d-5
  real(8), parameter :: aeps = 1.d-35, atmn = 1.d-12
! reaction frozen temperature
  real(8), parameter :: a400 = 400.d0
  real(8), parameter :: aerr = 1.d-12
  integer, parameter :: itr = 10 ! number of iteration in Newton method
!
! constant parameter for Troe (for reduction)
  real(8), parameter :: atrc8 =  -0.4d0 - 0.67d0*log10(0.8d0)
  real(8), parameter :: atrn8 = 0.75d0 - 1.27d0*log10(0.8d0)
  real(8), parameter :: atrc7 = -0.4d0 - 0.67d0*log10(0.7d0)
  real(8), parameter :: atrn7 = 0.75d0 - 1.27d0*log10(0.7d0)

! reaction parameter for stanford model
  real(8), parameter :: acf(1:lnr) = (/ 1.04d+14, 5.59d+13, 3.70d+19,         &
  &                                     5.59d+13, 5.69d+18, 5.59d+13,         &
  &                                     2.65d+19, 8.59d+14, 9.55d+15,         &
  &                                     1.74d+12, 7.59d+13, 2.89d+13,         &
  &                                     1.30d+11, 4.20d+14, 6.06d+27,         &
  &                                     1.00d+26, 3.57d+4,  3.82d+12,         &
  &                                     8.79d+14, 2.17d+8,  7.08d+13,         &
  &                                     1.45d+12, 3.66d+6,  1.63d+13,         &
  &                                     1.21d+7,  1.02d+13, 8.43d+11,         &
  &                                     5.84d+18, 9.03d+14, 4.58d+19,         &
  &                                     6.16d+15, 4.71d+18 /)
  real(8), parameter :: ant(1:lnr) = (/ 0.0d0,  0.2d0, -1.0d0,  0.2d0,        &
  &                                    -1.1d0,  0.2d0, -1.3d0,  0.0d0,        &
  &                                     0.0d0,  0.0d0,  0.0d0,  0.0d0,        &
  &                                     0.0d0,  0.0d0,-3.31d0,-2.44d0,        &
  &                                     2.4d0,  0.0d0,  0.0d0, 1.52d0,        &
  &                                     0.0d0,  0.0d0,2.087d0,  0.0d0,        &
  &                                     2.0d0,  0.0d0,  0.0d0, -1.1d0,        &
  &                                     0.0d0, -1.4d0, -0.5d0, -1.0d0 /)
  real(8), parameter :: aea(1:lnr) = (/ 15286d0,  0.0d0,   0.0d0,  0.0d0,     &
  &                                       0.0d0,  0.0d0,   0.0d0,48560d0,     &
  &                                     42203d0,  318d0,  7269d0, -500d0,     &
  &                                     -1603d0,11980d0,120770d0,120160d0,    &
  &                                     -2111d0, 7948d0, 19170d0,  3457d0,    &
  &                                       300d0,  0.0d0, -1450d0, -445d0,     &
  &                                      5200d0, 3577d0, 3970d0,104380d0,     &
  &                                     96070d0,104380d0, 0.0d0,  0.0d0 /)

!///////////////////////////////////////////////////////////////////////
!//////////////////////  stanford reaction model   /////////////////////
!///////////////////////////////////////////////////////////////////////
!
!        1 :  h    + o2         = oh   + o
!        2 :  h    + o2   + h2o = ho2  + h2o        !pressure-dependent
!        3 :  h    + o2   + o2  = ho2  + o2         !pressure-dependent
!        4 :  h    + o2   + M   = ho2  + M          !pressure-dependent
!        5 :  h2o2        + M   = 2oh  + M          !pressure-dependent
!        6 :  oh   + h2o2       = h2o  + ho2        !non-Arrhenius
!        7 :  oh   + ho2        = h2o  + o2
!        8 :  ho2  + ho2        = h2o2 + o2         !non-Arrhenius
!        9 :  h2o         + M   = h    + oh  + M
!        10:  h2o         + h2o = oh   + h   + h2o
!        11:  oh   + oh         = h2o  + o
!        12:  o    + h2         = h    + oh         !non-Arrhenius
!        13:  h2   + oh         = h2o  + h
!        14:  h    + ho2        = oh   + oh
!        15:  h    + ho2        = h2o  + o
!        16:  h    + ho2        = h2   + o2
!        17:  o    + ho2        = oh   + o2
!        18:  h2o2 + h          = ho2  + h2
!        19:  h2o2 + h          = h2o  + oh
!        20:  h2o2 + o          = oh   + ho2
!        21:  h2          + M   = h    + h   + M
!        22:  h2          + h2  = h    + h   + h2
!        23:  h2          + M   = h    + h   + M    ! M = O2, N2
!        24:  o    + o    + M   = o2         + M
!        25:  o    + h    + M   = oh         + M
!
!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////

!  bdbg = 'on'

!----------------------------------------------------------------------
  if (bdbg(1:2) == 'on') write(6,*) 'in imtss '
!-----------------------------------------------------------------------

! allocate
  allocate(adns(lsp),afm(lsp),ab(6), dml(lsp), dmlr(lsp), dmlru(lsp),      &
  &        dplh(7,lsp,3), atmw(lsp),adns0(lsp))


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
  ah      = 0.d0
  atmw(:) = aYi(:)*dmlr(:)      ! amlf(:)/atw
  druo    = dru*(sum(atmw(:)))
  itm = 1 + int(a05 + sign(a05,(atmp - dplt)))

  !dprs   = totaldens*druo*atmp
  totaldens = dprs*atmpr/druo

  adns(:) = aYi(:)*totaldens
  adns0(:) = adns(:)
  arho    = sum(adns(1:lsp-1))
  adnsr   = 1.d0 / arho

  do i = 1, lsp
    call calc_phti(i,itm,dtmp,phti,dmlr,dplh)
    ahti = phti
    ah   = ah  + adns(i)*ahti/arho
  end do
  aeng = dprs - ah*arho

  aflag = a05 + sign(a05,(atmp - a400))

! flag for frozen species :: aflgmi (i species)
  aflgm1 = 1.d0 * aflag
  aflgm2 = 1.d0 * aflag
  aflgm3 = 1.d0 * aflag
  aflgm4 = 1.d0 * aflag
  aflgm5 = 1.d0 * aflag
  aflgm6 = 1.d0 * aflag
  aflgm7 = 1.d0 * aflag
  aflgm8 = 1.d0 * aflag

  do m = 1, 1000000    !MTS loop

    atmpr = 1.d0 / atmp

    art = arcr*atmpr
!      art  =  1.0d0/(arc*atmp)

!     [E/R]
!     activation energy

    ae01   = min(aemx,max(aemn,-aea( 1)*art))
    ae02h  = min(aemx,max(aemn,-aea( 2)*art))
    ae02l  = min(aemx,max(aemn,-aea( 3)*art))
    ae03h  = min(aemx,max(aemn,-aea( 4)*art))
    ae03l  = min(aemx,max(aemn,-aea( 5)*art))
    ae04h  = min(aemx,max(aemn,-aea( 6)*art))
    ae04l  = min(aemx,max(aemn,-aea( 7)*art))
    ae05h  = min(aemx,max(aemn,-aea( 8)*art))
    ae05l  = min(aemx,max(aemn,-aea( 9)*art))
    ae06a  = min(aemx,max(aemn,-aea(10)*art))
    ae06b  = min(aemx,max(aemn,-aea(11)*art))
    ae07   = min(aemx,max(aemn,-aea(12)*art))
    ae08a  = min(aemx,max(aemn,-aea(13)*art))
    ae08b  = min(aemx,max(aemn,-aea(14)*art))
    ae09   = min(aemx,max(aemn,-aea(15)*art))
    ae10   = min(aemx,max(aemn,-aea(16)*art))
    ae11   = min(aemx,max(aemn,-aea(17)*art))
    ae12a  = min(aemx,max(aemn,-aea(18)*art))
    ae12b  = min(aemx,max(aemn,-aea(19)*art))
    ae13   = min(aemx,max(aemn,-aea(20)*art))
    ae14   = min(aemx,max(aemn,-aea(21)*art))
    ae15   = min(aemx,max(aemn,-aea(22)*art))
    ae16   = min(aemx,max(aemn,-aea(23)*art))
    ae17   = min(aemx,max(aemn,-aea(24)*art))
    ae18   = min(aemx,max(aemn,-aea(25)*art))
    ae19   = min(aemx,max(aemn,-aea(26)*art))
    ae20   = min(aemx,max(aemn,-aea(27)*art))
    ae21   = min(aemx,max(aemn,-aea(28)*art))
    ae22   = min(aemx,max(aemn,-aea(29)*art))
    ae23   = min(aemx,max(aemn,-aea(30)*art))
    ae24   = min(aemx,max(aemn,-aea(31)*art))
    ae25   = min(aemx,max(aemn,-aea(32)*art))


!     [kf]
!     forward reaction rate constant
!     Arrhenius-s form : kf=AT**n exp(-Ea/RT)

    akf01  = acf( 1)*exp(ae01)                ! n = 0.0d0
    akf02h = acf( 2)*atmp**ant( 2)*exp(ae02h) ! n = 0.2d0
    akf02l = acf( 3)*atmpr*exp(ae02l)         ! n = -1.d0
    akf03h = acf( 4)*atmp**ant( 4)*exp(ae03h) ! n = 0.2d0
    akf03l = acf( 5)*atmp**ant( 5)*exp(ae03l) ! n = -1.1d0
    akf04h = acf( 6)*atmp**ant( 6)*exp(ae04h) ! n = 0.2d0
    akf04l = acf( 7)*atmp**ant( 7)*exp(ae04l) ! n = -1.3d0
    akf05h = acf( 8)*exp(ae05h)               ! n = 0.d0
    akf05l = acf( 9)*exp(ae05l)               ! n = 0.d0
    akf06a = acf(10)*exp(ae06a)               ! n = 0.d0
    akf06b = acf(11)*exp(ae06b)               ! n = 0.d0
    akf07  = acf(12)*exp(ae07)                ! n = 0.d0
    akf08a = acf(13)*exp(ae08a)               ! n = 0.d0
    akf08b = acf(14)*exp(ae08b)               ! n = 0.d0
    akf09  = acf(15)*atmp**ant(15)*exp(ae09)  ! n = -3.31d0
    akf10  = acf(16)*atmp**ant(16)*exp(ae10)  ! n = -2.44d0
    akf11  = acf(17)*atmp**ant(17)*exp(ae11)  ! n = 2.4d0
    akf12a = acf(18)*exp(ae12a)               ! n = 0.d0
    akf12b = acf(19)*exp(ae12b)               ! n = 0.d0
    akf13  = acf(20)*atmp**ant(20)*exp(ae13)  ! n = 1.52d0
    akf14  = acf(21)*exp(ae14)                ! n = 0.d0
    akf15  = acf(22)*exp(ae15)                ! n = 0.d0
    akf16  = acf(23)*atmp**ant(23)*exp(ae16)  ! n = 2.087d0
    akf17  = acf(24)*exp(ae17)                ! n = 0.d0
    akf18  = acf(25)*atmp*atmp*exp(ae18)      ! n = 2.d0
    akf19  = acf(26)*exp(ae19)                ! n = 0.d0
    akf20  = acf(27)*exp(ae20)                ! n = 0.d0
    akf21  = acf(28)*atmp**ant(28)*exp(ae21)  ! n = -1.1d0
    akf22  = acf(29)*exp(ae22)                ! n = 0.d0
    akf23  = acf(30)*atmp**ant(30)*exp(ae23)  ! n = -1.4d0
    akf24  = acf(31)*sqrt(atmpr)*exp(ae24)    ! n = -0.5d0
    akf25  = acf(32)*atmpr*exp(ae25)          ! n = -1.d0


!      akf01  = acf( 1)*atmp**ant( 1)*exp(ae01)
!      akf02h = acf( 2)*atmp**ant( 2)*exp(ae02h)
!      akf02l = acf( 3)*atmp**ant( 3)*exp(ae02l)
!      akf03h = acf( 4)*atmp**ant( 4)*exp(ae03h)
!      akf03l = acf( 5)*atmp**ant( 5)*exp(ae03l)
!      akf04h = acf( 6)*atmp**ant( 6)*exp(ae04h)
!      akf04l = acf( 7)*atmp**ant( 7)*exp(ae04l)
!      akf05h = acf( 8)*atmp**ant( 8)*exp(ae05h)
!      akf05l = acf( 9)*atmp**ant( 9)*exp(ae05l)
!      akf06a = acf(10)*atmp**ant(10)*exp(ae06a)
!      akf06b = acf(11)*atmp**ant(11)*exp(ae06b)
!      akf07  = acf(12)*atmp**ant(12)*exp(ae07)
!      akf08a = acf(13)*atmp**ant(13)*exp(ae08a)
!      akf08b = acf(14)*atmp**ant(14)*exp(ae08b)
!      akf09  = acf(15)*atmp**ant(15)*exp(ae09)
!      akf10  = acf(16)*atmp**ant(16)*exp(ae10)
!      akf11  = acf(17)*atmp**ant(17)*exp(ae11)
!      akf12a = acf(18)*atmp**ant(18)*exp(ae12a)
!      akf12b = acf(19)*atmp**ant(19)*exp(ae12b)
!      akf13  = acf(20)*atmp**ant(20)*exp(ae13)
!      akf14  = acf(21)*atmp**ant(21)*exp(ae14)
!      akf15  = acf(22)*atmp**ant(22)*exp(ae15)
!      akf16  = acf(23)*atmp**ant(23)*exp(ae16)
!      akf17  = acf(24)*atmp**ant(24)*exp(ae17)
!      akf18  = acf(25)*atmp**ant(25)*exp(ae18)
!      akf19  = acf(26)*atmp**ant(26)*exp(ae19)
!      akf20  = acf(27)*atmp**ant(27)*exp(ae20)
!      akf21  = acf(28)*atmp**ant(28)*exp(ae21)
!      akf22  = acf(29)*atmp**ant(29)*exp(ae22)
!      akf23  = acf(30)*atmp**ant(30)*exp(ae23)
!      akf24  = acf(31)*atmp**ant(31)*exp(ae24)
!      akf25  = acf(32)*atmp**ant(32)*exp(ae25)

!     [Xi]
!     concentration of species Xi
!     mole density

      ah2   = adns(1) * dmlr(1)
      ao2   = adns(2) * dmlr(2)
      ah    = adns(3) * dmlr(3)
      ao    = adns(4) * dmlr(4)
      aoh   = adns(5) * dmlr(5)
      ah2o  = adns(6) * dmlr(6)
      aho2  = adns(7) * dmlr(7)
      ah2o2 = adns(8) * dmlr(8)
      an2   = adns(9) * dmlr(9)

!    [M]
!    Third Body

      am04  = ah2*1.5d0  +ao2*0.d0  +ah  +ao  +aoh  +ah2o*0.d0  +aho2         &
  &         + ah2o2      +an2
      am05  = ah2        +ao2      +ah  +ao  +aoh  +ah2o*9.0d0  +aho2         &
  &         + ah2o2      +an2*1.5d0
      am09  = ah2*3.0d0  +ao2*1.5d0  +ah  +ao  +aoh  +ah2o*0.0d0              &
  &         + aho2  + ah2o2    +an2*2.0d0
      am21  = ah2*0.0d0  +ao2*0.0d0  +ah  +ao  +aoh  +ah2o*14.4d0             &
  &         + aho2 + ah2o2    +an2*0.0d0
      am23  = ao2        +an2
      am24  = ah2*2.5d0  +ao2      +ah  +ao  +aoh  +ah2o*12d0   +aho2         &
  &         + ah2o2      +an2
      am25  = ah2*2.5d0  +ao2      +ah  +ao  +aoh  +ah2o*12d0   +aho2         &
  &         + ah2o2      +an2

!-----------------------------------------------------------
!  the reaction rate constant of #02

      atrp = akf02l*ah2o*1.0d-3/akf02h

!      atrc = -0.4d0 - 0.67d0*log10(0.8d0)
!      atrn = 0.75d0 - 1.27d0*log10(0.8d0)

      atrpl= log10(atrp+1.0d-30)
      atra = (atrpl+atrc8)/(atrn8-0.14d0*(atrpl+atrc8))
      atrb = 1.d0/(1.d0 + atra*atra)
      atrf = 0.8d0**atrb

      akf02 = akf02h*atrp/(1.d0+atrp)*atrf
!----------------------------------------------------------
!  the reaction rate constant of #03

      atrp = akf03l*ao2*1.0d-3/akf03h

!      atrc = -0.4d0 - 0.67d0*log10(0.7d0)
!      atrn = 0.75d0 - 1.27d0*log10(0.7d0)

      atrpl= log10(atrp+1.0d-30)
      atra = (atrpl+atrc7)/(atrn7-0.14d0*(atrpl+atrc7))
      atrb = 1.d0/(1.d0 + atra*atra)
      atrf = 0.7d0**atrb

      akf03 = akf03h*atrp/(1.d0+atrp)*atrf

!----------------------------------------------------------
!  the reaction rate constant of #04

      atrp = akf04l*am04*1.0d-3/akf04h

!      atrc = -0.4d0 - 0.67d0*log10(0.7d0)
!      atrn = 0.75d0 - 1.27d0*log10(0.7d0)

      atrpl= log10(atrp+1.0d-30)
      atra = (atrpl+atrc7)/(atrn7-0.14d0*(atrpl+atrc7))
      atrb = 1.d0/(1.d0 + atra*atra)
      atrf = 0.7d0**atrb

      akf04 = akf04h*atrp/(1.d0+atrp)*atrf

!----------------------------------------------------------
!  the reaction rate constant of #05

      atrp = akf05l*am05*1.0d-3/akf05h

!         atrc = -0.4 - 0.67*log10(0.8)
!         atrn = 0.75 - 1.27*log10(0.8)

!         atrpl= log10(atrp)
!         atra = (atrpl+atrc)/(atrn-0.14*(atrpl+atrc))
!         atrb = 1.0/(1.0 + atra*atra)
!         atrf = 0.7**atrb

!         akf05 = akf05h*atrp/(1.0+atrp)*atrf
      akf05 = akf05h*atrp/(1.d0+atrp)
!-----------------------------------------------------------
!  the reaction rate constant of #06 & #08 & #12

      akf06 = akf06a + akf06b
      akf08 = akf08a + akf08b
      akf12 = akf12a + akf12b

!------------------------------------------------------------

!    [kb]
!    backward reaction rate constant

! From the analysis in Intel Adviser,
! using user function in the loop takes longer time to calculate.
! In order to improve calculation performance,
! I refrain from using user function.
!
! Demerit is that the code is not beautiful without using user function

      itm = 1 + int(a05 + sign(a05,(atmp - dplt)))

      atmpl = log(atmp)

!..01
      as1 =     dplh(1,ls1,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls1,itm)     + atmp*(dplh(3,ls1,itm)*a12             &
  &       + atmp*(dplh(4,ls1,itm)*a13 + atmp*(dplh(5,ls1,itm)*a14))))         &
  &       + dplh(7,ls1,itm)
!..02
      as2 =     dplh(1,ls2,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls2,itm)     + atmp*(dplh(3,ls2,itm)*a12             &
  &       + atmp*(dplh(4,ls2,itm)*a13 + atmp*(dplh(5,ls2,itm)*a14))))         &
  &       + dplh(7,ls2,itm)
!..03
      as3 =     dplh(1,ls3,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls3,itm)     + atmp*(dplh(3,ls3,itm)*a12             &
  &       + atmp*(dplh(4,ls3,itm)*a13 + atmp*(dplh(5,ls3,itm)*a14))))         &
  &       + dplh(7,ls3,itm)
!..04
      as4 =     dplh(1,ls4,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls4,itm)     + atmp*(dplh(3,ls4,itm)*a12             &
  &       + atmp*(dplh(4,ls4,itm)*a13 + atmp*(dplh(5,ls4,itm)*a14))))         &
  &       + dplh(7,ls4,itm)
!..05
      as5 =     dplh(1,ls5,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls5,itm)     + atmp*(dplh(3,ls5,itm)*a12             &
  &       + atmp*(dplh(4,ls5,itm)*a13 + atmp*(dplh(5,ls5,itm)*a14))))         &
  &       + dplh(7,ls5,itm)
!..06
      as6 =     dplh(1,ls6,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls6,itm)     + atmp*(dplh(3,ls6,itm)*a12             &
  &       + atmp*(dplh(4,ls6,itm)*a13 + atmp*(dplh(5,ls6,itm)*a14))))         &
  &       + dplh(7,ls6,itm)
!..07
      as7 =     dplh(1,ls7,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls7,itm)     + atmp*(dplh(3,ls7,itm)*a12             &
  &       + atmp*(dplh(4,ls7,itm)*a13 + atmp*(dplh(5,ls7,itm)*a14))))         &
  &       + dplh(7,ls7,itm)
!..08
      as8 =     dplh(1,ls8,itm)*atmpl                                         &
  &       + atmp*(dplh(2,ls8,itm)     + atmp*(dplh(3,ls8,itm)*a12             &
  &       + atmp*(dplh(4,ls8,itm)*a13 + atmp*(dplh(5,ls8,itm)*a14))))         &
  &       + dplh(7,ls8,itm)

      acon = atmpr
!     acon = arur*atmpr
!     acon = 1.0d0/(dru*atmp)

!..01
      ah1 = (atmp*(dplh(1,ls1,itm)                                            &
  &         + atmp*(dplh(2,ls1,itm)*a12                                       &
  &         + atmp*(dplh(3,ls1,itm)*a13                                       &
  &         + atmp*(dplh(4,ls1,itm)*a14                                       &
  &         + atmp*(dplh(5,ls1,itm)*a15)))))                                  &
  &         +     dplh(6,ls1,itm))                                            &
  &         * acon
!..02
      ah2 = (atmp*(dplh(1,ls2,itm)                                            &
  &         + atmp*(dplh(2,ls2,itm)*a12                                       &
  &         + atmp*(dplh(3,ls2,itm)*a13                                       &
  &         + atmp*(dplh(4,ls2,itm)*a14                                       &
  &         + atmp*(dplh(5,ls2,itm)*a15)))))                                  &
  &         +     dplh(6,ls2,itm))                                            &
  &         * acon
!..03
      ah3 = (atmp*(dplh(1,ls3,itm)                                            &
  &         + atmp*(dplh(2,ls3,itm)*a12                                       &
  &         + atmp*(dplh(3,ls3,itm)*a13                                       &
  &         + atmp*(dplh(4,ls3,itm)*a14                                       &
  &         + atmp*(dplh(5,ls3,itm)*a15)))))                                  &
  &         +     dplh(6,ls3,itm))                                            &
  &         * acon
!..04
      ah4 = (atmp*(dplh(1,ls4,itm)                                            &
  &         + atmp*(dplh(2,ls4,itm)*a12                                       &
  &         + atmp*(dplh(3,ls4,itm)*a13                                       &
  &         + atmp*(dplh(4,ls4,itm)*a14                                       &
  &         + atmp*(dplh(5,ls4,itm)*a15)))))                                  &
  &         +     dplh(6,ls4,itm))                                            &
  &         * acon
!..05
      ah5 = (atmp*(dplh(1,ls5,itm)                                            &
  &         + atmp*(dplh(2,ls5,itm)*a12                                       &
  &         + atmp*(dplh(3,ls5,itm)*a13                                       &
  &         + atmp*(dplh(4,ls5,itm)*a14                                       &
  &         + atmp*(dplh(5,ls5,itm)*a15)))))                                  &
  &         +     dplh(6,ls5,itm))                                            &
  &         * acon
!..06
      ah6 = (atmp*(dplh(1,ls6,itm)                                            &
  &         + atmp*(dplh(2,ls6,itm)*a12                                       &
  &         + atmp*(dplh(3,ls6,itm)*a13                                       &
  &         + atmp*(dplh(4,ls6,itm)*a14                                       &
  &         + atmp*(dplh(5,ls6,itm)*a15)))))                                  &
  &         +     dplh(6,ls6,itm))                                            &
  &         * acon
!..07
      ah7 = (atmp*(dplh(1,ls7,itm)                                            &
  &         + atmp*(dplh(2,ls7,itm)*a12                                       &
  &         + atmp*(dplh(3,ls7,itm)*a13                                       &
  &         + atmp*(dplh(4,ls7,itm)*a14                                       &
  &         + atmp*(dplh(5,ls7,itm)*a15)))))                                  &
  &         +     dplh(6,ls7,itm))                                            &
  &         * acon
!..08
      ah8 = (atmp*(dplh(1,ls8,itm)                                            &
  &         + atmp*(dplh(2,ls8,itm)*a12                                       &
  &         + atmp*(dplh(3,ls8,itm)*a13                                       &
  &         + atmp*(dplh(4,ls8,itm)*a14                                       &
  &         + atmp*(dplh(5,ls8,itm)*a15)))))                                  &
  &         +     dplh(6,ls8,itm))                                            &
  &         * acon

      apt = arp*atmp
      aptr = 1.d0 / apt


!..01
      ahh01 = ah5+ah4-ah3-ah2
      ag    = (as5+as4-as3-as2)-ahh01
      akb01 = akf01/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..02
      ahh02 = ah7    -ah3-ah2
      ag    = (as7    -as3-as2)-ahh02
      akb02 = akf02/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * apt))
!..03
      ahh03 = ah7    -ah3-ah2
      ag    = (as7    -as3-as2)-ahh03
      akb03 = akf03/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * apt))
!..04
      ahh04 = ah7    -ah3-ah2
      ag    = (as7   -as3-as2)-ahh04
      akb04 = akf04/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * apt))
!..05
      ahh05 = ah5+ah5    -ah8
      ag    = (as5+as5    -as8)-ahh05
      akb05 = akf05/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * aptr))
!  &         / apt))
!..06
      ahh06 = ah6+ah7-ah5-ah8
      ag    = (as6+as7-as5-as8)-ahh06
      akb06 = akf06/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..07
      ahh07 = ah6+ah2-ah5-ah7
      ag    = (as6+as2-as5-as7)-ahh07
      akb07 = akf07/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..08
      ahh08 = ah8+ah2-ah7-ah7
      ag    = (as8+as2-as7-as7)-ahh08
      akb08 = akf08/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..09
      ahh09 = ah3+ah5    -ah6
      ag    = (as3+as5    -as6)-ahh09
      akb09 = akf09/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * aptr))
!  &         / apt))
!..10
      ahh10 = ah5+ah3    -ah6
      ag    = (as5+as3    -as6)-ahh10
      akb10 = akf10/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * aptr))
!  &         / apt))
!..11
      ahh11 = ah6+ah4-ah5-ah5
      ag    = (as6+as4-as5-as5)-ahh11
      akb11 = akf11/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..12
      ahh12 = ah3+ah5-ah4-ah1
      ag    = (as3+as5-as4-as1)-ahh12
      akb12 = akf12/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..13
      ahh13 = ah6+ah3-ah1-ah5
      ag    = (as6+as3-as1-as5)-ahh13
      akb13 = akf13/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..14
      ahh14 = ah5+ah5-ah3-ah7
      ag    = (as5+as5-as3-as7)-ahh14
      akb14 = akf14/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..15
      ahh15 = ah6+ah4-ah3-ah7
      ag    = (as6+as4-as3-as7)-ahh15
      akb15 = akf15/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..16
      ahh16 = ah1+ah2-ah3-ah7
      ag    = (as1+as2-as3-as7)-ahh16
      akb16 = akf16/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..17
      ahh17 = ah5+ah2-ah4-ah7
      ag    = (as5+as2-as4-as7)-ahh17
      akb17 = akf17/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..18
      ahh18 = ah7+ah1-ah8-ah3
      ag    = (as7+as1-as8-as3)-ahh18
      akb18 = akf18/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..19
      ahh19 = ah6+ah5-ah8-ah3
      ag    = (as6+as5-as8-as3)-ahh19
      akb19 = akf19/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..20
      ahh20 = ah5+ah7-ah8-ah4
      ag    = (as5+as7-as8-as4)-ahh20
      akb20 = akf20/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..21
      ahh21 = ah3+ah3    -ah1
      ag    = (as3+as3    -as1)-ahh21
      akb21 = akf21/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * aptr))
!  &         / apt))
!..22
      ahh22 = ah3+ah3    -ah1
      ag    = (as3+as3    -as1)-ahh22
      akb22 = akf22/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * aptr))
!  &         / apt))
!..23
      ahh23 = ah3+ah3    -ah1
      ag    = (as3+as3    -as1)-ahh23
      akb23 = akf23/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * aptr))
!  &         / apt))
!..24
      ahh24 = ah2    -ah4-ah4
      ag    = (as2    -as4-as4)-ahh24
      akb24 = akf24/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * apt))
!..25
      ahh25 = ah5    -ah4-ah3
      ag    = (as5    -as4-as3)-ahh25
      akb25 = akf25/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
  &         * apt))

!-----------------------------------------------------------------------
! adjust unit

      akf01 = akf01 * 1.d-3
      akf02 = akf02 * 1.d-3
      akf03 = akf03 * 1.d-3
      akf04 = akf04 * 1.d-3
!      akf05 = akf05 * 1.d-3
      akf06 = akf06 * 1.d-3
      akf07 = akf07 * 1.d-3
      akf08 = akf08 * 1.d-3
      akf09 = akf09 * 1.d-3
      akf10 = akf10 * 1.d-3
      akf11 = akf11 * 1.d-3
      akf12 = akf12 * 1.d-3
      akf13 = akf13 * 1.d-3
      akf14 = akf14 * 1.d-3
      akf15 = akf15 * 1.d-3
      akf16 = akf16 * 1.d-3
      akf17 = akf17 * 1.d-3
      akf18 = akf18 * 1.d-3
      akf19 = akf19 * 1.d-3
      akf20 = akf20 * 1.d-3
      akf21 = akf21 * 1.d-3
      akf22 = akf22 * 1.d-3
      akf23 = akf23 * 1.d-3
      akf24 = akf24 * 1.d-6
      akf25 = akf25 * 1.d-6

      akb01 = akb01 * 1.d-3
!      akb02 = akb02 * 1.d-3
!      akb03 = akb03 * 1.d-3
!      akb04 = akb04 * 1.d-3
      akb05 = akb05 * 1.d-3
      akb06 = akb06 * 1.d-3
      akb07 = akb07 * 1.d-3
      akb08 = akb08 * 1.d-3
      akb09 = akb09 * 1.d-6
      akb10 = akb10 * 1.d-6
      akb11 = akb11 * 1.d-3
      akb12 = akb12 * 1.d-3
      akb13 = akb13 * 1.d-3
      akb14 = akb14 * 1.d-3
      akb15 = akb15 * 1.d-3
      akb16 = akb16 * 1.d-3
      akb17 = akb17 * 1.d-3
      akb18 = akb18 * 1.d-3
      akb19 = akb19 * 1.d-3
      akb20 = akb20 * 1.d-3
      akb21 = akb21 * 1.d-6
      akb22 = akb22 * 1.d-6
      akb23 = akb23 * 1.d-6
      akb24 = akb24 * 1.d-3
      akb25 = akb25 * 1.d-3

!-----------------------------------------------------------------------

      ah2   = adns(1)*dmlr(1)
!      ao2   = dflw(j,k,l2)/dml(2)
!      ah    = dflw(j,k,l3)/dml(3)
!      ao    = dflw(j,k,l4)/dml(4)
!      aoh   = dflw(j,k,l5)/dml(5)
!      ah2o  = dflw(j,k,l6)/dml(6)
!      aho2  = dflw(j,k,l7)/dml(7)
!      ah2o2 = dflw(j,k,l8)/dml(8)
!      an2   = dflw(j,k,l9)/dml(9)

!    make  kf[ ][ ] : rection rate constant * mole density

      akfc01 = akf01 * ah    * ao2
      akfc02 = akf02 * ah    * ao2
      akfc03 = akf03 * ah    * ao2
      akfc04 = akf04 * ah    * ao2
      akfc05 = akf05 * ah2o2
      akfc06 = akf06 * aoh   * ah2o2
      akfc07 = akf07 * aoh   * aho2
      akfc08 = akf08 * aho2  * aho2
      akfc09 = akf09 * ah2o          * am09
      akfc10 = akf10 * ah2o          * ah2o
      akfc11 = akf11 * aoh   * aoh
      akfc12 = akf12 * ao    * ah2
      akfc13 = akf13 * ah2   * aoh
      akfc14 = akf14 * ah    * aho2
      akfc15 = akf15 * ah    * aho2
      akfc16 = akf16 * ah    * aho2
      akfc17 = akf17 * ao    * aho2
      akfc18 = akf18 * ah2o2 * ah
      akfc19 = akf19 * ah2o2 * ah
      akfc20 = akf20 * ah2o2 * ao
      akfc21 = akf21 * ah2           * am21
      akfc22 = akf22 * ah2           * ah2
      akfc23 = akf23 * ah2           * am23
      akfc24 = akf24 * ao    * ao    * am24
      akfc25 = akf25 * ao    * ah    * am25

!    make kb[ ][ ] : rection rate constant * mole density

      akbc01 = akb01 * aoh   * ao
      akbc02 = akb02 * aho2
      akbc03 = akb03 * aho2
      akbc04 = akb04 * aho2
      akbc05 = akb05 * aoh   * aoh
      akbc06 = akb06 * ah2o  * aho2
      akbc07 = akb07 * ah2o  * ao2
      akbc08 = akb08 * ah2o2 * ao2
      akbc09 = akb09 * ah    * aoh   * am09
      akbc10 = akb10 * aoh   * ah    * ah2o
      akbc11 = akb11 * ah2o  * ao
      akbc12 = akb12 * ah    * aoh
      akbc13 = akb13 * ah2o  * ah
      akbc14 = akb14 * aoh   * aoh
      akbc15 = akb15 * ah2o  * ao
      akbc16 = akb16 * ah2   * ao2
      akbc17 = akb17 * aoh   * ao2
      akbc18 = akb18 * aho2  * ah2
      akbc19 = akb19 * ah2o  * aoh
      akbc20 = akb20 * aoh   * aho2
      akbc21 = akb21 * ah    * ah    * am21
      akbc22 = akb22 * ah    * ah    * ah2
      akbc23 = akb23 * ah    * ah    * am23
      akbc24 = akb24 * ao2           * am24
      akbc25 = akb25 * aoh           * am25

!    reaction rate of each elementary reaction
!    kf[ ][ ] - kb[ ][ ]

      ar01 = akfc01 - akbc01
      ar02 = akfc02 - akbc02
      ar03 = akfc03 - akbc03
      ar04 = akfc04 - akbc04
      ar05 = akfc05 - akbc05
      ar06 = akfc06 - akbc06
      ar07 = akfc07 - akbc07
      ar08 = akfc08 - akbc08
      ar09 = akfc09 - akbc09
      ar10 = akfc10 - akbc10
      ar11 = akfc11 - akbc11
      ar12 = akfc12 - akbc12
      ar13 = akfc13 - akbc13
      ar14 = akfc14 - akbc14
      ar15 = akfc15 - akbc15
      ar16 = akfc16 - akbc16
      ar17 = akfc17 - akbc17
      ar18 = akfc18 - akbc18
      ar19 = akfc19 - akbc19
      ar20 = akfc20 - akbc20
      ar21 = akfc21 - akbc21
      ar22 = akfc22 - akbc22
      ar23 = akfc23 - akbc23
      ar24 = akfc24 - akbc24
      ar25 = akfc25 - akbc25

!    reaction rate of each species

      awh2   = dml(1) * (-ar12 -ar13 +ar16 +ar18 -ar21 -ar22 -ar23)
      awo2   = dml(2) * (-ar01 -ar02 -ar03 -ar04 +ar07 +ar08 +ar16            &
  &                      +ar17 +ar24)
      awh    = dml(3) * (-ar01 -ar02 -ar03 -ar04 +ar09 +ar10 +ar12            &
  &                      +ar13 -ar14 -ar15 -ar16 -ar18 -ar19 +ar21            &
  &                      +ar21 +ar22 +ar22 +ar23 +ar23 -ar25)
      awo    = dml(4) * (+ar01 +ar11 -ar12 +ar15 -ar17 -ar20 -ar24            &
  &                      -ar24 -ar25)
      awoh   = dml(5) * (+ar01 +ar05 +ar05 -ar06 -ar07 +ar09 +ar10            &
  &                      -ar11 -ar11 +ar12 -ar13 +ar14 +ar14 +ar17            &
  &                      +ar19 +ar20 +ar25)
      awh2o  = dml(6) * (+ar06 +ar07 -ar09 -ar10 +ar11 +ar13 +ar15            &
  &                      +ar19)
      awho2  = dml(7) * (+ar02 +ar03 +ar04 +ar06 -ar07 -ar08 -ar08            &
  &                      -ar14 -ar15 -ar16 -ar17 +ar18 +ar20)
      awh2o2 = dml(8) * (-ar05 -ar06 +ar08 -ar18 -ar19 -ar20)

!----------------------------------------------------------------------
!    MTS method

!    make reaction time scale
      ath2  = adns(1)*1.d-3*aflgm1/abs(awh2+aeps)                             &
  &         + 1.d0-aflgm1
      ato2  = adns(2)*1.d-3*aflgm2/abs(awo2+aeps)                             &
  &         + 1.d0-aflgm2
      ath   = adns(3)*1.d-3*aflgm3/abs(awh+aeps)                              &
  &         + 1.d0-aflgm3
      ato   = adns(4)*1.d-3*aflgm4/abs(awo+aeps)                              &
  &         + 1.d0-aflgm4
      atoh  = adns(5)*1.d-3*aflgm5/abs(awoh+aeps)                             &
  &         + 1.d0-aflgm5
      ath2o = adns(6)*1.d-3*aflgm6/abs(awh2o+aeps)                            &
  &         + 1.d0-aflgm6
      atho2 = adns(7)*1.d-3*aflgm7/abs(awho2+aeps)                            &
  &         + 1.d0-aflgm7
      ath2o2= adns(8)*1.d-3*aflgm8/abs(awh2o2+aeps)                           &
  &         + 1.d0-aflgm8

!    ignore uninvolved chemical species
      afg   = a05 + sign(a05, ath2-atmn)
      ath2  = ath2  *afg + 1.d0-afg
      afg   = a05 + sign(a05, ato2-atmn)
      ato2  = ato2  *afg + 1.d0-afg
      afg   = a05 + sign(a05, ath -atmn)
      ath   = ath   *afg + 1.d0-afg
      afg   = a05 + sign(a05, ato -atmn)
      ato   = ato   *afg + 1.d0-afg
      afg   = a05 + sign(a05, atoh-atmn)
      atoh  = atoh  *afg + 1.d0-afg
      afg   = a05 + sign(a05, ath2o-atmn)
      ath2o = ath2o *afg + 1.d0-afg
      afg   = a05 + sign(a05, atho2-atmn)
      atho2 = atho2 *afg + 1.d0-afg
      afg   = a05 + sign(a05, ath2o2-atmn)
      ath2o2= ath2o2*afg + 1.d0-afg


!    make minimum reaction time scale
      adlt1 = 0.8d0*min(ath2,ato2,ath,ato,atoh,ath2o,atho2,ath2o2)
      adlt2 = max(adlt1,atmn)
      atmx = delt - atime
      adlt = min(adlt2,atmx)

      awh2  = adlt*awh2
      awo2  = adlt*awo2
      awh   = adlt*awh
      awo   = adlt*awo
      awoh  = adlt*awoh
      awh2o = adlt*awh2o
      awho2 = adlt*awho2
      awh2o2= adlt*awh2o2

      adns(1) = adns(1) + aflgm1*awh2
      adns(2) = adns(2) + aflgm2*awo2
      adns(3) = adns(3) + aflgm3*awh
      adns(4) = adns(4) + aflgm4*awo
      adns(5) = adns(5) + aflgm5*awoh
      adns(6) = adns(6) + aflgm6*awh2o
      adns(7) = adns(7) + aflgm7*awho2
      adns(8) = adns(8) + aflgm8*awh2o2

!---------------------------------------------------------------------
!    revision

      asgm = 0.0d0

      adns(1:lsp-1) = max(0.d0,adns(1:lsp-1)) ! mass correction

      asgm = sum(adns(1:lsp-1))

      asgm = arho / asgm

      adns(1:lsp-1) = adns(1:lsp-1) * asgm


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

        afdd = 2.d0*ab(2) + at*( 6.d0*ab(3)                                  &
        &                 + at*(12.d0*ab(4) + at*(20.d0*ab(5))))
        afd  = ab(1) + at*(2.d0*ab(2) + at*(3.d0*ab(3)                       &
        &            + at*(4.d0*ab(4) + at*(5.d0*ab(5)))))
        af   = at*(ab(1) + at*(ab(2) + at*(ab(3) &
        &    + at*(ab(4) + at*ab(5)))))  + ab(6)

        addd = afdd*af/afd

        adt  = af/(afd-a05*addd)
!        adt  = af/afd

        at = at - adt
        if (abs(adt) < aerr) exit

      end do

      atmp = at

!-----------------------------------------
! update aflg(i) (, which represents frozen species)
! -----------------------------------------

    aaa = abs(awh2)-aber
    abb = abs(awo2)/(adns(2)-awo2+aeps) - arer
    aflgm2 = a05+ sign(a05,aaa*abb*min(aaa,abb))

    aaa = abs(awh)-aber
    abb = abs(awh)/(adns(3)-awh+aeps) - arer
    aflgm3 = a05+ sign(a05,aaa*abb*min(aaa,abb))

    aaa = abs(awo)-aber
    abb = abs(awo)/(adns(4)-awo+aeps) - arer
    aflgm4 = a05+ sign(a05,aaa*abb*min(aaa,abb))

    aaa = abs(awoh)-aber
    abb = abs(awoh)/(adns(5)-awoh+aeps) - arer
    aflgm5 = a05+ sign(a05,aaa*abb*min(aaa,abb))

    aaa = abs(awh2o)-aber
    abb = abs(awh2o)/(adns(6)-awh2o+aeps) - arer
    aflgm6 = a05+ sign(a05,aaa*abb*min(aaa,abb))

    aaa = abs(awho2)-aber
    abb = abs(awho2)/(adns(7)-awho2+aeps) - arer
    aflgm7 = a05+ sign(a05,aaa*abb*min(aaa,abb))

    aaa = abs(awh2o2)-aber
    abb = abs(awh2o2)/(adns(8)-awh2o2+aeps) - arer
    aflgm8 = a05+ sign(a05,aaa*abb*min(aaa,abb))

!-----------------------------------------------------

    atime = atime + adlt

    if (atime>=delt)exit

    end do

    omega_ave(:) = (adns(:)-adns0(:))/delt
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

  end subroutine imtss_omega

!=======================================================================
