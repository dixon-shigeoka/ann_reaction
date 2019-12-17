!---------------------------------------------------------------------------
  subroutine main
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     this program is created for estimating detonation profile
!     this code was written by hiroaki watanabe, 2014/11/11
!     converted to f90 and added wk model by shunsuke shigeoka 2018/5/
!
!     Target            : Hydrogen detonaion
!
!     Stanford model
!     Reference: Hong, Z., Davidson, D. F., Hanson, R. K.,
!                "An improved H2/O2 mechanism based on recent shocktube/
!                 laser absorption measurements,"
!                Combustion and Flame, 158 (2011) 633-644
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  real(8), parameter :: dplt = 1000.d0
  real(8), parameter :: dru  = 8314.4d0
  real(8), parameter :: a05  = 0.5d0
  real(8), parameter :: epsl = 1.d-3
  integer, parameter :: jmax = 1000000001
  integer, parameter :: l1=1, l2=2, l3=3, l4=4, l5=5, l6=6, l7=7, l8=8, l9=9, &
       &                     lfl=3, lsp=9, lh=10, la=11, lr=lsp+1, lu=lsp+2,       &
       &                     le=lsp+3, lt=lsp+4, lp=lsp+5, lg=lsp+6, imax=lg,      &
       &                     lmss=1, lmmt=2, leng=3, iplus=3, iimax=imax+iplus
  real(8), allocatable, save :: dflw(:), dcnv(:), dmh2(:), dtime(:), dml(:),  &
       &                             dplh(:,:,:)
  real(8), save :: dcvt, dcuv, ddltx
  real(8) :: ddgd, au0, dmc
  integer :: jend, lstp, lint
  character :: bdbg*4, brct*6, btyp*3
  !icond
  character :: crct*4
  real(8) :: axx
  integer :: ijj
  !ioutp
  real(8), intent(in) :: rmh2
  integer :: i
  character :: efile*4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     1d-steady-Euler-equation
!     this is the main subroutine
!     this quasi 1d programe can consider curvature effect
!
!     this subroutine was created by hiroaki watanabe, 2014/11/11
!       and reshaped for f90 and addition of curvature effect
!                                 by shunsuke shigeoka, 2018/5/31
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!---------------------------------------------------------------------------
!   for nitial condition
!---------------------------------------------------------------------------

  ddgd = 0.d0

  !---------------------------------------------------------------------------
  !   for nitial condition
  !---------------------------------------------------------------------------

  write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(6,*) ' '
  write(6,*) ' Input : Calculation Condition'
  write(6,*) ' '
  write(6,*) ' : Debug   on / off'
  write(6,*) ' '
  read(5,'(a)') bdbg

  if (bdbg(1:2) == 'on') write(6,*) ' Debug on in icond !!'



  write(6,*) ' : Grid width [nm] : End step : interval for output'
  write(6,*) ' '
  read(5,*) ddltx, jend, lint
  if (jend == 0)  jend = jmax
  if (lint == 0)  lint = jend + 2

  ddltx = ddltx * 1.d-9
  axx = ddltx * dble(jend-1)

  ijj = int((jend-1) / lint + 1)


  write(6,*) ' : Reaction type'
  write(6,*) '    H2-O2    reaction                     (H2o2)'
  write(6,*) '    H2-O2-He reaction                     (O2He)'
  write(6,*) '    H2-O2-Ar reaction                     (O2Ar)'
  write(6,*) '    H2-N2    reaction                     (H2N2)'
  write(6,*) '    H2-N2-He reaction                     (N2He)'
  write(6,*) '    H2-N2-Ar reaction                     (N2Ar)'
  write(6,*) ' '
  read(5,'(a)')  crct

  brct(1:3) = '#'//crct(1:2)
  brct(4:6) = '-'//crct(3:4)

  write(6,*) ' : Choose ZND model(znd) or WK model(wk)'
  read(5,*)  btyp(1:3)
  if(btyp(1:2) == 'wk') then
  write(6,*) ' : The curvature of shock       '
  read(5,*)  dcuv
  end if
  write(6,*) ' '
  write(6,*) ' model : ',btyp
  write(6,*) ' the curvature of shock : ',dcuv

  write(6,*) ' '
  write(6,*) ' -----   Debug               : ', bdbg(1:3)
  write(6,*) ' -----   Reaction type       : ', crct(1:4)
  write(6,*) ' '
  write(6,*) ' -----   End step            : ', jend
  write(6,*) ' -----   Step interval       : ', lint
  write(6,*) ' -----   Data number         : ', ijj
  write(6,*) ' '
  write(6,*) ' -----   Total grid length   : ', axx,'m'
  write(6,*) ' -----   Grid width          : ', ddltx,'m'
  write(6,*) ' '


!--------------------------------------------------------------------------

  if (lint < jend + 1) then
  write(6,*) ' files are opened !!'
  write(6,*) ' '

  open(10,file=brct(1:6)//'-1', form='formatted')
  open(11,file=brct(1:6)//'-2', form='formatted')
  open(12,file=brct(1:6)//'-3', form='formatted')
  open(13,file=brct(1:6)//'-4', form='formatted')
  open(14,file=brct(1:6)//'-5', form='formatted')
  end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------------
!   allocatable
!---------------------------------------------------------------------------

  allocate(dflw(imax), dcnv(iplus), dmh2(jmax), dtime(jmax))

!----------------------------------------------------------------------------

  allocate(dml(lsp), dplh(7,lsp,3))

  dml(1:lsp) = (/2.016d0,  32.000d0, 1.008d0, 16.000d0, 17.008d0, 18.016d0,     &
  &            33.008d0, 34.016d0, 28.016d0 /)

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
! H2  O
  &  4.19864056d+00, -2.03643410d-03,  6.52040211d-06, -5.48797062d-09,       &
  &  1.77197817d-12, -3.02937267d+04, -8.49032208d-01,                        &
! HO  2
  &  4.30179807d+00, -4.74912097d-03,  2.11582905d-05, -2.42763914d-08,       &
  &  9.29225225d-12,  2.64018485d+02,  3.71666220d+00,                        &
! H2O  2
  &  4.27611269d+00, -5.42822417d-04,  1.67335701d-05, -2.15770813d-08,       &
  &  8.62454363d-12, -1.77025821d+04,  3.43505074d+00,                        &
! N  2
  &  3.298677d+00,    1.4082404d-03,  -3.963222d-06,    5.641515d-09,         &
  & -2.444854d-12,   -1.0208999d+03,   3.950372d+00,                          &

! -- 1000-5000K --
! H  2
  &  3.33727920d+00, -4.94024731d-05,  4.99456778d-07, -1.79566394d-10,       &
  &  2.00255376d-14, -9.50158922d+02, -3.20502331d+00,                        &
! O  2
  &  3.28253784d+00,  1.48308754d-03, -7.57966669d-07,  2.09470555d-10,       &
  & -2.16717794d-14, -1.08845772d+03,  5.45323129d+00,                        &
!   H
  &  2.50000001d+00, -2.30842973d-11,  1.61561948d-14, -4.73515235d-18,       &
  &  4.98197357d-22,  2.54736599d+04, -4.46682914d-01,                        &
!   O
  &  2.56942078d+00, -8.59741137d-05,  4.19484589d-08, -1.00177799d-11,       &
  &  1.22833691d-15,  2.92175791d+04,  4.78433864d+00,                        &
! O  H
  &  3.09288767d+00,  5.48429716d-04,  1.26505228d-07, -8.79461556d-11,       &
  &  1.17412376d-14,  3.61585000d+03,  4.47669610d+00,                        &
! H2  O
  &  3.03399249d+00,  2.17691804d-03, -1.64072518d-07, -9.70419870d-11,       &
  &  1.68200992d-14, -3.00042971d+04,  4.96677010d+00,                        &
! HO  2
  &  4.17228741d+00,  1.88117627d-03, -3.46277286d-07,  1.94657549d-11,       &
  &  1.76256905d-16,  3.10206839d+01,  2.95767672d+00,                        &
! H2O  2
  &  4.16500285d+00,  4.90831694d-03, -1.90139225d-06,  3.71185986d-10,       &
  & -2.87908305d-14, -1.78617877d+04,  2.91615662d+00,                        &
! N  2
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

!---------------------------------------------------------------------------
!   calculation start
!---------------------------------------------------------------------------

!iinit
  call iinit(dflw,dcnv,dmh2,dtime,au0)




  call irhrl(au0)

  do lstp = 1, jend

    if (bdbg(1:2) == 'on') write(6,*) ' in main ',lstp, jend

    ddgd = ddgd + ddltx

    call iintg

    if (mod(lstp,lint) == 0) then
    call ioutp ('data',dmh2(lstp+2))
 !   write(6,*) dcvt, dcnv(lmss)
    endif

  end do

  call iehrl

!---------------------------------------------------------------------------
end subroutine main

!===========================================================================

subroutine ispec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-spec(properties of species)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!---------------------------------------------------------------------------
! block data for chemical reaction and thermal properties
!---------------------------------------------------------------------------



!--------------------------------------------------------------------------

end subroutine ispec

!===========================================================================

subroutine icond

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-cond(condition for calculation)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

!------------------------------------------------------------------------



!------------------------------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end

!========================================================================

subroutine iinit(dflw,dcnv,dmh2,dtime,au0)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-init(initial condition)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  real(8), intent(out) :: dflw(:), dcnv(:),dmh2(:), dtime(:), au0
  real(8), allocatable :: amlf(:), adns(:), acpi(:), ahhi(:)
  real(8) :: at0, ap00, ad0, ar0, ahh0, acp0, amss, ammt, aeng,               &
  &          amc0, agm0, ass0, ap0, atm, atw, ae0
  integer :: i, j, k, m, iflg

!------------------------------------------------------------------------

  if (bdbg(1:2) == 'on') write(6,*) ' in iinit '

!------------------------------------------------------------------------

  allocate(amlf(lsp), adns(lsp), acpi(lsp), ahhi(lsp))

  write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(6,*) ' '
  write(6,*) '         Input : Initial Condition'
  write(6,*) ' '
  write(6,*) ' : mole fraction of H2, O2, H2O, N2'
  write(6,*) ' '
  read(5,*) amlf(l1), amlf(l2), amlf(l6), amlf(l9)
  write(6,*) '   mole fraction of H2      : ',amlf(l1)
  write(6,*) '   mole fraction of O2      : ',amlf(l2)
  write(6,*) '   mole fraction of H2O     : ',amlf(l6)
  write(6,*) '   mole fraction of N2      : ',amlf(l9)
  write(6,*) ' '

  write(6,*) ' : Temprerature [K] : Pressure [atm] : Velocity [m/s]'
  write(6,*) ' '
  write(6,*) ' ** CJ velocity (1atm, 300K) for stoicheiometric'
  write(6,*) '     Reaction Type                CJ Vel.'
  write(6,*) '       2H2-O2                   2841.0 m/s'
  write(6,*) '       2H2-O2-Ar                2268.6 m/s'
  write(6,*) '       2H2-O2-He                3125.4 m/s'
  write(6,*) '       2H2-O2-3.76N2            1971.0 m/s'
  write(6,*) '       2H2-O2-3.76N2-Ar         1856.8 m/s'
  write(6,*) '       2H2-O2-3.76N2-He         2073.7 m/s'
  write(6,*) ' '
  write(6,*) ' ** Ref. AISTJAN **'
  write(6,*) '   URL : http://www.aist.go.jp/RIODB/ChemTherm/aistjan-.html'
  write(6,*) ' '
  read(5,*) at0, ap00, au0

  ap0 = ap00 * 101325.d0

  atm = 0.d0
  atw = 0.d0
  do i = 1, lsp
     atm = atm + amlf(i)
     atw = atw + amlf(i) * dml(i)
  end do

  ar0 = dru * atm / atw
  ad0 = ap0 / ar0 / at0

  iflg = 1 + int(a05 + sign(a05,(at0-dplt)))
  do i = 1, lsp
    adns(i) = ad0 * amlf(i) * dml(i) / atw
  end do

  call ithmc('cp',adns,ad0,at0,ar0,iflg,acp0,acpi)
  call ithmc('hh',adns,ad0,at0,ar0,iflg,ahh0,ahhi)

  ae0 = ad0 * ahh0 - ap0 + 0.5d0*ad0*au0*au0

  agm0 = acp0 / (acp0 - ar0)
  ass0 = sqrt(agm0 * ar0* at0)
  amc0 = au0 / ass0

!  initial mass fraction of reactant H2
  dmh2(1) = adns(1) / ad0

!  half-reaction time
  dtime(1) = 0.d0

  write(6,*) ' :  Initial Condition'
  write(6,*) ' ----   Pressure             : ',ap00,'atm'
  write(6,*) ' ----   Temperature          : ',at0, 'K'
  write(6,*) ' ----   Velocity             : ',au0, 'm/s'
  write(6,*) ' ----   Density              : ',ad0, 'kg/m3'
  write(6,*) ' ----   Gas Constant         : ',ar0, 'J/kg/K'
  write(6,*) ' ----   Total Energy         : ',ae0, 'J/m3'
  write(6,*) ' ----   Total Enthalpy       : ',ahh0,'J/kg'
  write(6,*) ' '
  write(6,*) ' ----   Sound Speed          : ',ass0,'m/s'
  write(6,*) ' ----   Specific Heat Ratio  : ',agm0
  write(6,*) ' ----   Mach Number          : ',amc0
  write(6,*) ' '
  write(6,*) ' ----   Mass Fraction of H2  : ',dmh2(1)
  write(6,*) ' '


!---------------------------------------------------------------------------
!  Conservation Values
!---------------------------------------------------------------------------
!  mass conservation
  amss = ad0  * au0

!  momentum conservation
  ammt = amss * au0 + ap0

!  energy conservation
!  aeng = (ae0+ap0) * au0
  aeng = ahh0 + 0.5d0*au0*au0


  write(6,*) ' :  Conservation Values'
  write(6,*) ' ----   mass                 : ',amss,'kg/m3/s'
  write(6,*) ' ----   momentum             : ',ammt,'Pa'
  write(6,*) ' ----   energy               : ',aeng,'Pa m/s'
  write(6,*) ' '


!--------------------------------------------------------------------------
! for data & flow
!--------------------------------------------------------------------------


  do i = l1, lsp
     dflw(i) = adns(i)
  end do

  dflw(lr) = ad0
  dflw(lu) = au0
  dflw(le) = ae0
  dflw(lp) = ap0
  dflw(lt) = at0
  dflw(lg) = ar0

  dcnv(lmss) = amss
  dcnv(lmmt) = ammt
  dcnv(leng) = aeng

  if (lint < jend + 1)  call ioutp('init',amc0)

!-------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine iinit

!=========================================================================

subroutine irhrl(au0)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-r(Rankine)-h(Hugoniot)-rl(relation)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  real(8), intent(in) :: au0
  real(8), allocatable :: acpi(:), ahhi(:), amsf(:), adns(:),                 &
  &                       paaa(:,:), aaa(:)
  real(8) :: pfnn, pfnd, ats, acnt, ads, aus, aps, aes, agms, asss, amcs,     &
  &          afnn, afnd, ap1, ap2, ap3, ap4, aq1, aq2, aq3, aq4, aq5, arrr,   &
  &          addd, avel, atmp, agm, agm1, agp1, asnd, amac, a1, a2, a3,       &
  &          amss, ammt, aeng, aerr, acp, acps, ahhs
  integer :: i, k, iflg, p

!--------------------------------------------------------------------------
!
!      if (bdbg(1:2) == 'on') write(6,*) ' in irhrl'
!
!--------------------------------------------------------------------------

  allocate(adns(lsp),amsf(lsp),acpi(lsp),ahhi(lsp))

  amss = dcnv(lmss)
  ammt = dcnv(lmmt)
  aeng = dcnv(leng)

!---------------------------------------------------------------------------

  aerr = 1.d-5

!---------------------------------------------------------------------------

  do i = 1, lsp
    adns(i) = dflw(i)
  end do

  addd = dflw(lr)
  avel = dflw(lu)
  atmp = dflw(lt)
  arrr = dflw(lg)

  iflg = 1 + int(a05 + sign(a05,(atmp-dplt)))
  call ithmc('cp',adns,addd,atmp,arrr,iflg,acp,acpi)

  agm = acp / (acp - arrr)
  agm1 = agm - 1.d0
  agp1 = agm + 1.d0

  !acnt = aeng / amss
  acnt = aeng

!---------------------------------------------------------------------------

  write(6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(6,*) ' '
  write(6,*) '     Rankine-Hugoniot Relation is being solved'
  write(6,*) ' '

  do i = 1, lsp
     amsf(i) = dflw(i) / dflw(lr)
  end do

  asnd = sqrt(agm*arrr*atmp)
  amac = avel / asnd
  a1 = 2.d0*agm*amac*amac - agm1
  a2 = agm1*amac*amac + 2.d0
  a3 = agp1*agp1 * amac*amac
  atmp = atmp * a1 * a2 / a3


  write(6,*) ' assumed temperature : ', atmp,'K '
  write(6,*) ' specific heat ratio : ', agm
  write(6,*) ' '


!============================
!      Newton Method
!============================

  allocate(paaa(6,2), aaa(6))

  do p = 1, 10000

    iflg = 1 + int(a05 + sign(a05,(atmp-dplt)))

    ! calculate sum(YRa)i
    do k = 1, 6
      paaa(k,iflg)  = ( amsf( 1) * dplh(k, 1,iflg) / dml( 1)                &
      &              +  amsf( 2) * dplh(k, 2,iflg) / dml( 2)                &
      &              +  amsf( 3) * dplh(k, 3,iflg) / dml( 3)                &
      &              +  amsf( 4) * dplh(k, 4,iflg) / dml( 4)                &
      &              +  amsf( 5) * dplh(k, 5,iflg) / dml( 5)                &
      &              +  amsf( 6) * dplh(k, 6,iflg) / dml( 6)                &
      &              +  amsf( 7) * dplh(k, 7,iflg) / dml( 7)                &
      &              +  amsf( 8) * dplh(k, 8,iflg) / dml( 8)                &
      &              +  amsf( 9) * dplh(k, 9,iflg) / dml( 9))               &
!      &              +  amsf(10) * dplh(ii,10,inum) / dml(10)                &
!      &              +  amsf(11) * dplh(ii,11,inum) / dml(11))               &
      &              *  dru

      aaa(k) = paaa(k,iflg)
    end do

! non-differential form
  ap1 = ammt*ammt - 4.d0*amss*amss*arrr*atmp    ! (p-rho*u^2)^2
  ap2 = ammt + sqrt(ap1)                        ! 2p
  ap3 = 0.5d0 * ap2 / arrr / atmp               ! rho
  ap4 = 0.5d0 * amss * amss / ap3 / ap3         ! 0.5*u^2
  afnn =       atmp * (aaa(1)                                              &
  &          + atmp * (aaa(2) / 2.d0                                       &
  &          + atmp * (aaa(3) / 3.d0                                       &
  &          + atmp * (aaa(4) / 4.d0                                       &
  &          + atmp * (aaa(5) / 5.d0)))))                                  &
  &          + aaa(6) + ap4 - acnt

! differential form of ap4
  aq1 = 2.d0 * arrr * atmp * atmp
  aq2 = 2.d0 * amss * amss * arrr * atmp
  aq3 = sqrt(ap1)
  aq4 = -1.d0 * (aq2/aq3 + ap2) / aq1
  aq5 = amss * amss / ap3 / ap3 / ap3
  afnd =               aaa(1)                                              &
  &          + atmp * (aaa(2)                                              &
  &          + atmp * (aaa(3)                                              &
  &          + atmp * (aaa(4)                                              &
  &          + atmp * (aaa(5)))))                                          &
  &          - aq4 * aq5

! Newton iteration
  atmp = atmp - afnn / afnd
  if (abs(afnn/afnd) < aerr) then
  if (atmp < 0.) then
    write(6,*) ' temperature is negative !!', atmp
    write(6,*) ' calculation stops', p
    stop
  end if
    exit
  end if

  end do


  write(6,*) ' temperature behind shock wave  : ',atmp,'K'
  write(6,*) ' '
  ats = atmp

  acnt = ammt*ammt - 4.d0*amss*amss*arrr*ats
  ads  = 0.5d0 * (ammt + sqrt(acnt)) / arrr / ats
  aus  = amss / ads
  aps  = ads * arrr * ats

  do i = 1, lsp
    adns(i) = amsf(i) * ads
  end do

  iflg = 1 + int(a05 + sign(a05,(ats-dplt)))
  call ithmc('cp',adns,ads,ats,arrr,iflg,acps,acpi)
  call ithmc('hh',adns,ads,ats,arrr,iflg,ahhs,ahhi)

  aes  = ads*ahhs - aps + 0.5d0*ads*aus*aus

  agms = acps / (acps - arrr)
  asss = sqrt(agms * arrr * ats)
  amcs = aus / asss

  dmh2(2) = adns(1) / ads

  if(btyp(1:2) == 'wk')  then
  write(6,*) au0, aus, ddltx, dcuv
  dcvt = ddltx*dcuv*(au0 - aus)
  end if

  write(6,*) ' :  Condition behind shock wave'
  write(6,*) ' '
  write(6,*) ' ----   pressure             : ', aps,'Pa'
  write(6,*) ' ----   temperature          : ', ats,'K'
  write(6,*) ' ----   velocity             : ', aus,'m/s'
  write(6,*) ' ----   denisity             : ', ads,'kg/m3'
  write(6,*) ' '
  write(6,*) ' ----   Mach number          : ', amcs
  write(6,*) ' ----   sound speed          : ', asss,'m/s'
  write(6,*) ' ----   Gas Constant         : ', arrr,'J/kg/K'
  write(6,*) ' '
  write(6,*) ' ----   mass fraction of H2  : ',dmh2(2)
  write(6,*) ' '


!--------------------------------------------------------------------------
! conservation values
!--------------------------------------------------------------------------
!  mass conservation
  amss = ads  * aus

!  momentum conservation
  ammt = amss * aus + aps

!  energy conservation
!  aeng = (aes+aps) * aus
  aeng = ahhs + 0.5d0*aus*aus

  write(6,*) ' :  Conservation Values'
  write(6,*) ' ----   mass                 : ',amss,'kg/m3/s'
  write(6,*) ' ----   momentum             : ',ammt,'Pa'
  write(6,*) ' ----   energy               : ',aeng,'Pa m/s'
  write(6,*) ' '


!--------------------------------------------------------------------------

  do i = 1, lsp
     dflw(i) = adns(i)
  end do

  dflw(lr) = ads
  dflw(lu) = aus
  dflw(le) = aes
  dflw(lt) = ats
  dflw(lp) = aps
  dflw(lg) = arrr

  if (lint < jend + 1) call ioutp('data',amcs)

  deallocate(paaa, aaa)

!---------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine irhrl

!==========================================================================

subroutine iintg

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-intg(integration)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

  real(8), allocatable :: amsf(:), adns(:), ahhi(:), acpi(:),                 &
  &                       pbbb(:,:), acc(:)
  real(8) :: aerr, arr, amss, ammt, aeng, t1, t2
  real(8) :: pfnn, pfnd, ats, acnt, ads, aus, aps, aes, agms, asss, amcs,     &
  &          afnn, afnd, ap1, ap2, ap3, ap4, aq1, aq2, aq3, aq4, aq5, arrr,   &
  &          addd, avel, atmp, agm, agm1, agp1, asnd, amac, a1, a2, a3, aprs, &
  &          axx, a11, a21, a22, a23, aee, adt, acp, ahh, ass, amc
  integer :: i, j, k, m, iflg
!--------------------------------------------------------------------------

!  pbbb(ii,itm) =  (amsf( 1) * dplh(ii, 1,itm) / dml( 1)
!  &             +  amsf( 2) * dplh(ii, 2,itm) / dml( 2)
!  &             +  amsf( 3) * dplh(ii, 3,itm) / dml( 3)
!  &             +  amsf( 4) * dplh(ii, 4,itm) / dml( 4)
!  &             +  amsf( 5) * dplh(ii, 5,itm) / dml( 5)
!  &             +  amsf( 6) * dplh(ii, 6,itm) / dml( 6)
!  &             +  amsf( 7) * dplh(ii, 7,itm) / dml( 7)
!  &             +  amsf( 8) * dplh(ii, 8,itm) / dml( 8)
!  &             +  amsf( 9) * dplh(ii, 9,itm) / dml( 9))
!!     -             +  amsf(10) * dplh(ii,10,itm) / dml(10)
!!     -             +  amsf(11) * dplh(ii,11,itm) / dml(11))
!  &             *  dru

!--------------------------------------------------------------------------
!---------------------------------------------------------------------------
!   for interface to C++
!---------------------------------------------------------------------------

  !interface
  !  subroutine ann_integrator(ann_t,ann_p,ann_m) bind(c,name='prediction')
  !    import
  !    real(c_double), intent(in) :: ann_t, ann_p
  !    real(c_double), intent(inout) :: ann_m(8)
  !  end subroutine ann_integrator
  !end interface
  !interface
  !  subroutine ann_integrator() bind(c,Name='prediction')
  !    import
  !  end subroutine ann_integrator
  !end interface

!--------------------------------------------------------------------------

  if (bdbg(1:2) == 'on') write(6,*) ' in iintg'

!--------------------------------------------------------------------------

  allocate(amsf(lsp),adns(lsp),ahhi(lsp),acpi(lsp),pbbb(6,2),acc(6))

  aerr = 1.d-5

!--------------------------------------------------------------------------
! integrate source term

  select case(btyp(1:2))

  case('zn')

    call isrce (amsf)
    !ann_temp = dflw(lt)
    !ann_pres = dflw(lp)
    !call ann_integrator(ann_temp,ann_pres,amsf)

  case('wk')

!    call cpu_time (t1)
!    do i = 1, 1000
    call isrce(amsf)
!    end do
!    call cpu_time (t2)
!    write(6,*) "cpu time: ", (t2-t1)*1000

  end select

  if (bdbg(1:2) == 'on') write(6,*) ' back from isrce'

    arr =  (amsf( 1) / dml( 1)                                               &
    &    +  amsf( 2) / dml( 2)                                               &
    &    +  amsf( 3) / dml( 3)                                               &
    &    +  amsf( 4) / dml( 4)                                               &
    &    +  amsf( 5) / dml( 5)                                               &
    &    +  amsf( 6) / dml( 6)                                               &
    &    +  amsf( 7) / dml( 7)                                               &
    &    +  amsf( 8) / dml( 8)                                               &
    &    +  amsf( 9) / dml( 9))                                              &
!     -    +  amsf(10) / dml(10)                                              &
!     -    +  amsf(11) / dml(11))                                             &
    &    *  dru

  select case(btyp(1:2))
  case('zn')
    amss = dcnv(lmss)
    ammt = dcnv(lmmt)
    aeng = dcnv(leng)

  case('wk')

!    if (dmc >= 1.d0-epsl) then
!    if (dmc <= 1.d0+epsl) then
!    amss = dcnv(lmss)
!    ammt = dcnv(lmmt)
!    aeng = dcnv(leng)
!    end if
!    end if

!    if (dmc <= 1.d0-epsl) then
!    amss = dcnv(lmss) - dflw(lr)*dcvt
!    ammt = dcnv(lmmt) - dflw(lr)*dflw(lu)*dcvt
!    aeng = dcnv(leng)
!    end if

!    if (dmc >= 1.d0+epsl) then
    amss = dcnv(lmss) - dflw(lr)*dcvt
    ammt = dcnv(lmmt) - dflw(lr)*dflw(lu)*dcvt
    aeng = dcnv(leng)
!    end if

  end select

!  acnt = aeng / amss
  acnt = aeng

  addd = dflw(lr)
  aprs = dflw(lp)
  atmp = dflw(lt)

  iflg = 1 + int(a05 + sign(a05,(atmp-dplt)))
  do k = 1, 6
    acc(k) =        (amsf( 1) * dplh(k, 1,iflg) / dml( 1)                     &
    &             +  amsf( 2) * dplh(k, 2,iflg) / dml( 2)                     &
    &             +  amsf( 3) * dplh(k, 3,iflg) / dml( 3)                     &
    &             +  amsf( 4) * dplh(k, 4,iflg) / dml( 4)                     &
    &             +  amsf( 5) * dplh(k, 5,iflg) / dml( 5)                     &
    &             +  amsf( 6) * dplh(k, 6,iflg) / dml( 6)                     &
    &             +  amsf( 7) * dplh(k, 7,iflg) / dml( 7)                     &
    &             +  amsf( 8) * dplh(k, 8,iflg) / dml( 8)                     &
    &             +  amsf( 9) * dplh(k, 9,iflg) / dml( 9))                    &
!       -             +  amsf(10) * dplh(ii,10,itm) / dml(10)                 &
!       -             +  amsf(11) * dplh(ii,11,itm) / dml(11))                &
    &             *  dru
  end do

  do m = 1, 100

  axx = amss * amss / addd / addd    ! u^2
  a11 = 0.5d0 * axx                  ! 0.5u^2
  a21 = 2.d0  * axx / addd           ! 2u^2/rho
  a22 = ammt / addd / addd           ! (p+rho*u^2)/rho^2
  a23 = axx / addd                   ! u^2/rho

! non-differential form
  afnn =     atmp * (acc(1)                                                   &
  &        + atmp * (acc(2) / 2.d0                                            &
  &        + atmp * (acc(3) / 3.d0                                            &
  &        + atmp * (acc(4) / 4.d0                                            &
  &        + atmp * (acc(5) / 5.d0)))))                                       &
  &        + acc(6) + a11 - acnt

! differential form
  afnd =             acc(1)                                                   &
  &        + atmp * (acc(2)                                                   &
  &        + atmp * (acc(3)                                                   &
  &        + atmp * (acc(4)                                                   &
  &        + atmp * (acc(5)))))
  afnd = afnd * (a21 - a22) / arr - a23

! Newton Method
  addd = addd - afnn / afnd
  if (abs(afnn/afnd) < aerr) then
  if (addd < 0.d0) then
     write(6,*) ' density is negative !!', addd, m
     stop
  end if
  exit
  end if

  end do


  do i = 1, lsp
    adns(i) = amsf(i) * addd
  end do

  avel = amss / addd
  aprs = ammt - amss * avel
  atmp = aprs / arr / addd

  iflg = 1 + int(a05 + sign(a05,(atmp-dplt)))
  call ithmc('cp',adns,addd,atmp,arr,iflg,acp,acpi)
  call ithmc('hh',adns,addd,atmp,arr,iflg,ahh,ahhi)

  agm = acp / (acp - arr)
  ass = sqrt(agm * arr * atmp)
  amc = avel / ass
  dmc = amc

  aee = addd * ahh - aprs + 0.5d0*addd*avel*avel

  do i = 1, lsp
    dflw(i) = adns(i)
  end do

  dflw(lr) = addd
  dflw(lu) = avel
  dflw(le) = aee
  dflw(lp) = aprs
  dflw(lt) = atmp
  dflw(lg) = arr


! mass and momentum conservation for wk model
  select case(btyp(1:2))
  case('wk')
    dcnv(lmss) = amss
    dcnv(lmmt) = ammt

  end select

! mass fraction of reactant H2
  dmh2(lstp+2) = dflw(1) / dflw(lr)

! half-reaction time
  adt = ddltx / avel
  dtime(lstp+1) = dtime(lstp) + adt

!---------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine iintg

!=======================================================================

subroutine isrce(rmsf)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-srce(source term)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  real(8), allocatable :: adns(:), rmsf(:), asoi(:), ahhi(:)
  real(8) :: at, afmt, amss, aeng, afdd, afd, af, addd, adt
  real(8) :: atime, atmp, arho, aflag, art, atmpr, atmpl, aptr, adnsr, arr,   &
  &          assm, ahsm, asum
  real(8) :: ar, au, av, awk1, awk2, awk3
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
  &          afg, adlt, adlt1, adlt2, asgm, atmx, aaa, abb, ajcb, aflgm1

  integer :: i, j, k, itm, m, iflg
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

!--------------------------------------------------------------------------
!
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
!        23:  h2          + M   = h2         + M
!        24:  o    + o    + M   = o2         + M
!        25:  o    + h    + M   = oh         + M
!
!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////
!
!
!--------------------------------------------------------------------------
!     if (bdbg(1:2) == 'on') write(lwrt,*) ' in isrce '
!--------------------------------------------------------------------------

! allocate
  allocate(adns(lsp),asoi(lsp),ahhi(lsp))

  do i = 1, lsp
    rmsf(i) = dflw(i) / dflw(lr)
    adns(i) = dflw(i)
  end do

  addd = dflw(lr)
  atmp = dflw(lt)
  arr  = dflw(lg)

!--------------------------------------------------------------------------

  art = 1.d0 / (arc*atmp)

!    [Ea/R]  activation energy

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

!    [kf] forward reaction rate constant
!      Arrhenius's form : kf=AT**n exp(-Ea/RT)

      akf01  = acf( 1)*atmp**ant( 1)*exp(ae01)
      akf02h = acf( 2)*atmp**ant( 2)*exp(ae02h)
      akf02l = acf( 3)*atmp**ant( 3)*exp(ae02l)
      akf03h = acf( 4)*atmp**ant( 4)*exp(ae03h)
      akf03l = acf( 5)*atmp**ant( 5)*exp(ae03l)
      akf04h = acf( 6)*atmp**ant( 6)*exp(ae04h)
      akf04l = acf( 7)*atmp**ant( 7)*exp(ae04l)
      akf05h = acf( 8)*atmp**ant( 8)*exp(ae05h)
      akf05l = acf( 9)*atmp**ant( 9)*exp(ae05l)
      akf06a = acf(10)*atmp**ant(10)*exp(ae06a)
      akf06b = acf(11)*atmp**ant(11)*exp(ae06b)
      akf07  = acf(12)*atmp**ant(12)*exp(ae07)
      akf08a = acf(13)*atmp**ant(13)*exp(ae08a)
      akf08b = acf(14)*atmp**ant(14)*exp(ae08b)
      akf09  = acf(15)*atmp**ant(15)*exp(ae09)
      akf10  = acf(16)*atmp**ant(16)*exp(ae10)
      akf11  = acf(17)*atmp**ant(17)*exp(ae11)
      akf12a = acf(18)*atmp**ant(18)*exp(ae12a)
      akf12b = acf(19)*atmp**ant(19)*exp(ae12b)
      akf13  = acf(20)*atmp**ant(20)*exp(ae13)
      akf14  = acf(21)*atmp**ant(21)*exp(ae14)
      akf15  = acf(22)*atmp**ant(22)*exp(ae15)
      akf16  = acf(23)*atmp**ant(23)*exp(ae16)
      akf17  = acf(24)*atmp**ant(24)*exp(ae17)
      akf18  = acf(25)*atmp**ant(25)*exp(ae18)
      akf19  = acf(26)*atmp**ant(26)*exp(ae19)
      akf20  = acf(27)*atmp**ant(27)*exp(ae20)
      akf21  = acf(28)*atmp**ant(28)*exp(ae21)
      akf22  = acf(29)*atmp**ant(29)*exp(ae22)
      akf23  = acf(30)*atmp**ant(30)*exp(ae23)
      akf24  = acf(31)*atmp**ant(31)*exp(ae24)
      akf25  = acf(32)*atmp**ant(32)*exp(ae25)

!-----------------------------------------------------------------------
! mole density
!-----------------------------------------------------------------------

      ah2   = adns(1) / dml(1)
      ao2   = adns(2) / dml(2)
      ah    = adns(3) / dml(3)
      ao    = adns(4) / dml(4)
      aoh   = adns(5) / dml(5)
      ah2o  = adns(6) / dml(6)
      aho2  = adns(7) / dml(7)
      ah2o2 = adns(8) / dml(8)
      an2   = adns(9) / dml(9)

! third body

      am04  = ah2*1.5  +ao2*0.0  +ah  +ao  +aoh  +ah2o*0.0  +aho2             &
      &       + ah2o2    +an2
      am05  = ah2      +ao2      +ah  +ao  +aoh  +ah2o*9.0  +aho2             &
      &       + ah2o2    +an2*1.5
      am09  = ah2*3.0  +ao2*1.5  +ah  +ao  +aoh  +ah2o*0.0  +aho2          &
      &       + ah2o2    +an2*2.0
      am21  = ah2*0.0  +ao2*0.0  +ah  +ao  +aoh  +ah2o*14.4 +aho2             &
      &       + ah2o2    +an2*0.0
      am23  = ao2      +an2
      am24  = ah2*2.5  +ao2      +ah  +ao  +aoh  +ah2o*12   +aho2             &
      &       + ah2o2    +an2
      am25  = ah2*2.5  +ao2      +ah  +ao  +aoh  +ah2o*12   +aho2             &
      &       + ah2o2    +an2

!-----------------------------------------------------------
!  the reaction rate constant of #02
!-----------------------------------------------------------

      atrp = akf02l*ah2o*1.0e-3/akf02h

!      atrc = -0.4 - 0.67*log10(0.8)
!      atrn = 0.75 - 1.27*log10(0.8)

      atrpl= log10(atrp+1.0e-30)
      atra = (atrpl+atrc8)/(atrn8-0.14*(atrpl+atrc8))
      atrb = 1.0/(1.0 + atra*atra)
      atrf = 0.8**atrb

      akf02 = akf02h*atrp/(1.0+atrp)*atrf
!----------------------------------------------------------
!  the reaction rate constant of #03
!----------------------------------------------------------

      atrp = akf03l*ao2*1.0e-3/akf03h

!      atrc = -0.4 - 0.67*log10(0.7)
!      atrn = 0.75 - 1.27*log10(0.7)

      atrpl= log10(atrp+1.0e-30)
      atra = (atrpl+atrc7)/(atrn7-0.14*(atrpl+atrc7))
      atrb = 1.0/(1.0 + atra*atra)
      atrf = 0.7**atrb

      akf03 = akf03h*atrp/(1.0+atrp)*atrf

!----------------------------------------------------------
!  the reaction rate constant of #04
!----------------------------------------------------------

      atrp = akf04l*am04*1.0e-3/akf04h

!      atrc = -0.4 - 0.67*log10(0.7)
!      atrn = 0.75 - 1.27*log10(0.7)

      atrpl= log10(atrp+1.0e-30)
      atra = (atrpl+atrc7)/(atrn7-0.14*(atrpl+atrc7))
      atrb = 1.0/(1.0 + atra*atra)
      atrf = 0.7**atrb

      akf04 = akf04h*atrp/(1.0+atrp)*atrf

!----------------------------------------------------------
!  the reaction rate constant of #05
!----------------------------------------------------------

      atrp = akf05l*am05*1.0e-3/akf05h

!      atrc = -0.4 - 0.67*log10(0.8)
!      atrn = 0.75 - 1.27*log10(0.8)
!
!      atrpl= log10(atrp)
!      atra = (atrpl+atrc8)/(atrn8-0.14*(atrpl+atrc8))
!      atrb = 1.0/(1.0 + atra*atra)
!      atrf = 0.7**atrb
!
!      akf05 = akf05h*atrp/(1.0+atrp)*atrf
      akf05 = akf05h*atrp/(1.0+atrp)

!-----------------------------------------------------------
!  the reaction rate constant of #06 & #08 & #12
!-----------------------------------------------------------

      akf06 = akf06a + akf06b
      akf08 = akf08a + akf08b
      akf12 = akf12a + akf12b

!------------------------------------------------------------
!    [kb] backward reaction rate constant
!------------------------------------------------------------

      iflg = 1 + int(a05 + sign(a05,(atmp-dplt)))

! entropy

      call ithmc('ss',adns,addd,atmp,arr,iflg,assm,asoi)

      as1 = asoi(1)
      as2 = asoi(2)
      as3 = asoi(3)
      as4 = asoi(4)
      as5 = asoi(5)
      as6 = asoi(6)
      as7 = asoi(7)
      as8 = asoi(8)

! enthalpy

      call ithmc('hh',adns,addd,atmp,arr,iflg,ahsm,ahhi)

      acon = 1.d0 / (dru*atmp)
      ah1 = ahhi(1)*dml(1)*acon
      ah2 = ahhi(2)*dml(2)*acon
      ah3 = ahhi(3)*dml(3)*acon
      ah4 = ahhi(4)*dml(4)*acon
      ah5 = ahhi(5)*dml(5)*acon
      ah6 = ahhi(6)*dml(6)*acon
      ah7 = ahhi(7)*dml(7)*acon
      ah8 = ahhi(8)*dml(8)*acon



      apt = arp*atmp
!..01
      ahh01 = ah5+ah4-ah3-ah2
      ag    = (as5+as4-as3-as2)-ahh01
      akb01 = akf01/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))))
!..02
      ahh02 = ah7    -ah3-ah2
      ag    = (as7    -as3-as2)-ahh02
      akb02 = akf02/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       * apt))
!..03
      ahh03 = ah7    -ah3-ah2
      ag    = (as7    -as3-as2)-ahh03
      akb03 = akf03/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       * apt))
!..04
      ahh04 = ah7    -ah3-ah2
      ag    = (as7   -as3-as2)-ahh04
      akb04 = akf04/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       * apt))
!..05
      ahh05 = ah5+ah5    -ah8
      ag    = (as5+as5    -as8)-ahh05
      akb05 = akf05/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       / apt))
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
      &       / apt))
!..10
      ahh10 = ah5+ah3    -ah6
      ag    = (as5+as3    -as6)-ahh10
      akb10 = akf10/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       / apt))
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
      &       / apt))
!..22
      ahh22 = ah3+ah3    -ah1
      ag    = (as3+as3    -as1)-ahh22
      akb22 = akf22/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       / apt))
!..23
      ahh23 = ah3+ah3    -ah1
      ag    = (as3+as3    -as1)-ahh23
      akb23 = akf23/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       / apt))
!..24
      ahh24 = ah2    -ah4-ah4
      ag    = (as2    -as4-as4)-ahh24
      akb24 = akf24/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       * apt))
!..25
      ahh25 = ah5    -ah4-ah3
      ag    = (as5    -as4-as3)-ahh25
      akb25 = akf25/min(akmx,max(akmn,exp(max(aemn,min(ag,aemx)))             &
      &       * apt))

!------------------------------------------------------------

      akf01 = akf01 * 1.e-3
      akf02 = akf02 * 1.e-3
      akf03 = akf03 * 1.e-3
      akf04 = akf04 * 1.e-3
!         akf05 = akf05 * 1.e-3
      akf06 = akf06 * 1.e-3
      akf07 = akf07 * 1.e-3
      akf08 = akf08 * 1.e-3
      akf09 = akf09 * 1.e-3
      akf10 = akf10 * 1.e-3
      akf11 = akf11 * 1.e-3
      akf12 = akf12 * 1.e-3
      akf13 = akf13 * 1.e-3
      akf14 = akf14 * 1.e-3
      akf15 = akf15 * 1.e-3
      akf16 = akf16 * 1.e-3
      akf17 = akf17 * 1.e-3
      akf18 = akf18 * 1.e-3
      akf19 = akf19 * 1.e-3
      akf20 = akf20 * 1.e-3
      akf21 = akf21 * 1.e-3
      akf22 = akf22 * 1.e-3
      akf23 = akf23 * 1.e-3
      akf24 = akf24 * 1.e-6
      akf25 = akf25 * 1.e-6

      akb01 = akb01 * 1.e-3
!      akb02 = akb02 * 1.e-3
!      akb03 = akb03 * 1.e-3
!      akb04 = akb04 * 1.e-3
      akb05 = akb05 * 1.e-3
      akb06 = akb06 * 1.e-3
      akb07 = akb07 * 1.e-3
      akb08 = akb08 * 1.e-3
      akb09 = akb09 * 1.e-6
      akb10 = akb10 * 1.e-6
      akb11 = akb11 * 1.e-3
      akb12 = akb12 * 1.e-3
      akb13 = akb13 * 1.e-3
      akb14 = akb14 * 1.e-3
      akb15 = akb15 * 1.e-3
      akb16 = akb16 * 1.e-3
      akb17 = akb17 * 1.e-3
      akb18 = akb18 * 1.e-3
      akb19 = akb19 * 1.e-3
      akb20 = akb20 * 1.e-3
      akb21 = akb21 * 1.e-6
      akb22 = akb22 * 1.e-6
      akb23 = akb23 * 1.e-6
      akb24 = akb24 * 1.e-3
      akb25 = akb25 * 1.e-3

!--------------------------------------------------------------------------
      ah2   = adns(1) / dml(1)
!      ao2   = adns(2) / dml(2)
!      ah    = adns(3) / dml(3)
!      ao    = adns(4) / dml(4)
!      aoh   = adns(5) / dml(5)
!      ah2o  = adns(6) / dml(6)
!      aho2  = adns(7) / dml(7)
!      ah2o2 = adns(8) / dml(8)
!      an2   = adns(9) / dml(9)
!
!--------------------------------------------------------------------------
!
!    kf[ ][ ] : rection rate constant * mole density
!
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

!   kb[ ][ ] : rection rate constant * mole density

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


!    RR = kf[ ][ ] - kb[ ][ ]

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

!   RR of each species

      awh2   = dml(1) * (-ar12 -ar13 +ar16 +ar18 -ar21 -ar22 -ar23)
      awo2   = dml(2) * (-ar01 -ar02 -ar03 -ar04 +ar07 +ar08 +ar16            &
      &                  +ar17 +ar24)
      awh    = dml(3) * (-ar01 -ar02 -ar03 -ar04 +ar09 +ar10 +ar12            &
      &                  +ar13 -ar14 -ar15 -ar16 -ar18 -ar19 +ar21            &
      &                  +ar21 +ar22 +ar22 +ar23 +ar23 -ar25)
      awo    = dml(4) * (+ar01 +ar11 -ar12 +ar15 -ar17 -ar20 -ar24            &
      &                  -ar24 -ar25)
      awoh   = dml(5) * (+ar01 +ar05 +ar05 -ar06 -ar07 +ar09 +ar10            &
      &                  -ar11 -ar11 +ar12 -ar13 +ar14 +ar14 +ar17            &
      &                  +ar19 +ar20 +ar25)
      awh2o  = dml(6) * (+ar06 +ar07 -ar09 -ar10 +ar11 +ar13 +ar15            &
      &                  +ar19)
      awho2  = dml(7) * (+ar02 +ar03 +ar04 +ar06 -ar07 -ar08 -ar08            &
      &                  -ar14 -ar15 -ar16 -ar17 +ar18 +ar20)
      awh2o2 = dml(8) * (-ar05 -ar06 +ar08 -ar18 -ar19 -ar20)


!--------------------------------------------------------------------------
! integration for each species
!--------------------------------------------------------------------------

      amss = dcnv(lmss)

      select case(btyp(1:2))
      case('zn')
        rmsf(1) = rmsf(1) + ddltx * awh2   / amss
        rmsf(2) = rmsf(2) + ddltx * awo2   / amss
        rmsf(3) = rmsf(3) + ddltx * awh    / amss
        rmsf(4) = rmsf(4) + ddltx * awo    / amss
        rmsf(5) = rmsf(5) + ddltx * awoh   / amss
        rmsf(6) = rmsf(6) + ddltx * awh2o  / amss
        rmsf(7) = rmsf(7) + ddltx * awho2  / amss
        rmsf(8) = rmsf(8) + ddltx * awh2o2 / amss

      case('wk')
        rmsf(1) = rmsf(1) + (awh2   / amss - rmsf(1)*dcvt)*ddltx
        rmsf(2) = rmsf(2) + (awo2   / amss - rmsf(2)*dcvt)*ddltx
        rmsf(3) = rmsf(3) + (awh    / amss - rmsf(3)*dcvt)*ddltx
        rmsf(4) = rmsf(4) + (awo    / amss - rmsf(4)*dcvt)*ddltx
        rmsf(5) = rmsf(5) + (awoh   / amss - rmsf(5)*dcvt)*ddltx
        rmsf(6) = rmsf(6) + (awh2o  / amss - rmsf(6)*dcvt)*ddltx
        rmsf(7) = rmsf(7) + (awho2  / amss - rmsf(7)*dcvt)*ddltx
        rmsf(8) = rmsf(8) + (awh2o2 / amss - rmsf(8)*dcvt)*ddltx

      end select

      asum = 0.
      do i = 1, lsp
        asum = asum + rmsf(i)
        if (rmsf(i) < 0.d0) then
          write(6,*) ' mass fraction is negative!!',rmsf(i),amss,awo2
          stop
        end if
      end do

      do i = 1, lsp
        rmsf(i) = rmsf(i) / asum
      end do

!-------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine isrce

!==========================================================================

subroutine iehrl

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-e(estimate)-hrl(half-reaction length L1/2)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  real(8), allocatable :: adns(:), acpi(:), ahhi(:)
  real(8) :: addd, atmp, arrr, avel, aprs, acp, agm, ass,              &
  &          amc, ah12, ax12, at12, aeng, amss, ammt, ahh
  integer :: j, i, iflg, ij12

!--------------------------------------------------------------------------

  allocate(adns(lsp), acpi(lsp), ahhi(lsp))

  addd = dflw(lr)
  atmp = dflw(lt)
  arrr = dflw(lg)
  avel = dflw(lu)
  aprs = dflw(lp)

  do i = 1, lsp
     adns(i) = dflw(i)
  end do

  iflg = 1 + int(a05 + sign(a05,(atmp-dplt)))
  call ithmc('cp',adns,addd,atmp,arrr,iflg,acp,acpi)
  call ithmc('hh',adns,addd,atmp,arrr,iflg,ahh,ahhi)


  agm = acp / (acp - arrr)
  ass = sqrt(agm * arrr * atmp)
  amc = avel / ass

  ah12 = 0.5d0 * (dmh2(1) + dmh2(jend+2))
  do j = 1, jend+2
     if (dmh2(j) < ah12) then
        ij12 = j
        go to 21
     end if
  end do

21   continue

  ax12 = dble(ij12-1) * ddltx
  ax12 = ax12 * 1.d6

  at12 = dtime(ij12 - 1)* 1.d9

  write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(6,*) ' '
  write(6,*) ' L1/2 condition'
  write(6,*) '    Half reaction length L1/2      : ', ax12,'um'
  write(6,*) '    Half reaction time t1/2        : ', at12,'nsec'
  write(6,*) '    Mass fraction of H2 at front   : ', dmh2(1)
  write(6,*) '    Mass fraction of H2 at end     : ', dmh2(jend+2)
  write(6,*) '    Mass fraction of H2 at L1/2    : ', ah12
  write(6,*) ' '
  write(6,*) ' '
  write(6,*) ' at Detonation End'
  write(6,*) '       Mach number              : ', amc
  write(6,*) '       Temperature              : ', dflw(lt),'K'
  write(6,*) '       Pressure                 : ', dflw(lp),'Pa'
  write(6,*) '       Gas Constant             : ', dflw(lg),'J/kg/K'
  write(6,*) ' '

!......  conservation values

  aeng = addd * ahh - aprs + 0.5d0*addd*avel*avel

  amss = addd * avel
  ammt = aprs + amss*avel
  aeng = (aeng + aprs) * avel

  write(6,*) ' :  Conservation Values'
  write(6,*) ' ----   mass                 : ',amss,'kg/m3/s'
  write(6,*) ' ----   momentum             : ',ammt,'Pa'
  write(6,*) ' ----   energy               : ',aeng,'Pa m/s'
  write(6,*) ' '
  write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(6,*) ' '

!--------------------------------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine iehrl

!==========================================================================

subroutine ithmc (efile,rdns,rddd,rtmp,rrr,kflg,rfsm,rfnc)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-thm(thermal)-c(condition)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  real(8), intent(in) :: rdns(:)
  real(8), intent(in) :: rddd, rtmp, rrr
  real(8), allocatable, intent(out) :: rfnc(:)
  integer, intent(in) :: kflg
  character, intent(in) :: efile*2
  real(8), intent(out) :: rfsm
  integer :: i

!--------------------------------------------------------------------------

  if (bdbg(1:2) == 'on') write(6,*) 'in ithmc ', efile

!--------------------------------------------------------------------------

  allocate(rfnc(lsp))

!--------------------------------------------------------------------------
  if (efile(1:2) == 'hh') then
!--------------------------------------------------------------------------

  rfsm = 0.d0
  do i = 1, lsp
    rfnc(i) = (           dplh(1,i,kflg)                                     &
    &            +  rtmp*(dplh(2,i,kflg) / 2.d0                              &
    &            +  rtmp*(dplh(3,i,kflg) / 3.d0                              &
    &            +  rtmp*(dplh(4,i,kflg) / 4.d0                              &
    &            +  rtmp*(dplh(5,i,kflg) / 5.d0))))                          &
    &            +        dplh(6,i,kflg) / rtmp)                             &
    &            *  rtmp* dru / dml(i)
    rfsm = rfsm + rfnc(i) * rdns(i) / rddd
  end do

!--------------------------------------------------------------------------
  elseif (efile(1:2) == 'cp') then
!--------------------------------------------------------------------------

  rfsm = 0.d0
  do i = 1, lsp
    rfnc(i) = (           dplh(1,i,kflg)                                     &
    &            +  rtmp*(dplh(2,i,kflg)                                     &
    &            +  rtmp*(dplh(3,i,kflg)                                     &
    &            +  rtmp*(dplh(4,i,kflg)                                     &
    &            +  rtmp*(dplh(5,i,kflg))))))                                &
    &            *  dru / dml(i)
    rfsm = rfsm + rfnc(i) * rdns(i) / rddd
  end do

!--------------------------------------------------------------------------
  elseif (efile(1:2) == 'ss') then
!--------------------------------------------------------------------------

  rfsm = 0.d0
  do i = 1, lsp
    rfnc(i) =             dplh(1,i,kflg) * log(rtmp)                        &
    &            +  rtmp*(dplh(2,i,kflg)                                    &
    &            +  rtmp*(dplh(3,i,kflg) / 2.d0                             &
    &            +  rtmp*(dplh(4,i,kflg) / 3.d0                             &
    &            +  rtmp*(dplh(5,i,kflg) / 4.d0))))                         &
    &            +        dplh(7,i,kflg)
    rfsm = rfsm + rfnc(i) * rdns(i) / rddd
  end do

!--------------------------------------------------------------------------
  else
!--------------------------------------------------------------------------

   write(6,*) ' error of argument  ', efile
   write(6,*) ' '
   stop

!--------------------------------------------------------------------------
  end if
!--------------------------------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine ithmc

!==========================================================================

subroutine ioutp(efile,rmh2)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     i(subroutine)-outp(output calculation data)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

end subroutine ioutp

!==========================================================================

end module steadystanford
