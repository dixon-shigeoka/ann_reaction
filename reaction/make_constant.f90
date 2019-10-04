subroutine make_constant(dtmp,dprs,aYi,totaldens,aeng)

  use mod_globals
  use mod_parameter, only : imax, ldim, lr, lu, lv, le, l1, ls1, lp1, l2, ls2,     &
       &                         lp2, l3, ls3, lp3, l4, ls4, lp4, l5, ls5, lp5, l6,&
       &                         ls6, lp6, l7, ls7, lp7, l8, ls8, lp8, l9, ls9,    &
       &                         lp9, lfl, lpr, lpe, lpp, a05, a12, a13, a14,      &
       &                         a15, dru, dplt
  use mod_function
  implicit none
  integer :: i, itm
  real(8), intent(in) :: dtmp, dprs, aYi(9)
  real(8), intent(out) :: totaldens, aeng
  real(8), allocatable :: adns(:), afm(:), ab(:), dml(:), dmlr(:), dmlru(:), dplh(:,:,:), atmw(:)
  real(8) :: at, afmt, afdd, afd, af, addd, adt, druo, phti, ahti, arho, atmpr, ah
  real(8), parameter :: aerr = 1.d-12
  integer, parameter :: itr = 10 ! number of iteration in Newton method
  integer, parameter :: lsp = 9

  ! allocate
  allocate(adns(lsp),afm(lsp),ab(6), dml(lsp), dmlr(lsp), dmlru(lsp), dplh(7,lsp,3), atmw(lsp))


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


  ah = 0d0
  atmpr = 1/dtmp
  at = dtmp
  itm = 1 + int(a05 + sign(a05,(at - dplt)))
  atmw(:) = aYi(:)*dmlr(:)      ! amlf(:)/atw
  druo    = dru*(sum(atmw(:)))
  totaldens = dprs*atmpr/druo

  adns(:) = aYi(:)*totaldens
  arho    = sum(adns(1:lsp-1))

  do i = 1, lsp
    call calc_phti(i,itm,dtmp,phti,dmlr,dplh)
    ahti = phti
    ah   = ah  + adns(i)*ahti/arho
  end do

  aeng = dprs - ah*arho

end subroutine make_constant
