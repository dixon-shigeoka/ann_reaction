module mod_parameter
!-----------------------------------------------------------------------
!  this program was created by akiko matsuo,                2011/06/24
!                  modified by Ryota Murakami for g2d-code, 2014/06/26
!                  rewrited by Ryo Mikoshiba for nc2d-code, 2015/03/10
!
!        notation:
!        +------------------------------------+
!        |           | local | args. | module |
!        |-----------+-------+-------+--------|
!        | integer   |   i   |   k   |   l    |
!        | real(*8)  |   a   |   r   |   d    |
!        | (real*4)  |   w   |   f   |   s    |
!        | character |   c   |   e   |   b    |
!        +------------------------------------+
!-----------------------------------------------------------------------
 
  implicit none

! Navier-Stokes Equation
  integer, parameter :: ldim = 2 ! dimension
  integer, parameter :: lfl  = 4 ! number of governing equation
  integer, save      :: imax     ! imax = lfl + lsp

  integer, parameter :: lr = 1   ! density 
  integer, parameter :: lu = 2   ! x direction momentum
  integer, parameter :: lv = 3   ! y direction momentum
  integer, parameter :: le = 4   ! internal energy
  integer, parameter :: lp = 4   ! pressure
  integer, parameter :: l1 = 5   ! species 1
  integer, parameter :: l2 = 6   ! species 2
  integer, parameter :: l3 = 7   ! species 3
  integer, parameter :: l4 = 8   ! species 4
  integer, parameter :: l5 = 9   ! species 5
  integer, parameter :: l6 = 10  ! species 5
  integer, parameter :: l7 = 11  ! species 5
  integer, parameter :: l8 = 12  ! species 5
  integer, parameter :: l9 = 13  ! species 5
  integer, parameter :: l10= 14  

! for species
  integer, parameter :: ls1 = 1
  integer, parameter :: ls2 = 2
  integer, parameter :: ls3 = 3
  integer, parameter :: ls4 = 4
  integer, parameter :: ls5 = 5
  integer, parameter :: ls6 = 6
  integer, parameter :: ls7 = 7
  integer, parameter :: ls8 = 8
  integer, parameter :: ls9 = 9
  integer, parameter :: ls10= 10
! for dppd
  integer, save      :: lpp
  !integer, parameter :: lpp = lsp + 2
  integer, parameter :: lpr = 1
  integer, parameter :: lpe = 2
  integer, parameter :: lp1 = 3
  integer, parameter :: lp2 = 4
  integer, parameter :: lp3 = 5
  integer, parameter :: lp4 = 6
  integer, parameter :: lp5 = 7
  integer, parameter :: lp6 = 8
  integer, parameter :: lp7 = 9
  integer, parameter :: lp8 = 10
  integer, parameter :: lp9 = 11
  integer, parameter :: lp10= 12

! for imrhs
  integer, save      :: iimax 
  integer, save      :: lm        ! enthalpy (shima's scheme)
  integer, save      :: lc        ! soundspeed (shima's scheme)
  !integer, parameter :: iimax = imax + 2
  !integer, parameter :: lm    = lfl + lsp + 1   ! enthalpy (shima's scheme)
  !integer, parameter :: lc    = lfl + lsp + 2   ! soundspeed (shima's scheme)

  integer, parameter :: lbcmx = 10 ! upper limit of boundary conditions
  integer, parameter :: llnmx = 10 ! upper limit of line history

! 2 equation turbulence model
  integer, parameter :: itmax = 2 ! number of turbulence model equation
  integer, parameter :: lkt   = 1 ! turbulent energy
  integer, parameter :: let   = 2 ! dissipation factor

! Physical constant
  real(8), parameter :: dprt = 0.9d0       ! turbulent prandtl number
  real(8), parameter :: dsct = 0.9d0       ! turbulent schmidt number
  real(8), parameter :: dru  = 8314.4d0    ! universal gas constant[J/kmol-K]
  real(8), parameter :: dplt = 1000.d0     ! for thermoprop

! for reference values
  real(8), parameter :: aeps2= 1.d-6
  real(8), parameter :: ainf = 1.d30
  real(8), parameter :: api  = acos(-1.d0) ! circular constant, pi
  real(8), parameter :: a05  = 0.5d0
  real(8), parameter :: a12  = 1.d0/2.d0
  real(8), parameter :: a13  = 1.d0/3.d0
  real(8), parameter :: a14  = 1.d0/4.d0
  real(8), parameter :: a15  = 1.d0/5.d0

! for inog 
  integer, parameter :: zmax = 2

! for file IO
  integer, parameter :: iioflw= 11 !#flow number
  integer, parameter :: iiogrd= 20 !#grid number
  integer, parameter :: iicom1= 86 !cannnot use 5 in MPI calculation
  integer, parameter :: iiocel= 88 ! mod_cell IO
  integer, parameter :: iioshk= 89 ! mod_shock IO

!#ifdef para_mpi
!  integer, parameter :: iiompi = 80
!#endif
!#ifdef timer
!  integer, parameter :: iiotim = 95
!  integer, parameter :: timer_trans = 1
!  integer, parameter :: timer_rhs   = 2
!  integer, parameter :: timer_react = 3
!  integer, parameter :: timer_therm = 4
!  integer, parameter :: timer_IO    = 5
!  integer, parameter :: timer_num   = 5
!  character(50) :: timer_list='transport, RHS, reaction, therm, IO'
!#endif

end module mod_parameter
