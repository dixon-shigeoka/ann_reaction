module mod_globals

  use mod_parameter, only : ldim, lbcmx

  implicit none
  character(2), save :: baod, bsa2
  character(3), save :: bdbg
  character(4), save :: bgcd
  character(4), save :: bcsp, brct, bpyg, bbpy, btic, btrn, bpot
  character(4), save :: bshm, blmt, btst, btin, btti, bflw, bwtm, badt, bini
  character(4), save :: bvis, bvisj, bbvs, bbcd(ldim*2,lbcmx), busr, bthm,    &
    &                   bdif, bltp, btcd(ldim*2,lbcmx), bwal
  character(4), save :: bstp, bsc2, blm2, bcct
  character(4), save :: bofm
  character(4), save :: bflc
  character(4), save :: bscd ! start_condition defined in iscnd
  character(100), save :: bgfn

  integer, save :: lbcd(ldim*2,lbcmx)
  integer, save :: jmax, kmax, ljbgn, ljend, lkbgn, lkend, ljnum, lknum, &
  &                ljbgn_dif, ljend_dif, lkbgn_dif, lkend_dif,           &
  &                ljbgn_int, ljend_int, lkbgn_int, lkend_int,           &
  &                ljbgn_com, ljend_com, lkbgn_com, lkend_com,           &
  &                ljbgn_swp, ljend_swp, lkbgn_swp, lkend_swp,           &
  &                jmax_file, kmax_file, linj, link
  !integer, save :: nmax
  ! bgn_dif : for difference method
  ! bgn_int : for interpolation      ex) MUSCL, WENO

  integer, save :: ljbgnj, ljendj, lkbgnj, lkendj ! for jet mapping
  integer, save :: lbgn, lend, lerr, lmvi, llne
  integer, save :: lrk, lrkmx ! for Runge-Kutta
  integer, save :: lwnb       ! number of wall block
  integer, save :: lfsa       ! for Fortified Solution Algorithm (FSA)
  integer, save :: lprf       ! for thermo dynamic coefficient
  integer, save :: lshp, led  ! for ishop on/off
  integer, save :: lrct       ! for reaction
  integer, save :: linp       ! for avoiding iinit2 in inog
  integer, save :: lmstp      ! for MUSTA
  integer, save :: lsub ! for subiteration
  integer, save :: lrcmx ! max step to record
  integer, save :: lstage ! advection & diffusion or chemcal source term stage

  real(8), save :: daoa, dre, dcfl, dxt, dyt, dercn, dmin, dtin, dpin  &
  &               ,drnd, dund, dend, dmua, dkpa, dfsn
  real(8), save :: ddlt, dlud ! Entropy operator, LU-ADI dissipation
  ! for Boundary condition
  real(8), save :: dprso, dprss, dtins, dps
  real(8), save :: dmxt
  real(8), save :: dwtmp
  ! for time interval #flow.
  real(8), save :: ddfl, dflt
  ! for subiteration
  real(8) :: derr1
  ! for calculation in a rotating frame of reference
  real(8),save :: domg ! angular velocity for rotating frame of reference

  integer, allocatable, save :: lwbc(:,:,:) ! wall block position
  integer, allocatable, save :: lswi(:,:)   ! for HLL-HLLC
  real(8), allocatable, save :: dpsd(:,:)   ! for SD-SLAU
  real(8), allocatable, save :: dchi(:,:), dfsa(:,:)   ! for FSA and MUSCL
  real(8), allocatable, save :: dgrd(:,:,:),dmtx(:,:,:),djcb(:,:),djcr(:,:),  &
  &                             dhlf(:,:,:,:), dmxsq(:,:,:,:)
  real(8), allocatable, save :: dbmap(:,:,:)
  real(8), allocatable, save :: dwwo(:,:,:)
  real(8), allocatable, save :: dgzp(:,:,:) ! for WENO-Z+

  integer, save :: lread, lwrt

  ! for ablation
  real(8), allocatable, save :: dttm(:,:)
  integer, save :: nmx, lcpl
  real(8), save :: dler

  !for Detonation code reaction_jachimowski
  real(8), allocatable, save :: drst(:,:)

  !mpi
  integer, save :: ierr_mpi, myrank_mpi, nprocs_mpi

  ! wall temperature
  real(8), allocatable, save :: dwtm(:,:)

!#ifdef counter
!  integer, allocatable, save :: lmts(:,:)
!#endif

end module mod_globals
