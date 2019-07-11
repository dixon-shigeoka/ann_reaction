module mod_function

  use mod_parameter

  contains

!=======================================================================

  subroutine calc_phti(i,itm,at,phti,dmlr,dplh)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none
  integer, intent(in) :: i, itm
  real(8), intent(in) :: at, dmlr(:), dplh(:,:,:)
  real(8), intent(out) :: phti
  real(8), parameter :: a13 = 1.d0 / 3.d0

! function
  phti =    (dplh(1,i,itm)                                                    &
  &    + at*(dplh(2,i,itm)*0.5d0 + at*(dplh(3,i,itm)*a13                      &
  &    + at*(dplh(4,i,itm)*0.25d0 + at*(dplh(5,i,itm)*0.2d0))))               &
  &    + dplh(6,i,itm)/at)                                                    &
  &    * dru * dmlr(i) * at

! original form
!  phti =    (dplh(1,i,itm)                                                    &
!  &    + at*(dplh(2,i,itm)/2.d0 + at*(dplh(3,i,itm)/3.d0                      &
!  &    + at*(dplh(4,i,itm)/4.d0 + at*(dplh(5,i,itm)/5.d0))))                  &
!  &    + dplh(6,i,itm)/at)                                                    &
!  &    * dru / dml(i) * at

!-----------------------------------------------------------------------

  end subroutine calc_phti


  end module
