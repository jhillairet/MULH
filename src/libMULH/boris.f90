module boris

contains
subroutine boris_solver(dt,Ef,Bf,v2,v1,x2,x1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!		Boris solver
!
! Calculates the position of the particle at time n+1 and the
! velocity of the particle at n+1/2 using the Boris scheme.
!
! INPUTS
!   - dt: real time step
!   - Ef: real 1D array w/ E field felt by particle
!   - Bf: real 1D array w/ B field felt by particle
!   - v1: real 1D array w/ particle velocity at time n-1/2
!   - x1: real 1D array w/ particle position at time n
!
! OUTPUTS
!   - v2: real 1D array w/ particle velocity at time n+1/2
!   - x2: real 1D array w/ particle position at time n+1
!

use precision_def
use def_types
use constants
implicit none
! Subroutine arguments
real(long), intent(inout) :: v2(3)
real(long), intent(inout), optional :: x2(3)
real(long), intent(in) :: dt, Ef(3), Bf(3), v1(3)
real(long), intent(in), optional :: x1(3)
! Dummy variables
real(long) :: stepf, vm(3), t(3), vp(3), vs(3), s(3)

! Useful factor
stepf = eme * dt/2

! Step 1: obtain v-
vm = v1 + stepf * Ef
    
! Step 2: obtain v'. Time average B field to obtain B at integer time
t = Bf*stepf
    
vp(1) = vm(1) + vm(2)*t(3) - vm(3)*t(2)
vp(2) = vm(2) + vm(3)*t(1) - vm(1)*t(3)
vp(3) = vm(3) + vm(1)*t(2) - vm(2)*t(1)

! Step 3: obtain v*
s = (t*2)/(1+t(1)*t(1)+t(2)*t(2)+t(3)*t(3))

vs(1) = vm(1) + vp(2)*s(3) - vp(3)*s(2)
vs(2) = vm(2) + vp(3)*s(1) - vp(1)*s(3)
vs(3) = vm(3) + vp(1)*s(2) - vp(2)*s(1)
    
! Step 4: update pv by half stepping due to E field
v2 = vs + stepf * Ef

! Step particles / update position if necessary
if (present(x1)) then
  x2 = x1 + dt*v2
endif

end subroutine boris_solver

end module boris
