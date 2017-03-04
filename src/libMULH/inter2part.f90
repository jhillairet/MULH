module inter2part

contains
function interB2particle(Hx,Hy,Hz,Hii,Hip,Hri,pwall) result(pB)
!*************************************************************************
!			interB2particle
!   Interpolate all three components of the magnetic field to
!   the desired location of a particle. Outputs the fields felt by the
!   particles.
!   - It outputs the magnetic field B, not H.
!
! INPUTS
!   - Hx, Hy, Hz: components of H-fields required to interpolate to
!   particle
!   - Hii: H node to the left (towards the 1st node) of the particle along
!   each direction
!   - Hip: H node to the right (towards the last node) of the particle along
!   each direction
!   - Hri: where the particle is on the H grid
!   - pwall: boolean to identify if the particle is between a wall and nodes
!   along each direction
!
! OUTPUTS
!   - pB: magnetic field felt by the particle
!

use precision_def
use constants
implicit none
! Subroutine arguments
real(long), allocatable, intent(in) :: Hx(:,:,:), Hy(:,:,:), Hz(:,:,:)
integer, intent(inout) :: Hii(3,3),Hip(3,3)
real(long), intent(in) :: Hri(3,3)
integer, intent(inout) :: pwall(4)
real(long) :: pB(3)
! Dummy variables
integer :: i, j, k
real(long) :: fii(3), fip(3)

! Interpolate Hx to the particle
pB = 0.
fii = 0.
fip = 0.

do i = 1,3
  do j = 1,3
    fip(j) = Hri(i,j) - Hii(i,j)
    fii(j) = 1. - fip(j)
  enddo
    
  if (i == 2) then

    if (pwall(3) == 1) then
      Hii(2,1) = Hip(2,1)
      fip(1) = fip(1) + fii(1)
      fii(1) = 0.
    endif
        
    if (pwall(1) == 1) then
      Hip(2,1) = Hii(2,1)
      fii(1) = fip(1) + fii(1)
      fip(1) = 0.
    endif

    pB(i) = Hy(Hii(i,1),Hii(i,2),Hii(i,3))*fii(1)*fii(2)*fii(3) + Hy(Hip(i,1),Hii(i,2),Hii(i,3))*fip(1)*fii(2)*fii(3) + &
		Hy(Hii(i,1),Hip(i,2),Hii(i,3))*fii(1)*fip(2)*fii(3) + Hy(Hii(i,1),Hii(i,2),Hip(i,3))*fii(1)*fii(2)*fip(3) + &
		Hy(Hip(i,1),Hip(i,2),Hii(i,3))*fip(1)*fip(2)*fii(3) + Hy(Hip(i,1),Hii(i,2),Hip(i,3))*fip(1)*fii(2)*fip(3) + &
		Hy(Hii(i,1),Hip(i,2),Hip(i,3))*fii(1)*fip(2)*fip(3) + Hy(Hip(i,1),Hip(i,2),Hip(i,3))*fip(1)*fip(2)*fip(3)

  elseif (i == 1) then
    
    if (pwall(2) == 1) then
      Hip(1,2) = Hii(1,2)
      fii(2) = fip(2) + fii(2)
      fip(2) = 0.
    endif
        
    if (pwall(4) == 1) then
      Hii(1,2) = Hip(1,2)
      fip(2) = fip(2) + fii(2)
      fii(2) = 0.
    endif

    pB(i) = Hx(Hii(i,1),Hii(i,2),Hii(i,3))*fii(1)*fii(2)*fii(3) + Hx(Hip(i,1),Hii(i,2),Hii(i,3))*fip(1)*fii(2)*fii(3) + &
		Hx(Hii(i,1),Hip(i,2),Hii(i,3))*fii(1)*fip(2)*fii(3) + Hx(Hii(i,1),Hii(i,2),Hip(i,3))*fii(1)*fii(2)*fip(3) + &
		Hx(Hip(i,1),Hip(i,2),Hii(i,3))*fip(1)*fip(2)*fii(3) + Hx(Hip(i,1),Hii(i,2),Hip(i,3))*fip(1)*fii(2)*fip(3) + &
		Hx(Hii(i,1),Hip(i,2),Hip(i,3))*fii(1)*fip(2)*fip(3) + Hx(Hip(i,1),Hip(i,2),Hip(i,3))*fip(1)*fip(2)*fip(3)

  elseif (i == 3) then

    if (pwall(3) == 1) then
      Hii(3,1) = Hip(3,1)
      fip(1) = fip(1) + fii(1)
      fii(1) = 0.
    endif
    
    if (pwall(1) == 1) then
      Hip(3,1) = Hii(3,1)
      fii(1) = fip(1) + fii(1)
      fip(1) = 0.
    endif
    
    if (pwall(2) == 1) then
      Hip(3,2) = Hii(3,2)
      fii(2) = fip(2) + fii(2)
      fip(2) = 0.
    endif
        
    if (pwall(4) == 1) then
      Hii(3,2) = Hip(3,2)
      fip(2) = fip(2) + fii(2)
      fii(2) = 0.
    endif

    pB(i) = Hz(Hii(i,1),Hii(i,2),Hii(i,3))*fii(1)*fii(2)*fii(3) + Hz(Hip(i,1),Hii(i,2),Hii(i,3))*fip(1)*fii(2)*fii(3) + &
		Hz(Hii(i,1),Hip(i,2),Hii(i,3))*fii(1)*fip(2)*fii(3) + Hz(Hii(i,1),Hii(i,2),Hip(i,3))*fii(1)*fii(2)*fip(3) + &
		Hz(Hip(i,1),Hip(i,2),Hii(i,3))*fip(1)*fip(2)*fii(3) + Hz(Hip(i,1),Hii(i,2),Hip(i,3))*fip(1)*fii(2)*fip(3) + &
		Hz(Hii(i,1),Hip(i,2),Hip(i,3))*fii(1)*fip(2)*fip(3) + Hz(Hip(i,1),Hip(i,2),Hip(i,3))*fip(1)*fip(2)*fip(3)

  endif

enddo
    
pB = pB * mu_0

end function interB2particle

function interE2particle(Ex,Ey,Ez,Eii,Eip,Eri,pwall) result(pE)
!*************************************************************************
!			interE2particle
!   Interpolate all three components of the electric field to
!   the desired location of a particle. Outputs the fields felt by the
!   particles.
!
! INPUTS
!   - Ex, Ey, Ez: components of electric field required to interpolate to
!   particle
!   - Eii: E node to the left (towards the 1st node) of the particle along
!   each direction
!   - Eip: E node to the right (towards the last node) of the particle along
!   each direction
!   - Eri: where the particle is on the E grid
!   - pwall: boolean to identify if the particle is between a wall and nodes
!   along each direction
!
! OUTPUTS
!   - pE: magnetic field felt by the particle
!

use precision_def
use def_types
use constants
implicit none
! Subroutine arguments
real(long), allocatable, intent(in) :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:)
integer, intent(inout) :: Eii(3,3),Eip(3,3)
real(long), intent(in) :: Eri(3,3)
integer, intent(inout) :: pwall(4)
real(long) :: pE(3)
! Dummy variables
integer :: i, j, k
real(long) :: fii(3), fip(3)

! Interpolate E to the particle
pE = 0.
fii = 0.
fip = 0.

do i = 1,3
  do j = 1,3        
    fip(j) = Eri(i,j) - Eii(i,j)
    fii(j) = 1 - fip(j)
  enddo
    
  if (pwall(3)==1 .AND. i==1) then
    Eii(1,1) = Eip(1,1)
    fip(1) = fip(1) + fii(1)
    fii(1) = 0.
  endif
    
  if (pwall(1) == 1 .AND. i==1) then
    Eip(1,1) = Eii(1,1)
    fii(1) = fip(1) + fii(1)
    fip(1) = 0.
  endif
    
  if (pwall(2) == 1 .AND. i==2) then
    Eip(2,2) = Eii(2,2)
    fii(2) = fip(2) + fii(2)
    fip(2) = 0.
  endif
    
  if (pwall(4) == 1 .AND. i==2) then
    Eii(2,2) = Eip(2,2)
    fip(2) = fip(2) + fii(2)
    fii(2) = 0.
  endif 

  if (i==1) then
    pE(i) = Ex(Eii(i,1),Eii(i,2),Eii(i,3))*fii(1)*fii(2)*fii(3) + Ex(Eip(i,1),Eii(i,2),Eii(i,3))*fip(1)*fii(2)*fii(3) + &
		Ex(Eii(i,1),Eip(i,2),Eii(i,3))*fii(1)*fip(2)*fii(3) + Ex(Eii(i,1),Eii(i,2),Eip(i,3))*fii(1)*fii(2)*fip(3) + &
		Ex(Eip(i,1),Eip(i,2),Eii(i,3))*fip(1)*fip(2)*fii(3) + Ex(Eip(i,1),Eii(i,2),Eip(i,3))*fip(1)*fii(2)*fip(3) + &
		Ex(Eii(i,1),Eip(i,2),Eip(i,3))*fii(1)*fip(2)*fip(3) + Ex(Eip(i,1),Eip(i,2),Eip(i,3))*fip(1)*fip(2)*fip(3)
  elseif (i==2) then
    pE(i) = Ey(Eii(i,1),Eii(i,2),Eii(i,3))*fii(1)*fii(2)*fii(3) + Ey(Eip(i,1),Eii(i,2),Eii(i,3))*fip(1)*fii(2)*fii(3) + &
		Ey(Eii(i,1),Eip(i,2),Eii(i,3))*fii(1)*fip(2)*fii(3) + Ey(Eii(i,1),Eii(i,2),Eip(i,3))*fii(1)*fii(2)*fip(3) + &
		Ey(Eip(i,1),Eip(i,2),Eii(i,3))*fip(1)*fip(2)*fii(3) + Ey(Eip(i,1),Eii(i,2),Eip(i,3))*fip(1)*fii(2)*fip(3) + &
		Ey(Eii(i,1),Eip(i,2),Eip(i,3))*fii(1)*fip(2)*fip(3) + Ey(Eip(i,1),Eip(i,2),Eip(i,3))*fip(1)*fip(2)*fip(3)
  elseif (i==3) then
    pE(i) = Ez(Eii(i,1),Eii(i,2),Eii(i,3))*fii(1)*fii(2)*fii(3) + Ez(Eip(i,1),Eii(i,2),Eii(i,3))*fip(1)*fii(2)*fii(3) + &
		Ez(Eii(i,1),Eip(i,2),Eii(i,3))*fii(1)*fip(2)*fii(3) + Ez(Eii(i,1),Eii(i,2),Eip(i,3))*fii(1)*fii(2)*fip(3) + &
		Ez(Eip(i,1),Eip(i,2),Eii(i,3))*fip(1)*fip(2)*fii(3) + Ez(Eip(i,1),Eii(i,2),Eip(i,3))*fip(1)*fii(2)*fip(3) + &
		Ez(Eii(i,1),Eip(i,2),Eip(i,3))*fii(1)*fip(2)*fip(3) + Ez(Eip(i,1),Eip(i,2),Eip(i,3))*fip(1)*fip(2)*fip(3)
  endif

enddo


end function interE2particle

end module inter2part
