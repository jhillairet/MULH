module p_locate

contains
subroutine plocate(px,d4x,x3,geo,pwall,Eii,Eip,Eri,Hii,Hip,Hri)
!**************************************************************
!				plocate
! Locates the particle in the E and H grids and produces ii,ip,wall arrays
! indicating where the particle is in those grids. Used to interpolate fields to
! particle

! INPUTS
!   - px: location of the particle
!   - d4x: array containing the step sizes dx, dy, dz and dt
!   - x3: array containing the spatial arrays x, y and z
!   - geo: user defined array containing geometry specs
!   - pwall: boolean to identify if the particle is between a wall and nodes
!   along each direction

! OUTPUTS
!   - ii: node to the left (towards the 1st node) of the particle along
!   each direction
!   - ip: node to the right (towards the last node) of the particle along
!   each direction
!   - ri: where the particle is on the grid
!   - pwall: boolean to identify if the particle is between a wall and nodes
!   along each direction
!
!

use precision_def
use def_types
implicit none
! Subroutine arguments
real(long), intent(in) :: px(3)
real(long), intent(in) :: d4x(4)
real(long), allocatable, intent(in) :: x3(:)
type(geo_t), intent(in) :: geo
integer, intent(out) :: pwall(4)
integer, intent(out) :: Eii(3,3),Eip(3,3),Hii(3,3),Hip(3,3)
real(long), intent(out) :: Eri(3,3), Hri(3,3)
! Unpacked variables
real(long) :: a,b
! Dummy variables
real(long) :: dx2,dy2
integer :: i,j,NOC(3)

a = geo%a
b = geo%b
NOC = geo%NOC

! Divide spatial steps by 2 for ease of coding and to decrease number of operations performed
dx2 = d4x(1)/2    
dy2 = d4x(2)/2

! Eri is where the particle is on the Egrid.  Eri is a real 
! number, and Eii and Eip are the integer Egrid positions to the left and 
! right of the particle location. (From Denton's espc.m)
Eri = 0.
Hri = 0.

do i = 1,3
  do j = 1,3

    if (j == 1) then
      Eri(i,j) = ( px(j) - x3(1) )/d4x(j) + 1.5
      Hri(i,j) = ( px(j) - x3(1) )/d4x(j) + 1.
    elseif (j == 2) then
      Eri(i,j) = ( px(j) - x3(NOC(1)+1) )/d4x(j) + 1.5
      Hri(i,j) = ( px(j) - x3(NOC(1)+1) )/d4x(j) + 1.
    else
      Eri(i,j) = ( px(j) - x3(NOC(1)+NOC(2)+1) )/d4x(j) + 1.5
      Hri(i,j) = ( px(j) - x3(NOC(1)+NOC(2)+1) )/d4x(j) + 1.
    endif

    if (i == j) then
      Eri(i,j) = Eri(i,j) - 0.5
      Hri(i,j) = Hri(i,j) + 0.5
    endif

  enddo
enddo

do i = 1,3
  do j = 1,3
    Eii(i,j) = floor(Eri(i,j))
    Eip(i,j) = Eii(i,j) + 1

    Hii(i,j) = floor(Hri(i,j))
    Hip(i,j) = Hii(i,j) + 1
  enddo
enddo

! The following are flags to identify if the particle is between a wall and a node.
pwall = 0
        
if (px(1) <= dx2) then
  pwall(3) = 1	! Particle is between wall and node along x
endif

if (px(1) >= (a-dx2)) then
  pwall(1) = 1	! Particle is between wall and node along x
endif
        
if (px(2) <= dy2) then
  pwall(4) = 1	! Particle is between wall and node along y
endif

if (px(2) >= (b-dy2)) then
  pwall(2) = 1
endif

end subroutine plocate

end module p_locate
