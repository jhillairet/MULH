module FIELDScalc

contains
subroutine pFIELDSanalytic(geo,d4x,px,wave,sB,n,pB,pE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	ANALYTIC FIELDS FELT BY PARTICLE
!
! Calculate the fields felt by particle at their exact
! location by using the analytic formulas for the TE10
! mode inside a rectangular waveguide
!
! INPUTS
!   - geo: user defined type array with info regarding geometry of waveguide and PML
!   - d4x: real array with the discretization steps of the 4 dimensions
!   - px: 1D array with particle location
!   - wave: real array with information about RF EM wave
!   - n: integer time counter variable
!   - sB: real array with information about static B fields
!
! OUTPUTS
!   - pB: 1D array w/ B field (tesla) felt by particle
!   - pE: 1D array w/ E field (V/m) felt by particle
!

use precision_def
use def_types
use constants
implicit none
! Subroutine arguments
real(long), intent(inout) :: pB(3)
real(long), intent(inout), optional :: pE(3)
type(geo_t), intent(in) :: geo
real(long), intent(in) :: d4x(4), wave(4), px(3), sB(5)
integer, intent(in) :: n
! Unpacked variables
real(long) :: lw, b, dt, E0, w, beta, r0, rLH

!!!!!!!!! Unpack variables !!!!!!!!!!!!!
b = geo%b
lw = geo%lw
dt = d4x(4)
w = wave(1)
beta = wave(2)
E0 = wave(3)
r0 = sB(4)
rLH = sB(5)

! Calculate pEx
if (present(pE)) then
  pE(1) = E0*sin(pi*px(2)/b)*sin(w*n*dt-beta*px(3))
endif
                
! Advance Bx
pB(1) = sB(1)*r0/(lw-px(3)+rLH)

! Advance By
pB(2) = (beta/w)*E0*sin(pi*px(2)/b)*sin(w*(n+0.5)*dt-beta*px(3)) + sB(2)

! Advance Bz
pB(3) = (pi/(w*b))*E0*cos(pi*px(2)/b)*cos(w*(n+0.5)*dt-beta*px(3)) + sB(3)

end subroutine pFIELDSanalytic


subroutine fdtd_fields(geo,wave,x3,d4x,ramp,CE,DH,CEsigz,DHsigz, &
	n,Ex,Ey,Ez,Hx,Hy,Hz,ExPML,EyPML,EzPML,HxPML,HyPML,HzPML,Jx,Jy,Jz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!				FDTD FIELD SOLVER
!
! Calculate the E and H fields according to the FDTD method and implementing two Berenger split-PMLs
! at each end of the waveguide.
!
! INPUTS
!   - geo: user defined type array with info regarding geometry of waveguide and PML
!   - wave: real array with information about RF EM wave
!   - x3: real array with x, y and z arrays stored in it
!   - d4x: real array with the discretization steps of the 4 dimensions
!   - ramp: integer, number of periods over which fields are increased up to amplitude of E0
!   - C arrays: real 3D arrays with material/medium information of each node for computation of E fields
!   - D arrays: real 3D arrays with material/medium information of each node for computation of H fields
!   - EPML arrays: real 4D arrays with split E-field components
!   - HPML arrays: real 4D arrays with split H-field components
!   - E arrays: real 3D arrays with E field
!   - H arrays: real 3D arrays with H field
!   - J arrays: real arrays allocated for the computation of currents due to particles (optional)
!
! OUTPUTS
!   - EPML arrays: real 4D arrays with split E-field components
!   - HPML arrays: real 4D arrays with split H-field components
!   - E arrays: real 3D arrays with E field
!   - H arrays: real 3D arrays with H field
!
!

use precision_def
use constants
use def_types
! Subroutine arguments
real(long), intent(in) :: wave(4),d4x(4)
real(long), intent(in), allocatable :: CE(:,:,:,:,:),DH(:,:,:,:,:),CEsigz(:,:,:,:),DHsigz(:,:,:,:),x3(:)
real(long), intent(inout), allocatable :: Ex(:,:,:), Ey(:,:,:), Ez(:,:,:), Hx(:,:,:), Hy(:,:,:), Hz(:,:,:), &
					ExPML(:,:,:,:), EyPML(:,:,:,:), EzPML(:,:,:,:), HxPML(:,:,:,:), &
					HyPML(:,:,:,:), HzPML(:,:,:,:), Jx(:,:,:), Jy(:,:,:), Jz(:,:,:)
integer, intent(in) :: n, ramp
type (geo_t), intent(in) :: geo
! Unpacked variables
integer :: NOC(3), NOC_PML
real(long) :: lambda,dx,dy,dz,dt,E0,w,beta,b
real(long), allocatable :: y(:), z(:)
! Dummy variables
integer :: i,j,k
real(long) :: F, rampf

! Unpack variables
w = wave(1)
beta = wave(2)
E0 = wave(3)
lambda = wave(4)
dx = d4x(1)
dy = d4x(2)
dz = d4x(3)
dt = d4x(4)
b = geo%b
NOC = geo%NOC
NOC_PML = geo%NOC_PML
allocate (y(NOC(2)), z(NOC(3)))
y = x3(NOC(1)+1:NOC(1)+NOC(2))
z = x3(NOC(1)+NOC(2)+1:NOC(1)+NOC(2)+NOC(3))

! At first increase fields slowly in time with exponential factor over ramp periods
rampf = ramp * ceiling(lambda/(c*dt))
if (n < rampf) then
  F = 1./(1. + exp(2*0.05*(rampf*0.5 - n)))
else
  F = 1.
endif

!!!!!!!! Field solver !!!!!!!!
! The limits of every spatial loop are chosen so as to keep appropriate components unchanged (zero) to obey BCs

! Advance Hx
do i = 2,NOC(1)
  do j = 1,NOC(2)

    ! Inside the waveguide
    do k = 1,NOC(3)-NOC_PML-1
      Hx(i,j,k) = DH(i,j,k,1,1)*Hx(i,j,k) + DH(i,j,k,1,2)*( (Ey(i,j,k+1)-Ey(i,j,k))/dz - (Ez(i,j+1,k)-Ez(i,j,k))/dy)
    enddo

    ! Plasma-side PML-vacuum boundary Hx nodes
    k = NOC(3)-NOC_PML
    Hx(i,j,k) = DH(i,j,k,1,1)*Hx(i,j,k) + &
		DH(i,j,k,1,2)*( ((EyPML(i,j,k+1,1)+EyPML(i,j,k+1,2))-Ey(i,j,k))/dz - (Ez(i,j+1,k)-Ez(i,j,k))/dy)

    ! Plasma-side PML Hx nodes
    do k = NOC(3)-NOC_PML+1,NOC(3)
      HxPML(i,j,k,1) = DHsigz(i,j,k,1)*HxPML(i,j,k,1) - &
		DHsigz(i,j,k,2)*( (EzPML(i,j+1,k,1)+EzPML(i,j+1,k,2) - EzPML(i,j,k,1)-EzPML(i,j,k,2)) / dy)
      HxPML(i,j,k,2) = DH(i,j,k,1,1)*HxPML(i,j,k,2) + &
		DH(i,j,k,1,2)*( (EyPML(i,j,k+1,1)+EyPML(i,j,k+1,2) - EyPML(i,j,k,1)-EyPML(i,j,k,2)) / dz)
      Hx(i,j,k) = HxPML(i,j,k,1) + HxPML(i,j,k,2)
    enddo

  enddo
enddo
Hx = F * Hx
    
! Advance Hy
do i = 1,NOC(1)
  do j = 2,NOC(2)

    ! Inside the waveguide
    do k = 1,NOC(3)-NOC_PML-1
      Hy(i,j,k) = DH(i,j,k,2,1)*Hy(i,j,k) + DH(i,j,k,2,2)*( (Ez(i+1,j,k)-Ez(i,j,k))/dx - (Ex(i,j,k+1)-Ex(i,j,k))/dz)
    enddo

    ! Plasma-side PML-vacuum boundary Hy nodes
    k = NOC(3)-NOC_PML
    Hy(i,j,k) = DH(i,j,k,2,1)*Hy(i,j,k) + &
		DH(i,j,k,2,2)*( (Ez(i+1,j,k)-Ez(i,j,k))/dx - ((ExPML(i,j,k+1,1)+ExPML(i,j,k+1,2))-Ex(i,j,k)) / dz)

    ! Plasma-side PML Hy nodes
    do k = NOC(3)-NOC_PML+1,NOC(3)
      HyPML(i,j,k,1) = DH(i,j,k,2,1)*HyPML(i,j,k,1) - &
		DH(i,j,k,2,2)*( (ExPML(i,j,k+1,1)+ExPML(i,j,k+1,2) - ExPML(i,j,k,1)-ExPML(i,j,k,2)) / dz)
      HyPML(i,j,k,2) = DHsigz(i,j,k,1)*HyPML(i,j,k,2) + &
		DHsigz(i,j,k,2)*( (EzPML(i+1,j,k,1)+EzPML(i+1,j,k,2) - EzPML(i,j,k,1)-EzPML(i,j,k,2)) / dx)
      Hy(i,j,k) = HyPML(i,j,k,1) + HyPML(i,j,k,2)
    enddo

  enddo
enddo
Hy = F * Hy
    
! Advance Hz
do i = 1,NOC(1)
  do j = 1,NOC(2)

    ! Inside the waveguide
    do k = 1,NOC(3)-NOC_PML
      Hz(i,j,k) = DH(i,j,k,3,1)*Hz(i,j,k) + DH(i,j,k,3,2)*( (Ex(i,j+1,k)-Ex(i,j,k))/dy - (Ey(i+1,j,k)-Ey(i,j,k))/dx)
    enddo

    ! Plasma-side PML Hz nodes
    do k = NOC(3)-NOC_PML+1,NOC(3)
      HzPML(i,j,k,1) = DHsigz(i,j,k,1)*HzPML(i,j,k,1) - &
		DHsigz(i,j,k,2)*( (EyPML(i+1,j,k,1)+EyPML(i+1,j,k,2) - EyPML(i,j,k,1)-EyPML(i,j,k,2)) / dx)
      HzPML(i,j,k,2) = DHsigz(i,j,k,1)*HzPML(i,j,k,2) + &
		DHsigz(i,j,k,2)*( (ExPML(i,j+1,k,1)+ExPML(i,j+1,k,2) - ExPML(i,j,k,1)-ExPML(i,j,k,2)) / dy)
      Hz(i,j,k) = HzPML(i,j,k,1) + HzPML(i,j,k,2)
    enddo

  enddo
enddo
Hz = F * Hz

! Advance Ex
do i = 1,NOC(1)
  do j = 2,NOC(2)

    Ex(i,j,1) = E0*sin(pi*(y(j)-dy/2)/b)*sin(w*(n*dt+dt)-beta*(z(1)-dz/2))    ! Impose layer driving TE_10 mode
    ! Inside the waveguide
    do k = 2,NOC(3)-NOC_PML
      Ex(i,j,k) = CE(i,j,k,1,1)*Ex(i,j,k) + CE(i,j,k,1,2)*( (Hz(i,j,k)-Hz(i,j-1,k))/dy - (Hy(i,j,k)-Hy(i,j,k-1))/dz - Jx(i,j,k))
    enddo

    ! Plasma-side PML-vacuum boundary Ex nodes
    k = NOC(3)-NOC_PML+1
    ExPML(i,j,k,1) = CEsigz(i,j,k,1)*ExPML(i,j,k,1) + &
		CEsigz(i,j,k,2)*( (HzPML(i,j,k,1)+HzPML(i,j,k,2) - HzPML(i,j-1,k,1)-HzPML(i,j-1,k,2)) / dy)
    ExPML(i,j,k,2) = CE(i,j,k,1,1)*ExPML(i,j,k,2) - CE(i,j,k,1,2)*( (HyPML(i,j,k,1)+HyPML(i,j,k,2) - Hy(i,j,k-1)) / dz)
    Ex(i,j,k) = ExPML(i,j,k,1) + ExPML(i,j,k,2)

    ! Plasma-side PML Ex nodes
    do k = NOC(3)-NOC_PML+2,NOC(3)
      ExPML(i,j,k,1) = CEsigz(i,j,k,1)*ExPML(i,j,k,1) + &
		CEsigz(i,j,k,2)*( (HzPML(i,j,k,1)+HzPML(i,j,k,2) - HzPML(i,j-1,k,1)-HzPML(i,j-1,k,2)) / dy)
      ExPML(i,j,k,2) = CE(i,j,k,1,1)*ExPML(i,j,k,2) - &
		CE(i,j,k,1,2)*( (HyPML(i,j,k,1)+HyPML(i,j,k,2) - HyPML(i,j,k-1,1)-HyPML(i,j,k-1,2)) / dz)
      Ex(i,j,k) = ExPML(i,j,k,1) + ExPML(i,j,k,2)
    enddo

  enddo
enddo
Ex = F * Ex
    
! Advance Ey
do i = 2,NOC(1)
  do j = 1,NOC(2)

    ! Inside the waveguide
    do k = 2,NOC(3)-NOC_PML
      Ey(i,j,k) = CE(i,j,k,2,1)*Ey(i,j,k) + CE(i,j,k,2,2)*( (Hx(i,j,k)-Hx(i,j,k-1))/dz - (Hz(i,j,k)-Hz(i-1,j,k))/dx - Jy(i,j,k))
    enddo

    ! Plasma-side PML-vacuum boundary Ey nodes
    k = NOC(3)-NOC_PML+1
    EyPML(i,j,k,1) = CE(i,j,k,2,1)*EyPML(i,j,k,1) + &
		CE(i,j,k,2,2)*( (HxPML(i,j,k,1)+HxPML(i,j,k,2) - Hx(i,j,k-1)) / dz)
    EyPML(i,j,k,2) = CEsigz(i,j,k,1)*EyPML(i,j,k,2) - &
		CEsigz(i,j,k,2)*( (HzPML(i,j,k,1)+HzPML(i,j,k,2) - HzPML(i-1,j,k,1)-HzPML(i-1,j,k,2)) / dx)
    Ey(i,j,k) = EyPML(i,j,k,1) + EyPML(i,j,k,2)

    ! Plasma-side PML Ey nodes
    do k = NOC(3)-NOC_PML+2,NOC(3)
      EyPML(i,j,k,1) = CE(i,j,k,2,1)*EyPML(i,j,k,1) + &
		CE(i,j,k,2,2)*( (HxPML(i,j,k,1)+HxPML(i,j,k,2) - HxPML(i,j,k-1,1)-HxPML(i,j,k-1,2)) / dz)
      EyPML(i,j,k,2) = CEsigz(i,j,k,1)*EyPML(i,j,k,2) - &
		CEsigz(i,j,k,2)*( (HzPML(i,j,k,1)+HzPML(i,j,k,2) - HzPML(i-1,j,k,1)-HzPML(i-1,j,k,2)) / dx)
      Ey(i,j,k) = EyPML(i,j,k,1) + EyPML(i,j,k,2)
    enddo

  enddo
enddo
Ey = F * Ey
    
! Advance Ez
do i = 2,NOC(1)
  do j = 2,NOC(2)

    ! Inside the waveguide
    do k = 1,NOC(3)-NOC_PML
      Ez(i,j,k) = CE(i,j,k,3,1)*Ez(i,j,k) + CE(i,j,k,3,2)*( (Hy(i,j,k)-Hy(i-1,j,k))/dx - (Hx(i,j,k)-Hx(i,j-1,k))/dy - Jz(i,j,k))
    enddo

    ! Source-side PML Ez nodes
    do k = NOC(3)-NOC_PML+1,NOC(3)
      EzPML(i,j,k,1) = CEsigz(i,j,k,1)*EzPML(i,j,k,1) + &
		CEsigz(i,j,k,2)*( (HyPML(i,j,k,1)+HyPML(i,j,k,2) - HyPML(i-1,j,k,1)-HyPML(i-1,j,k,2)) / dx)
      EzPML(i,j,k,2) = CEsigz(i,j,k,1)*EzPML(i,j,k,2) - &
		CEsigz(i,j,k,2)*( (HxPML(i,j,k,1)+HxPML(i,j,k,2) - HxPML(i,j-1,k,1)-HxPML(i,j-1,k,2)) / dy)
      Ez(i,j,k) = EzPML(i,j,k,1) + EzPML(i,j,k,2)
    enddo
  enddo
enddo
Ez = F * Ez

end subroutine fdtd_fields

end module FIELDScalc
