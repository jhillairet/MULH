module prep_fields

contains
subroutine prep_ex_fields(geo,wave,ndt,d4x,x3,P_I,sB,Exn,Eyn,Ezn,Hxn,Hyn,Hzn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!				PREPARE EXTERNAL FIELDS
!
! This subroutine prepares the En and Hn arrays with the fields at all times coming from an
! external solver (like COMSOL).
!
! At the moment it adapts the fields obtained from Melanie Preyna's Piccolo COMSOL Model. These
! fields have been processed with the MATLAB script COMSOLfields.m. If other fields from another
! solver need to be used, some of this code needs to be changed. However, the basic output is to
! produce the En and Hn arrays which contain the fields at all nodes in the geometry at all times
! during one period of the EM wave.
!
! INPUTS
!   - geo: user defined type array with info regarding geometry of waveguide and PML
!   - wave: real array with information about RF EM wave
!   - ndt: real; the period of the RF EM fields stored
!   - d4x: real array with the discretization steps of the 4 dimensions
!   - x3: real array with x, y and z arrays stored in it
!   - P_I: integer power of RF wave
!   - sB: real array with information about static B fields
!
! OUTPUTS
!   - En arrays: E fields at all times
!   - Hn arrays: H fields at all times
!

use precision_def
use constants
use def_types
! Subroutine arguments
integer, intent(in):: ndt, P_I
real(long), intent(in) :: sB(5),wave(4),d4x(4)
real(long), intent(in), allocatable :: x3(:)
real(long), intent(out), allocatable :: Exn(:,:,:,:),Eyn(:,:,:,:),Ezn(:,:,:,:),Hxn(:,:,:,:),Hyn(:,:,:,:),Hzn(:,:,:,:)
type (geo_t), intent(in) :: geo
! Unpacked variables
integer :: NOC(3)
real(long) :: a,b,lw,dy,dz,dt,beta,E0,w,r0,rLH
real(long), allocatable :: y(:), z(:)
! Dummy variables
integer :: i,j,k,n
real(long) :: t
real(long), allocatable :: enorm(:,:),hnorm(:,:),phie(:,:),phih(:,:)
 character(4) :: P,NOCc

!!!!!!!! Unpack variables !!!!!!!!!
w = wave(1)
beta = wave(2)
E0 = wave(3)
dy = d4x(2)
dz = d4x(3)
dt = d4x(4)
NOC = geo%NOC
a = geo%a
b = geo%b
lw = geo%lw
r0 = sB(4)
rLH = sB(5)
allocate (y(NOC(2)), z(NOC(3)))
y = x3(NOC(1)+1:NOC(1)+NOC(2))
z = x3(NOC(1)+NOC(2)+1:NOC(1)+NOC(2)+NOC(3))

! Allocate arrays where E and H fields will be stored for all times
allocate (Exn(NOC(1),NOC(2)+1,NOC(3)+1,ndt))
allocate (Eyn(NOC(1)+1,NOC(2),NOC(3)+1,ndt))
allocate (Ezn(NOC(1)+1,NOC(2)+1,NOC(3),ndt))
Exn = 0.
Eyn = 0.
Ezn = 0.
allocate (Hxn(NOC(1)+1,NOC(2),NOC(3),ndt))
allocate (Hyn(NOC(1),NOC(2)+1,NOC(3),ndt))
allocate (Hzn(NOC(1),NOC(2),NOC(3)+1,ndt))
Hxn = 0.
Hyn = 0.
Hzn = 0.

!***** USE FIELDS CALCULATED w/ PICCOLO (COMSOL model) *****!
!
! For simplicity and speed, fields are first processed with MATLAB using the script
! comsolfieds.m. Such program produces enorm, hnorm, phie and phih, which are used below
! to generate Exn and Hyn.
! The parameters in comsolfields.m (NOC, NOC_lambda, etc) have to agree w/ those used in MULH.
!

allocate ( enorm(NOC(3)+1,NOC(1)), phie(NOC(3)+1,NOC(1)), hnorm(NOC(3),NOC(1)), phih(NOC(3),NOC(1)) )

!write(P,'(I4)') P_I
!write(NOCc, '(I4)') NOC_lambda

! Load norm and angle of E and H fields
!open(unit=25,file='../data/P'//trim(P)//'/NOC'//trim(NOCc)//'/enorm.txt',status='old')
!open(unit=26,file='../data/P'//trim(P)//'/NOC'//trim(NOCc)//'/phie.txt',status='old')
!open(unit=27,file='../data/P'//trim(P)//'/NOC'//trim(NOCc)//'/hnorm.txt',status='old')
!open(unit=28,file='../data/P'//trim(P)//'/NOC'//trim(NOCc)//'/phih.txt',status='old')
open(unit=25,file='../data/enorm.txt',status='old')
open(unit=26,file='../data/phie.txt',status='old')
open(unit=27,file='../data/hnorm.txt',status='old')
open(unit=28,file='../data/phih.txt',status='old')

do k = 1,NOC(3)+1

  read(25,*) enorm(k,:)
  read(26,*) phie(k,:)

  if (k < NOC(3)+1) then
    read(27,*) hnorm(k,:)
    read(28,*) phih(k,:)
  endif

enddo

 close(unit=25)
 close(unit=26)
 close(unit=27)
 close(unit=28)

do n = 1,ndt

  t = n*dt
  do i = 1,NOC(1)
    do j = 1,NOC(2)+1
      do k = 1,NOC(3)+1
	Exn(i,j,k,n) = E0 * enorm(k,i) * sin(pi*(j-1)*dy/b) * cos(w*t+phie(k,i))
      enddo
    enddo
  enddo

  do i = 1,NOC(1)
    do j = 1,NOC(2)+1
      do k = 1,NOC(3)
	Hyn(i,j,k,n) = (beta/(w*mu_0)) * E0 * hnorm(k,i) * sin(pi*(j-1)*dy/b) * cos(w*t+phih(k,i)) + sB(2)/mu_0
      enddo
    enddo
  enddo

  ! Hx according to the static toroidal magnetic field
  do i = 1,NOC(1)+1
    do j = 1,NOC(2)
      do k = 1,NOC(3)
	Hxn(i,j,k,n) = sB(1)*r0/(mu_0*(lw-z(k)+rLH))
      enddo
    enddo
  enddo

  ! Artificial Hzn to account for the rectangular geometry instead of PICCOLO's parallel plate geometry
  do i = 1,NOC(1)
    do j = 1,NOC(2)
      do k = 1,NOC(3)+1
        if (k > 1 .AND. k < NOC(3)+1) then
	  Hzn(i,j,k,n) = (pi/(w*b*mu_0)) * E0 * cos(pi*y(j)/b) * cos(w*(n+0.5)*dt-beta*(z(k)-dz*0.5)) + sB(3)/mu_0
        else
          Hzn(i,j,k,n) = sB(3)/mu_0
        endif
      enddo
    enddo
  enddo

enddo

end subroutine prep_ex_fields



subroutine prep_fdtd_fields(geo,wave,d4x,x3,CE,DH,CEsigz,DHsigz, &
	ExPML,EyPML,EzPML,HxPML,HyPML,HzPML,Ex2,Hy2,Hz2,Jx,Jy,Jz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!				PREPARE FDTD FIELD SOLVER ARRAYS
!
! This subroutine prepares the C and D arrays used in the FDTD field solver and sets up the 
! initial conditions on the E and H fields.
!
! This is necessary if space charge effects are taken into account.
!
! INPUTS
!   - geo: user defined type array with info regarding geometry of waveguide and PML
!   - wave: real array with information about RF EM wave
!   - d4x: real array with the discretization steps of the 4 dimensions
!   - x3: real array with x, y and z arrays stored in it
!
! OUTPUTS
!   - C arrays: real 3D arrays with material/medium information of each node for computation of E fields
!   - D arrays: real 3D arrays with material/medium information of each node for computation of H fields
!   - EPML arrays: real 4D arrays with split E-field components
!   - HPML arrays: real 4D arrays with split H-field components
!   - Ex2, Hy2, Hz2: intial conditions
!   - J arrays: real arrays allocated for the computation of currents due to particles (optional)
!
!
use precision_def
use constants
use def_types
! Subroutine arguments
real(long), intent(in) :: wave(4),d4x(4)
real(long), intent(in), allocatable :: x3(:)
real(long), intent(out), allocatable :: CE(:,:,:,:,:), DH(:,:,:,:,:),CEsigz(:,:,:,:),DHsigz(:,:,:,:), &
		ExPML(:,:,:,:), EyPML(:,:,:,:), EzPML(:,:,:,:), HxPML(:,:,:,:), HyPML(:,:,:,:), HzPML(:,:,:,:)
real(long), intent(inout), allocatable :: Ex2(:,:,:), Hy2(:,:,:), Hz2(:,:,:), Jx(:,:,:), Jy(:,:,:), Jz(:,:,:)
type (geo_t), intent(in) :: geo
! Unpacked variables
integer :: NOC(3), NOC_PML, m_PML
real(long) :: a,b,lw,dy,dz,dt,beta,E0,w,R_max
real(long), allocatable :: y(:), z(:)
! Dummy variables
integer :: i,j,k
real(long) :: t, sigma_zmax, PMLfactor, z0, sigma_z, F, eta, n

!!!!!!!! Unpack variables !!!!!!!!!
w = wave(1)
beta = wave(2)
E0 = wave(3)
dy = d4x(2)
dz = d4x(3)
dt = d4x(4)
NOC = geo%NOC
NOC_PML = geo%NOC_PML
a = geo%a
b = geo%b
lw = geo%lw
allocate (y(NOC(2)), z(NOC(3)))
y = x3(NOC(1)+1:NOC(1)+NOC(2))
z = x3(NOC(1)+NOC(2)+1:NOC(1)+NOC(2)+NOC(3))

eta = sqrt(mu_0/eps_0)           ! Vacuum impedance
!!!!!!!! Allocate J arrays for computing current due to electrons !!!!!!!!
allocate (Jx(NOC(1),NOC(2)+1,NOC(3)+1),Jy(NOC(1)+1,NOC(2),NOC(3)+1),Jz(NOC(1)+1,NOC(2)+1,NOC(3)))
Jx = 0.
Jy = 0.
Jz = 0.

!!!!!!!! Allocate additional E and H arrays required for subcomponents in PML !!!!!!!
! In Berenger's split PML one needs to compute 
! Exy and Exz to obtain Ex
! Eyz and Eyx to obtain Ey
! Ezx and Ezy to obtain Ez
! Hxy and Hxz to obtain Hx
! Hyz and Hyx to obtain Hy
! Hzx and Hzy to obtain Hz
! The first of each pair will be stored in the corresponding E/H array. The second has an array of its own.
allocate (ExPML(NOC(1),NOC(2)+1,NOC(3)+1,2),EyPML(NOC(1)+1,NOC(2),NOC(3)+1,2),EzPML(NOC(1)+1,NOC(2)+1,NOC(3),2))
allocate (HxPML(NOC(1)+1,NOC(2),NOC(3),2),HyPML(NOC(1),NOC(2)+1,NOC(3),2),HzPML(NOC(1),NOC(2),NOC(3)+1,2))
ExPML = 0.
EyPML = 0.
EzPML = 0.
HxPML = 0.
HyPML = 0.
HzPML = 0.

!!!!!!!! Allocate C and D arrays. Material properties at field component locations (include dt) !!!!!!!
! First 3 dimensions of CE and DH are for the spatial components
! Fourth dimension specifies whether x, y or z component of E/H
! 5th dimension specifies whether to use the Ca/Da or the Cb/Db factor
allocate (CE(NOC(1)+1,NOC(2)+1,NOC(3)+1,3,2),DH(NOC(1)+1,NOC(2)+1,NOC(3)+1,3,2))
allocate (CEsigz(NOC(1)+1,NOC(2)+1,NOC(3)+1,2),DHsigz(NOC(1)+1,NOC(2)+1,NOC(3)+1,2))

 CE(:,:,:,:,1) = 1.
DH(:,:,:,:,1) = 1.
 CEsigz(:,:,:,1) = 1.
DHsigz(:,:,:,1) = 1.
 CE(:,:,:,:,2) = dt/eps_0
DH(:,:,:,:,2) = dt/mu_0
 CEsigz(:,:,:,2) = dt/eps_0	! Used because sigma_x and sigma_y are zero for a PML at the end of waveguide along z
DHsigz(:,:,:,2) = dt/mu_0	! Used because sigma_x and sigma_y are zero for a PML at the end of waveguide along z

R_max = geo%R_max
m_PML = geo%m_PML
!!!!!!!! Calculate parameters for polynomial-graded Berenger's PML and populate material property arrays C and D !!!!!!!
PMLfactor = - log(R_max)/( (2.**(m_PML+2)) * eta * dz * (NOC_PML**(m_PML+1)) )  

do i = 1,NOC(1)+1
  do j = 1,NOC(2)+1
        
    n = 0.
    do k = NOC(3)-NOC_PML+1,NOC(3)+1	! Plasma-side PML

      if (k == NOC(3)-NOC_PML+1) then
	sigma_z = PMLfactor		! Fornodes on Vacuum-PML
      else
	n = n + 0.5
	sigma_z = PMLfactor * ( (2*n+1)**(m_PML+1) - (2*n-1)**(m_PML+1) )
      endif
            
      ! Populate C and D according to Taflove1995 exponential decay FDTD-PML
      if (i < NOC(1)+1) then
        CE(i,j,k,1,1) = exp(-sigma_z*dt/eps_0)
        CE(i,j,k,1,2) = (1-exp(-sigma_z*dt/eps_0))/sigma_z
      endif
      if (j < NOC(2)+1) then
        CE(i,j,k,2,1) = exp(-sigma_z*dt/eps_0)
        CE(i,j,k,2,2) = (1-exp(-sigma_z*dt/eps_0))/sigma_z
      endif
      if ((i < NOC(1)+1) .AND. (j < NOC(2)+1)) then
        DH(i,j,k,3,1) = exp(-sigma_z*dt/eps_0)
        DH(i,j,k,3,2) = (1-exp(-sigma_z*dt/eps_0))*eps_0/(mu_0*sigma_z)
      endif
      

      if (k < NOC(3)+1) then
	n = n + 0.5
        sigma_z = PMLfactor * ( (2*n+1)**(m_PML+1) - (2*n-1)**(m_PML+1) )
    
        CE(i,j,k,3,1) = exp(-sigma_z*dt/eps_0)
        CE(i,j,k,3,2) = (1-exp(-sigma_z*dt/eps_0))/sigma_z
        if (j < NOC(2)+1) then
          DH(i,j,k,1,1) = exp(-sigma_z*dt/eps_0)
          DH(i,j,k,1,2) = (1-exp(-sigma_z*dt/eps_0))*eps_0/(mu_0*sigma_z)
        endif
        if (i < NOC(1)+1) then
          DH(i,j,k,2,1) = exp(-sigma_z*dt/eps_0)
          DH(i,j,k,2,2) = (1-exp(-sigma_z*dt/eps_0))*eps_0/(mu_0*sigma_z)
        endif
      endif
    
    enddo

  enddo
enddo

!!!! Initial conditions !!!!
! Use solutions for TE_10 in rectangular waveguide from Cheng's Field and
! Wave electromagnetics + "Multipactor in rectangular waveguides" Semenov
! Phys. Plasma 14 (2007).

! Impose initial conditions over all relevant components in Yee cell layer
! at source-side end of the waveguide. Bear in mind that spatial coords are not defined at
! same places as fields, that's why factors of dx and dz are subtracted.
! Notice that E and H fields are initialized w/ a dt/2 time difference.
do i = 1,NOC(1)
    do j = 1,NOC(2)
        if (j /= 1) then	! i = 1,NOC_Y+1 remain zeros because of PEC BCs
             Ex2(i,j,1) = E0 * sin(pi*(y(j)-dy/2)/b) * sin(w*dt-beta*(z(1)-dz/2))
        endif
        
    enddo
enddo

F = 1/( 1 + exp(2*0.05*5) )        ! Smoothing function to increment fields slowly
!F = 1

Ex2 = F * Ex2
!Hy2 = F * Hy2
!Hz2 = F * Hz2

end subroutine prep_fdtd_fields

end module prep_fields
