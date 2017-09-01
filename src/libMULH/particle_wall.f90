module particle_wall

contains
subroutine particlewall(px2t,pv2t,s,d4x,x3,geo,particle,sey,fields,wave,n,sB,SEYvsE,atype, &
		Hx0,Hy0,Hz0,Hx1,Hy1,Hz1,Ex1,Ey1,Ez1,Hx2,Hy2,Hz2,Ex2,Ey2,Ez2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!					PARTICLE-WALL INTERACTION
!
! Numerics concerning the interaction between particles and walls upon collision. The ejection and collection algorithms
! are based on Cheng G., Lie L. Comp. Phys. Com. 182 (2011) 1295-1303.
!
! INPUTS:
!   - d4x: real array with the discretization steps of the 4 dimensions 
!   - x3: real array with x, y and z arrays stored in it
!   - geo: user defined type array with info regarding geometry of waveguide and PML
!   - particle: user defined array containing all info about colliding particle
!   - sey: user defined array containig info about SEE model
!   - fields: integer indicating how EM fields are calculated
!   - wave: real array with information about RF EM wave
!   - n: integer, time iteration counter
!   - sB: array containing DC B fields and location of LH antenna
!   - SEYvsE: array used for storing values of SEY at each impacting Energy
!   - atype: integer, type of analysis being performed
!   - H arrays: real 3D arrays containg H fields
!   - E arrays: real 3D arrays containg E fields
!
! OUTPUTS:
!   - px2t: real 2D array with location of electrons coming off the wall
!   - pv2t: real 2D array with velocity of electrons coming off the wall
!   - s: number of electrons coming off the wall
!
!

use precision_def
use def_types
use constants
use p_locate
use inter2part
use FIELDScalc
use boris
use seem
implicit none
! Subroutine arguments
real(long), allocatable, intent(out) :: px2t(:,:), pv2t(:,:)
integer, intent(out) :: s, fields, n
integer, intent(in) :: atype
real(long), intent(in) :: d4x(4), wave(4), sB(5)
real(long), allocatable, intent(in) :: x3(:)
type(geo_t), intent(in) :: geo
type(particle_t), intent(in) :: particle
type(sey_t), intent(in) :: sey
real(long),allocatable, intent(inout) :: SEYvsE(:,:,:)
real(long), allocatable, intent(in), optional :: Ex2(:,:,:),Ey2(:,:,:),Ez2(:,:,:),Hx2(:,:,:),Hy2(:,:,:),Hz2(:,:,:), &
		Ex1(:,:,:),Ey1(:,:,:),Ez1(:,:,:),Hx1(:,:,:),Hy1(:,:,:),Hz1(:,:,:),Hx0(:,:,:),Hy0(:,:,:),Hz0(:,:,:)
! Unpacked variables
real(long) :: dt, a, b, lw, px1(3), pv1(3), pv2(3), pB0(3), pB1(3), pB2(3), pE1(3), pE2(3)
integer :: NOC(3), pcol(4)
real(long), allocatable :: z(:)
! Dummy variables
real(long) :: pvn(3), pBn(3), xb, f1dt, f1dtx, f1dty, f1, pvf1d2(3), pxf1(3), wB1(3), wE1(3), wB2(3), wE2(3), &
	gradEn(3,3), gradBn(3,3), dp, pEf1(3), pBf1(3), pvf1(3), pvd2f(3), wB0(3), eB1(3), eB2(3), eE2(3), dtp, &
	pEf(3), pBf(3), an(3), Eri(3,3), Hri(3,3)
real(long), allocatable :: pvf(:,:), px2temp(:,:), pv2temp(:,:)
integer :: i, j, l, o, k, pwall(4), Eii(3,3),Eip(3,3),Hii(3,3),Hip(3,3)
integer, allocatable :: s3(:)


!!!!!! Unpack variables !!!!!!
a = geo%a
b = geo%b
lw = geo%lw
NOC = geo%NOC
allocate (z(NOC(3)))
z = x3(NOC(1)+NOC(2)+1:NOC(1)+NOC(2)+NOC(3))
dt = d4x(4)
pcol = particle%pcol
px1 = particle%px1
pv1 = particle%pv1
pv2 = particle%pv2
pB0 = particle%pB0
pB1 = particle%pB1
pB2 = particle%pB2
pE1 = particle%pE1
pE2 = particle%pE2

! Initialize variables
wB1 = 0.
wB2 = 0.
wE1 = 0.
wE2 = 0.

! Velocity and time at time n
pvn = (pv1+pv2)*0.5
pBn = (pB1+pB2)*0.5

! Particle acceleration at time n
an(1) = eme*(pE2(1) + pvn(2)*pBn(3) - pvn(3)*pBn(2))
an(2) = eme*(pE2(2) + pvn(3)*pBn(1) - pvn(1)*pBn(3))
an(3) = eme*(pE2(3) + pvn(1)*pBn(2) - pvn(2)*pBn(1))

!!!!!! Detect whether particle crossed one or two walls !!!!!!
if (count(pcol .EQ. 1) == 1) then
  
  ! Detect which wall the particle collided with
  if (pcol(1) == 1) then
    xb = a
    i = 1
  elseif (pcol(2) == 1) then
    xb = b
    i = 2
  elseif (pcol(3) == 1) then
    xb = 0.
    i = 1
  elseif (pcol(4) == 1) then
    xb = 0.
    i = 2
  endif
  
  f1dt = 2*abs(px1(i)-xb) / (abs(pvn(i)) * sqrt( pvn(i)*pvn(i) - 2*an(i)*(px1(i)-xb)))	! Fractional time step when collision occurs

else

  ! Detect which wall was crossed first by comparing the fractional time step when collision occurs with each boundary
  if (pcol(1)==1 .AND. pcol(2)==1) then

    f1dtx = 2*abs(px1(1)-a) / (abs(pvn(1)) * sqrt( pvn(1)*pvn(1) - 2*an(1)*(px1(1)-a) ))
    f1dty = 2*abs(px1(2)-b) / (abs(pvn(2)) * sqrt( pvn(2)*pvn(2) - 2*an(2)*(px1(2)-b) ))

    if (f1dtx < f1dty) then
      f1dt = f1dtx
      pcol(2) = 0
    elseif (f1dty < f1dtx) then
      f1dt = f1dty
      pcol(1) = 0
    endif

  elseif (pcol(1)==1 .AND. pcol(4)==1) then

    f1dtx = 2*abs(px1(1)-a) / (abs(pvn(1)) * sqrt( pvn(1)*pvn(1) - 2*an(1)*(px1(1)-a) ))
    f1dty = 2*abs(px1(2)) / (abs(pvn(2)) * sqrt( pvn(2)*pvn(2) - 2*an(2)*px1(2) ))

    if (f1dtx < f1dty) then
      f1dt = f1dtx
      pcol(4) = 0
    elseif (f1dty < f1dtx) then
      f1dt = f1dty
      pcol(1) = 0
    endif

  elseif (pcol(3)==1 .AND. pcol(2)==1) then

    f1dtx = 2*abs(px1(1)) / (abs(pvn(1)) * sqrt( pvn(1)*pvn(1) - 2*an(1)*px1(1) ))
    f1dty = 2*abs(px1(2)-b) / (abs(pvn(2)) * sqrt( pvn(2)*pvn(2) - 2*an(2)*(px1(2)-b) ))

    if (f1dtx < f1dty) then
      f1dt = f1dtx
      pcol(2) = 0
    elseif (f1dty < f1dtx) then
      f1dt = f1dty
      pcol(3) = 0
    endif

  elseif (pcol(3)==1 .AND. pcol(4)==1) then

    f1dtx = 2*abs(px1(1)) / (abs(pvn(1)) * sqrt( pvn(1)*pvn(1) - 2*an(1)*px1(1) ))
    f1dty = 2*abs(px1(2)) / (abs(pvn(2)) * sqrt( pvn(2)*pvn(2) - 2*an(2)*px1(2) ))

    if (f1dtx < f1dty) then
      f1dt = f1dtx
      pcol(4) = 0
    elseif (f1dty < f1dtx) then
      f1dt = f1dty
      pcol(3) = 0
    endif

  endif

endif

f1 = f1dt/dt	! Fraction of time step when collision occurrs

if (isnan(f1)) then
  write(*,*) 'f1 is wrong. terminating.'
  call exit
endif

!!!!!! Particle Collection !!!!!!

! Particle velocity at time n+f1/2
pvf1d2(1) = pvn(1) + an(1) * f1dt * 0.5
pvf1d2(2) = pvn(2) + an(2) * f1dt * 0.5
pvf1d2(3) = pvn(3) + an(3) * f1dt * 0.5

! Particle position at time n+f1 (collision position)
pxf1(1) = px1(1) + pvf1d2(1) * f1dt
pxf1(2) = px1(2) + pvf1d2(2) * f1dt
pxf1(3) = px1(3) + pvf1d2(3) * f1dt

if (fields == 1 .OR. fields == 3) then

  ! Interpolate fields from grid to collision point
  call plocate(pxf1,d4x,x3,geo,pwall,Eii,Eip,Eri,Hii,Hip,Hri)

  ! Interpolate fields from nodes to particles
  wB1 = interB2particle(Hx1,Hy1,Hz1,Hii,Hip,Hri,pwall,sB,geo,px1(3))
  wB2 = interB2particle(Hx2,Hy2,Hz2,Hii,Hip,Hri,pwall,sB,geo,px1(3))
  wE2 = interE2particle(Ex2,Ey2,Ez2,Eii,Eip,Eri,pwall)

elseif (fields == 2) then

  ! Calculate fields using analytic formulas
  call pFIELDSanalytic(geo,d4x,pxf1,wave,sB,n-1,wB1,wE1)
  call pFIELDSanalytic(geo,d4x,pxf1,wave,sB,n,wB2,wE2)

endif

! Obtain fields at time n+f1/2-1/4
gradEn = 0.
gradBn = 0.
do j = 1,3
  do l = 1,3
    gradEn(j,l) = (wE2(j) - pE2(j))/(pxf1(l) - px1(l))
    gradBn(j,l) = (wB2(j)+wB1(j) - pB2(j)-pB1(j))/(2*(pxf1(l)-px1(l)))

  enddo
enddo

dp = (f1 - 0.5) * 0.5
pEf1 = 0.
pBf1 = 0.
do j = 1,3
  if ((pxf1(1) - px1(1))==0. .AND. (pxf1(2) - px1(2))==0.) then
    pEf1(j) = pE2(j) + ( pE2(j)-pE1(j) + dt*pvn(3)*gradEn(j,3) ) * dp
    pBf1(j) = pBn(j) + ( (pB2(j)-pB0(j))*0.5 + dt*pvn(3)*gradBn(j,3) ) * dp
  elseif ((pxf1(2) - px1(2)) == 0. .AND. (pxf1(3) - px1(3))==0.) then
    pEf1(j) = pE2(j) + ( pE2(j)-pE1(j) + dt*pvn(1)*gradEn(j,1) ) * dp
    pBf1(j) = pBn(j) + ( (pB2(j)-pB0(j))*0.5 + dt*pvn(1)*gradBn(j,1) ) * dp
  elseif ((pxf1(3) - px1(3)) == 0. .AND. (pxf1(1) - px1(1))==0.) then
    pEf1(j) = pE2(j) + ( pE2(j)-pE1(j) + dt*pvn(2)*gradEn(j,2) ) * dp
    pBf1(j) = pBn(j) + ( (pB2(j)-pB0(j))*0.5 + dt*pvn(2)*gradBn(j,2) ) * dp
  elseif ((pxf1(1) - px1(1)) == 0.) then
    pEf1(j) = pE2(j) + ( pE2(j)-pE1(j) + dt*( pvn(2)*gradEn(j,2) + pvn(3)*gradEn(j,3) ) ) * dp
    pBf1(j) = pBn(j) + ( (pB2(j)-pB0(j))*0.5 + dt*( pvn(2)*gradBn(j,2) + pvn(3)*gradBn(j,3) ) ) * dp
  elseif ((pxf1(2) - px1(2)) == 0.) then
    pEf1(j) = pE2(j) + ( pE2(j)-pE1(j) + dt*( pvn(1)*gradEn(j,1) + pvn(3)*gradEn(j,3) ) ) * dp
    pBf1(j) = pBn(j) + ( (pB2(j)-pB0(j))*0.5 + dt*( pvn(1)*gradBn(j,1) + pvn(3)*gradBn(j,3) ) ) * dp
  elseif ((pxf1(3) - px1(3)) == 0.) then
    pEf1(j) = pE2(j) + ( pE2(j)-pE1(j) + dt*( pvn(1)*gradEn(j,1) + pvn(2)*gradEn(j,2) ) ) * dp
    pBf1(j) = pBn(j) + ( (pB2(j)-pB0(j))*0.5 + dt*( pvn(1)*gradBn(j,1) + pvn(2)*gradBn(j,2) ) ) * dp
  else
    pEf1(j) = pE2(j) + ( pE2(j)-pE1(j) + dt*( pvn(1)*gradEn(j,1) + pvn(2)*gradEn(j,2) + pvn(3)*gradEn(j,3) ) ) * dp
    pBf1(j) = pBn(j) + ( (pB2(j)-pB0(j))*0.5 + dt*( pvn(1)*gradBn(j,1) + pvn(2)*gradBn(j,2) + pvn(3)*gradBn(j,3) ) ) * dp
  endif
enddo

! Obtain collision velocity at time n+f1 using the Boris solver    
 call boris_solver(dt*(f1+0.5),pEf1,pBf1,pvf1,pv1)
do k = 1,3
  if (isNaN(pvf1(k))) then
    write(*,*) 'Error in particle_wall.f90: pvf1 is wrong.'
    write(*,*) 'pxf1 = ', pxf1
    write(*,*) 'px1 = ', px1
    write(*,*) 'pvf1 = ', pvf1
    write(*,*) 'pv1 = ', pv1
    write(*,*) 'gradEn = ',gradEn(1,:)
    write(*,*) gradEn(2,:)
    write(*,*) gradEn(3,:)
    write(*,*) 'gradBn = ',gradBn(1,:)
    write(*,*) gradBn(2,:)
    write(*,*) gradBn(3,:)
    write(*,*) 'pEf1 = ', pEf1
    write(*,*) 'pBf1 = ', pBf1
    write(*,*) 'Terminating.'
    call exit
  endif
enddo

!!!!!! Determine number of emission electrons and their velocities !!!!!!
 call see(atype,pvf,s,pcol,pvf1,sey,SEYvsE)

if (.NOT. (allocated(pvf))) then
  write(*,*) 'error in subroutine see. pvf not allocated. TERMINATING.'
  call exit
endif

if (s > 1) then
  allocate(px2t(s,3))
  allocate(pv2t(s,3))
  px2t = 0.      ! Temporary array for storing position of emitted electrons
  pv2t = 0.      ! Temporary array for storing velocity of emitted electrons
elseif (s == 1 .OR. s == 0) then
  allocate(px2t(1,3))
  allocate(pv2t(1,3))
  px2t = 0.
  pv2t = 0.
endif

if (s>0) then
  do l = 1,s
    if (pvf(l,1)==0. .AND. pvf(l,2)==0. .AND. pvf(l,3)==0.) then
      write(*,*) 'stationary secondary'
      write(*,*) 's = ',s
      write(*,*) 'pvf = ',pvf(l,:)
    endif
  enddo
endif

!!!!!! Particle Ejection !!!!!!
if (s > 0) then
        
  o = 0                      ! Counter to keep track of secondary electrons that escape
  allocate(s3(s))            ! Store k indices of secondary electrons lost
  s3 = 0

  pBn = 0.5*(wB1 + wB2)
  dp = (1-f1) * dt       ! Useful factor  
  dtp = dt*(0.5-f1)      ! Useful factor
  
  do k = 1,s

    ! Particle acceleration at time n+1/4+3f/4
    an(1) = eme*( wE2(1) + pvf(k,2)*pBn(3) - pvf(k,3)*pBn(2) )
    an(2) = eme*( wE2(2) + pvf(k,3)*pBn(1) - pvf(k,1)*pBn(3) )
    an(3) = eme*( wE2(3) + pvf(k,1)*pBn(2) - pvf(k,2)*pBn(1) )

    ! Particle velocity at time n+1/2+f/2
    pvd2f(1) = pvf(k,1) + an(1) * dp * 0.5
    pvd2f(2) = pvf(k,2) + an(2) * dp * 0.5
    pvd2f(3) = pvf(k,3) + an(3) * dp * 0.5

    ! Particle position at time n+1 (ejection position). Temporary variable
    px2t(k,1) = pxf1(1) + pvd2f(1) * dp
    px2t(k,2) = pxf1(2) + pvd2f(2) * dp
    px2t(k,3) = pxf1(3) + pvd2f(3) * dp

    ! Check if new location is inside physical waveguide
    if ((px2t(k,1)>0. .AND. px2t(k,1)<a) .AND. (px2t(k,2)>0. .AND. px2t(k,2)<b) &
            .AND. (px2t(k,3)>0. .AND. px2t(k,3)<Lw)) then
            
      if (fields == 1 .OR. fields == 3) then
        ! Interpolate fields from grid to collision point
	wB0 = interB2particle(Hx0,Hy0,Hz0,Hii,Hip,Hri,pwall,sB,geo,px1(3))
	wE1 = interE2particle(Ex1,Ey1,Ez1,Eii,Eip,Eri,pwall)

	! Interpolate fields from grid to ejection point
	call plocate(px2t(k,:),d4x,x3,geo,pwall,Eii,Eip,Eri,Hii,Hip,Hri)

	! Interpolate fields from nodes to particles
	eB1 = interB2particle(Hx1,Hy1,Hz1,Hii,Hip,Hri,pwall,sB,geo,px1(3))
	eB2 = interB2particle(Hx2,Hy2,Hz2,Hii,Hip,Hri,pwall,sB,geo,px1(3))
	eE2 = interE2particle(Ex2,Ey2,Ez2,Eii,Eip,Eri,pwall)

      elseif (fields == 2) then
        ! Calculate field at collision point for 2 iterations ago
	call pFIELDSanalytic(geo,d4x,pxf1,wave,sB,n-2,wB0)
            
        ! Interpolate fields from grid to ejection point
	call pFIELDSanalytic(geo,d4x,px2t(k,:),wave,sB,n-1,eB1)
	call pFIELDSanalytic(geo,d4x,px2t(k,:),wave,sB,n,eB2,eE2)
      endif
            
      ! Fields felt by the particle at time n+f/2+1/4
      gradEn = 0.	! Gradient of E-field felt at ejection position at time n
      gradBn = 0.	! Gradient of B-field felt at ejection position at time n
      do j = 1,3
        do l = 1,3
          gradEn(j,l) = (eE2(j) - wE2(j))/(px2t(k,l) - pxf1(l))
          gradBn(j,l) = (eB2(j)+eB1(j) - wB2(j)-wB1(j))/(2*(px2t(k,l)-pxf1(l)))

	  if (px2t(k,l)-pxf1(l)==0.) then
	    write(*,*) 'difference is zero'
	    write(*,*) 'log of px2t = ',log(px2t(k,l))
	    write(*,*) 'log of pxf1 = ',log(pxf1(l))
	    write(*,*) 'diff of terms*10e10 = ',px2t(k,l)*1e10-pxf1(l)*1e10
	  endif
        enddo
      enddo

      pEf = 0.
      pBf = 0.
      do j = 1,3
        if ((px2t(k,1) - pxf1(1)) == 0. .AND. (px2t(k,2) - pxf1(2)) == 0.) then
          pEf(j) = wE1(j) + ( (0.5+f1) * (wE2(j) - wE1(j)) + dtp*pvf(k,3)*gradEn(j,3) ) * 0.5
          pBf(j) = pBn(j) + ( (0.5+f1) * (wB2(j) - wB0(j)) * 0.5 + dtp*pvf(k,3)*gradBn(j,3) ) * 0.5
        elseif ((px2t(k,2) - pxf1(2)) == 0. .AND. (px2t(k,3) - pxf1(3)) == 0.) then
          pEf(j) = wE1(j) + ( (0.5+f1) * (wE2(j) - wE1(j)) + dtp*pvf(k,1)*gradEn(j,1) ) * 0.5
          pBf(j) = pBn(j) + ( (0.5+f1) * (wB2(j) - wB0(j)) * 0.5 + dtp*pvf(k,1)*gradBn(j,1) ) * 0.5
        elseif ((px2t(k,3) - pxf1(3)) == 0. .AND. (px2t(k,1) - pxf1(1)) == 0.) then
          pEf(j) = wE1(j) + ( (0.5+f1) * (wE2(j) - wE1(j)) + dtp*pvf(k,2)*gradEn(j,2) ) * 0.5
          pBf(j) = pBn(j) + ( (0.5+f1) * (wB2(j) - wB0(j)) * 0.5 + dtp*pvf(k,2)*gradBn(j,2) ) * 0.5
        elseif ((px2t(k,1) - pxf1(1)) == 0.) then
          pEf(j) = wE1(j) + ( (0.5+f1) * (wE2(j) - wE1(j)) + dtp*( pvf(k,2)*gradEn(j,2) + &
			pvf(k,3)*gradEn(j,3) ) ) * 0.5
          pBf(j) = pBn(j) + ( (0.5+f1) * (wB2(j) - wB0(j)) * 0.5 + dtp*( pvf(k,2)*gradBn(j,2) + &
			pvf(k,3)*gradBn(j,3) ) ) * 0.5
        elseif ((px2t(k,2) - pxf1(2)) == 0.) then
          pEf(j) = wE1(j) + ( (0.5+f1) * (wE2(j) - wE1(j)) + dtp*( pvf(k,1)*gradEn(j,1) + &
			pvf(k,3)*gradEn(j,3) ) ) * 0.5
          pBf(j) = pBn(j) + ( (0.5+f1) * (wB2(j) - wB0(j)) * 0.5 + dtp*( pvf(k,1)*gradBn(j,1) + &
			pvf(k,3)*gradBn(j,3) ) ) * 0.5
        elseif ((px2t(k,3) - pxf1(3)) == 0.) then
          pEf(j) = wE1(j) + ( (0.5+f1) * (wE2(j) - wE1(j)) + dtp*( pvf(k,1)*gradEn(j,1) + &
			pvf(k,2)*gradEn(j,2) ) ) * 0.5
          pBf(j) = pBn(j) + ( (0.5+f1) * (wB2(j) - wB0(j)) * 0.5 + dtp*( pvf(k,1)*gradBn(j,1) + &
			pvf(k,2)*gradBn(j,2) ) ) * 0.5
        else
          pEf(j) = wE1(j) + ( (0.5+f1) * (wE2(j) - wE1(j)) + dtp*( pvf(k,1)*gradEn(j,1) + &
			pvf(k,2)*gradEn(j,2) + pvf(k,3)*gradEn(j,3) ) ) * 0.5
          pBf(j) = pBn(j) + ( (0.5+f1) * (wB2(j) - wB0(j)) * 0.5 + dtp*( pvf(k,1)*gradBn(j,1) + &
			pvf(k,2)*gradBn(j,2) + pvf(k,3)*gradBn(j,3) ) ) * 0.5
        endif
      enddo

      ! Obtain emission velocity at time n+1/2 using the Boris solver
      call boris_solver(dtp,pEf,pBf,pv2t(k,:),pvf(k,:))

      if (px2t(k,1)==0. .AND. px2t(k,2)==0. .AND. px2t(k,3)==0.) then
        write(*,*) 'stationary secondary'
        write(*,*) 's = ',s
        write(*,*) 'o = ',o
        write(*,*) 'px2t = ',px2t(k,:)
!        write(*,*) 'pv2t = ',pv2t
      endif

    else
           
      ! Count and record which secondaries escape the domain 
      o = o + 1
      s3(o) = k
            
    endif

  enddo

  ! Get rid of secondaries that escaped the domain
  if (o > 0 .AND. s > 1 .AND. o/=s) then 

    allocate(px2temp(s-o,3))	! Array where location of active secondaries will be temporarily stored
    allocate(pv2temp(s-o,3))	! Array where velocities of active secondaries will be temporarily stored

    i = 1
    j = 1
    do l = 1,s
        
        if (l /= s3(i)) then
          px2temp(j,:) = px2t(l,:)
          pv2temp(j,:) = pv2t(l,:)
	  j = j + 1
	elseif (l == s3(i)) then
	  i = i + 1
        endif

    enddo
    
    s = s - o

    if (s /= j-1) then
	write(*,*) 'counting is wrong'
    endif

    deallocate(px2t)
    deallocate(pv2t)

    allocate(px2t(s,3))
    allocate(pv2t(s,3))

    px2t = px2temp
    pv2t = pv2temp

    do l = 1,s
      do k = 1,3
        if (isnan(px2t(l,k))) then
          write(*,*) 'px2t is wrong. terminating.'
          call exit
        endif
        if (isnan(pv2t(l,k))) then
          write(*,*) 'pv2t is wrong. terminating.'
          call exit
        endif
      enddo
    enddo

    deallocate(px2temp)
    deallocate(pv2temp)

  elseif (o == s) then

    s = 0

  endif

  deallocate(s3)
        
endif

if (o > 0 .AND. s > 0) then

!  write(*,*) 's = ',s
!  write(*,*) 'px2t = '
!  do k = 1,s
!    write(*,*) (px2t(k,j), j=1,3)
!  enddo
!  write(*,*) 'pv2t = '
!  do k = 1,s
!    write(*,*) (pv2t(k,j), j=1,3)
!  enddo
!  call exit

endif

!write(*,*) 's = ',s


deallocate(pvf)

end subroutine particlewall

end module particle_wall
