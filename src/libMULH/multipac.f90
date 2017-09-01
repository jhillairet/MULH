module multipac

contains
subroutine multipactor(atype,dirname,P_I,outputs,pro,simno,complete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!			Multipactor
!
! This routine contains most of the skeleton holding the numerics
! together.
!
! ARGUMENTS
!   - atype: integer, type of analysis being performed
!   - dirname: string containing the directory name where results are written
!   - P_I: integer, input power of RF wave
!   - outputs: defined type, indicate which things to write/save/output
!   - pro: integer flag indicating whether to proceed or not (due to error)
!   - simno: integer, simulation number
!   - complete: integer flag indicating whether multipactor was found or not
!
!

use constants
use data_type
use mod_Vaughan
use prep_fields
use particle_load
use reshape_array
use particle_launch
use p_locate
use inter2part
use FIELDScalc
use boris
use particle_wall
use MP_detector
use secondary_launch

implicit none

! Subroutine arguments
 character(42), intent(in) :: dirname
integer, intent(in) :: atype, P_I, simno
integer, intent(inout) :: complete
type(output_s), intent(in) :: outputs
logical, intent(inout) :: pro
!real(long), allocatable, intent(inout) :: pvm(:,:)


open (unit=11,file=trim(dirname)//'/scenario.txt',status='old')

read(11,*) a
read(11,*) b
read(11,*) lw
read(11,*) NOC_lambda
read(11,*) CLfactor
read(11,*) ramp
read(11,*) NOC_PML
read(11,*) R_max
read(11,*) m_PML
read(11,*) f_I
read(11,*) fields
read(11,*) sB(1)
read(11,*) sB(2)
read(11,*) sB(3)
read(11,*) Np
read(11,*) vth
read(11,*) vra
read(11,*) px_i
read(11,*) fmax
read(11,*) nsamplev
read(11,*) launch
read(11,*) seec
read(11,*) mat
read(11,*) ReRr
read(11,*) E1
read(11,*) Emax
read(11,*) delta_b
read(11,*) deltamax
read(11,*) p_n
read(11,*) kse
read(11,*) ks
read(11,*) Eom
read(11,*) Ns
read(11,*) a_lara
read(11,*) z_lara

 close (unit=11)

psave = outputs%psave

!!!!!!!!!!!!! Parameters derived from user inputs !!!!!!!!!!
eta = mu_0*c	! Wave impedance
w = 2*pi*f_I   	! Angular frequency of incident wave
lambda = c/f_I 	! Wavelength of incident wave
T = lambda/c    ! Period
beta = sqrt( mu_0*eps_0*w*w - (pi/b)**2 )	! Perpendicular wave number

dz = lambda/NOC_lambda

NOC(3) = ceiling(lw/dz)             ! Number of Cells along z (in physical waveguide area)

dz = lw/NOC(3)			! Size of cells along z

d = dz * NOC_PML               	! Depth of PML in meters

NOC(3) = ceiling((real(lw,4)+real(d,4))/real(dz,4))	! Number of Cells along z
NOC(1) = NOC_lambda           ! Number of Cells along x
NOC(2) = NOC_lambda             ! Number of Cells along y

dx = a/NOC(1)                  	! Size of cells along x
dy = b/NOC(2)                  	! Size of cells along y
dz = lw/(NOC(3)-NOC_PML)     	! Size of cells along z

d = dz*NOC_PML			! Adjust d with new dz

! Time step using Courant-Levy condition
gamm = (( 1./(dx*dx) + 1./(dy*dy) + 1./(dz*dz) )**(-0.5))/c
dt = CLfactor*gamm
ndt = floor(T/dt)		! Use integer number of steps per period
dt = T/ndt

if (dt > gamm) then
   print*, 'Reduce Courant-Levy condition factor'
   stop 1
endif

E0 = sqrt(4*eta*P_I/(a*sqrt(b*b-(lambda*lambda/4))))	! Amplitude of waves for design power ooutput (V/m)
pstart = INT(T/dt) * ramp			! Iteration when particle dyamics starts
Nt = Np + Np*Ns			! Total number of particles tracked

! Indentify antenna and allocate plasma center radius and antenna radius
if (INT(a*1000)==8 .AND. INT(b*1000)==70) then
    r0 = 2.96		! Radial position of plasma center
    rLH = 3.93		! Radial position (in tokamak) of grill mouth
    if(n==0) print*, 'ToreSupra_C3 detected'
else if (INT(a*1000)==9 .AND. INT(b*1000)==72) then
    r0 = 2.96
    rLH = 3.93
    if(n==0) print*, 'JET detected'
end if

! Modify Vaughan model to match E1
if ((seec == 1 .OR. seec == 3) .AND. deltamax > 1) then
  call modVaughan(E_0,E_0p,E1,Emax,deltamax)
elseif (deltamax < 1) then
  E_0 = E1
  E_0p = E1 - 1
endif

! Package info
sey = sey_t(Emax,seec,ReRr,E1,E_0,E_0p,deltamax,kse,ks,Eom,p_n,a_lara,delta_b,z_lara)	! Particle-wall
d4x = (/ dx, dy, dz, dt /) 	! Space-time steps
wave = (/ w, beta, E0, lambda/)	! Electromagnetic wave
geo = geo_t(a, b, lw, R_max, d, NOC_PML, NOC, m_PML) 	! Geometry
seed = seed_t(Np, Nt, vth, fmax, nsamplev, vra)	! Seed particles
! Static DC magnetic fields and plasma and LH antenna radial location
sB(4) = r0
sB(5) = rLH

!!!!!!!!!!!!!!!!!! Populate arrays of spatial coordinates !!!!!!!!!!!!!!!
allocate (x(NOC(1)), y(NOC(2)), z(NOC(3)))

do i = 1,NOC(1)
  x(i) = (i-0.5)*dx
enddo

do i = 1,NOC(2)
  y(i) = (i-0.5)*dy
enddo

do i = 1,NOC(3)
  z(i) = (i-0.5)*dz
enddo

allocate (x3(sum(NOC)))

x3 = (/ x, y, z /)	! Pack x,y,z arrays

!!!!!!!!!!!!! Storing E-field at all x,y,z nodes !!!!!!!!!!!!!!
!allocate (Ext(NOC(1),NOC(2),NOC(3)))
!do i = 1,NOC(1)
!  do j = 1,NOC(2)
!    do k = 1,NOC(3)

!      call pFIELDSanalytic(geo,d4x,(/x(i),y(j),z(k)/),wave,sB,0,pB2,pE2)
!      Ext(i,j,k) = pE2(1)

!    enddo
!  enddo
!enddo

! store y-z plane at x(1)
!open(unit=63,file=trim(dirname)//'/ext.txt',status='unknown')
!do k = 1,NOC(3)
!do j = 1,NOC(2)
!  do i = 1,NOC(1)
!    write(63,*) (Ext(1,j,k), k=1,NOC(3))
!  enddo
!enddo
! close(unit=63)
psave = 0
!!!!!! Save setup !!!!!!
open(unit=13,file=trim(dirname)//'/setup.txt',status='unknown')
if (psave /= 0) then
  write(13,10) it,pstart,psave,NOC,NOC_PML,dt,dx,dy,dz,a,b,f_I,P_I,E1,Emax,delta_b,deltamax
else
  write(13,11) pstart,psave,NOC,NOC_PML,dt,dx,dy,dz,a,b,f_I,P_I,E1,Emax,delta_b,deltamax
endif
 close(unit=13)

!open(unit=63,file='ReRr.txt',status='unknown',position='append')
!write(63,*) P_I
!close(unit=63) 


!!!!!!!!!!!!!!!!! Prepare fields from external solver !!!!!!!!!!!!!!
if (fields == 1 .OR. fields == 3) then
  ! E field components
  allocate (Ex2(NOC(1),NOC(2)+1,NOC(3)+1))
  allocate (Ey2(NOC(1)+1,NOC(2),NOC(3)+1))
  allocate (Ez2(NOC(1)+1,NOC(2)+1,NOC(3)))
  Ex2 = 0.
  Ey2 = 0.
  Ez2 = 0.
  allocate (Ex1(NOC(1),NOC(2)+1,NOC(3)+1))
  allocate (Ey1(NOC(1)+1,NOC(2),NOC(3)+1))
  allocate (Ez1(NOC(1)+1,NOC(2)+1,NOC(3)))
  Ex1 = 0.
  Ey1 = 0.
  Ez1 = 0.

  ! H field components
  allocate (Hx2(NOC(1)+1,NOC(2),NOC(3)))
  allocate (Hy2(NOC(1),NOC(2)+1,NOC(3)))
  allocate (Hz2(NOC(1),NOC(2),NOC(3)+1))
  Hx2 = 0.
  Hy2 = 0.
  Hz2 = 0.
  allocate (Hx1(NOC(1)+1,NOC(2),NOC(3)))
  allocate (Hy1(NOC(1),NOC(2)+1,NOC(3)))
  allocate (Hz1(NOC(1),NOC(2),NOC(3)+1))
  Hx1 = 0.
  Hy1 = 0.
  Hz1 = 0.
  allocate (Hx0(NOC(1)+1,NOC(2),NOC(3)))
  allocate (Hy0(NOC(1),NOC(2)+1,NOC(3)))
  allocate (Hz0(NOC(1),NOC(2),NOC(3)+1))
  Hx0 = 0.
  Hy0 = 0.
  Hz0 = 0.

  if (fields == 1) then
    call prep_fdtd_fields(geo,wave,d4x,x3,CE,DH,CEsigz,DHsigz, &
	ExPML,EyPML,EzPML,HxPML,HyPML,HzPML,Ex2,Hy2,Hz2,Jx,Jy,Jz)
  else
    call prep_ex_fields(geo,wave,ndt,d4x,x3,P_I,sB,Exn,Eyn,Ezn,Hxn,Hyn,Hzn)
  endif

endif


!!!!!!!!!!!!!!!!! Iinitialize particle arrays and variables !!!!!!!!!!!!!!!!!!!

allocate (px1(Nt,3), px2(Nt,3), pv1(Nt,3), pv2(Nt,3), pstat(Nt))
allocate (pphases(Np,2), phi0s(Np,2))
allocate (pactts(1))
allocate (pactps(1))

pcol = 0	! Flags, one element for each wall particle collided with
part = 0.	! Used to pack all info about a particle
pB0 = 0.	! B felt by particle @ time n-3/2
pB1 = 0.	! B felt by particle @ time n-1/2
pB2 = 0.	! B felt by particle @ time n+1/2
pE1 = 0.	! E felt by particle @ time n-1
pE2 = 0.	! E felt by particle @ time n
px1 = 0.	! Particle position @ time n
px2 = 0.	! Particle position @ time n+1
pv1 = 0.	! Particle velocity @ time n-1/2
pv2 = 0.	! Particle velocity @ time n+1/2
pstat = 0	! Status of particles: 0= dormant, 1=launched, 2=collided, 3=escaped, 4=absorbed
pphases = 0.	! Phases each particle will be launched @ (1st column) & flag if launched @ that phase already (2nd column)
phi0s = 0.	! Variable that stores phases at which particles were launched/activated.
pwall = 0	! Flag to identify if particle is between node and wall along a given direction
ps_c = 0	! Counter keeping track of number of secondary electrons activated
lost = 0	! Number of escape particles (1st column) and number of absorbed particles (2nd column)
launched = 0	! Flag to indicate all initial particles have been launched
complete = 0	! Flag what was the result of completing the simulation
pactc = 1	! Count entries in array of time vs. active particles which will be stored
pact = 0	! Number of active particles at each iteration
pactts = 0.	! Times at which the number of particles is saved in pactps. Used for writing results in MULH_time_results.
pactps = 0	! Number of particles at pactts times. Used for writing results in MULH_time_results
full = .FALSE.

! Open files used for loading random values and controlling stochasticity in order to produce identical runs
if (atype == 8 .OR. atype == 9) then
  open(unit=14,file='../data/px2.txt',status='old',action='read')
  open(unit=15,file='../data/pv2.txt',status='old',action='read')
  open(unit=16,file='../data/phases0.txt',status='old',action='read')
  rewind(unit=14)
  rewind(unit=15)
  rewind(unit=16)

  if (atype == 9) then
    open(unit=17,file='../data/srand.txt',status='old',action='read')
    open(unit=18,file='../data/drand.txt',status='old',action='read')
    open(unit=19,file='../data/irand.txt',status='old',action='read')
    open(unit=20,file='../data/secrand.txt',status='old',action='read')
  else
    open(unit=21,file='../data/rstock.txt',status='old',action='read')
  endif

endif

! Load seed electrons
 call particleload(atype,d4x,geo,seed,px_i,fields,px2,pv2)
if (atype == 8 .OR. atype == 9) then
 close(unit=14)
 close(unit=15)
endif

! Store initial particle locations
if (psave /= 0) then
  makedir = 'mkdir -p '//trim(dirname)// "/" // 'px'
  call system(makedir)
  open(unit=3,file=trim(dirname)//'/px/'//'px0.txt',status='new')
  do i = 1,Np
    write(3,*) (px2(i,j), j=1,3)
  enddo
  close(unit=3)
  makedir = 'mkdir -p '//trim(dirname)// "/" // 'pv'
  call system(makedir)
  open(unit=3,file=trim(dirname)//'/pv/'//'pv0.txt',status='new')
  do i = 1,Np
    write(3,*) (pv2(i,j), j=1,3)
  enddo
endif

!!!!!!!!!!!!!!!!!!!!!!! MAIN TIME LOOP !!!!!!!!!!!!!!!!!!!!!!
n = 0
write(*,'(A21,F8.3,A2)') ' Current power used: ',real(P_I)/1000,'kW'
do while (complete == 0)

  n = n+1
  time = n*dt

  if (fields == 1 .OR. fields == 3) then
    if (n > 1) then
      if (n > 2) then
	! Store fields from two iterations ago
	Hx0 = Hx1
	Hy0 = Hy1
	Hz0 = Hz1
      endif

      ! Store fields from previous iteration
      Hx1 = Hx2
      Hy1 = Hy2
      Hz1 = Hz2
      Ex1 = Ex2
      Ey1 = Ey2
      Ez1 = Ez2
    endif
    
    ! Updates fields
    if (fields == 3) then
      if (mod(n,ndt) == 0) then
        Ex2 = Exn(:,:,:,ndt)
        Ey2 = Eyn(:,:,:,ndt)
        Ez2 = Ezn(:,:,:,ndt)
        Hx2 = Hxn(:,:,:,ndt)
        Hy2 = Hyn(:,:,:,ndt)
        Hz2 = Hzn(:,:,:,ndt)
      else
        Ex2 = Exn(:,:,:,mod(n,ndt))
        Ey2 = Eyn(:,:,:,mod(n,ndt))
        Ez2 = Ezn(:,:,:,mod(n,ndt))
        Hx2 = Hxn(:,:,:,mod(n,ndt))
        Hy2 = Hyn(:,:,:,mod(n,ndt))
        Hz2 = Hzn(:,:,:,mod(n,ndt))
      endif
    else
      call fdtd_fields(geo,wave,x3,d4x,ramp,CE,DH,CEsigz,DHsigz, &
	n,Ex2,Ey2,Ez2,Hx2,Hy2,Hz2,ExPML,EyPML,EzPML,HxPML,HyPML,HzPML,Jx,Jy,Jz)
    endif
  endif
  
  !!!!!!!!!!!!!!!!!!!!! PARTICLES !!!!!!!!!!!!!!!!!!!!!
  if (n > pstart) then
    
    ! Indicate progress
    !if (mod(n,10000) == 0) then
     ! call date_and_time(TIME=runtime)
      !write(*,8) 'Time:',runtime,'| SimTime:',time,'s | Electrons in waveguide:',pact
    !endif

    px1 = px2	! Store previous particle positions
    pv1 = pv2	! Store previous particle velocities
    ps_c0 = 0	! Counter keeping track of secondaries activated within one time iteration
    pact0 = pact

    ! Activate particles only if they haven't all been activated already
    if (launched == 0) then

	! IMPORTANT NOTE: This code assumes that all seed particles are
        ! launched before any seed collides with a wall. If not true
        ! then need to change the way particles are added/absorbed at
        ! earlier times and possibly bubblesort.
	call particlelaunch(atype,geo,d4x,pact,pphases,pstart,pstat,phi0s,Np,wave,px1,pv1,launch,n)

	if (pact == Np) then
	  launched = 1
	endif
    endif

    do m = 1,pact

      if (pstat(m)==1) then
	
	pcol = 0	

	if (fields == 1 .OR. fields == 3) then

	  ! Locate particle in E and H grids
	  call plocate(px1(m,:),d4x,x3,geo,pwall,Eii,Eip,Eri,Hii,Hip,Hri)

	  ! Interpolate fields from nodes to particles
	  pB1 = interB2particle(Hx1,Hy1,Hz1,Hii,Hip,Hri,pwall,sB,geo,px1(m,3))
	  pB2 = interB2particle(Hx2,Hy2,Hz2,Hii,Hip,Hri,pwall,sB,geo,px1(m,3))
	  pE2 = interE2particle(Ex2,Ey2,Ez2,Eii,Eip,Eri,pwall)

	else if (fields == 2) then
	  ! Calculate analytic fields inside rectangular waveguide at times n and n-1
	  call pFIELDSanalytic(geo,d4x,px1(m,:),wave,sB,n,pB2,pE2)
	  call pFIELDSanalytic(geo,d4x,px1(m,:),wave,sB,n-1,pB1,pE1)

	endif

	! Step particles with Boris explicit trapezoidal scheme
	call boris_solver(dt,pE2,(pB1+pB2)*0.5,pv2(m,:),pv1(m,:),px2(m,:),px1(m,:))

	! Detect if particle underwent collision with wall or escaped
        if (px2(m,1) <= 0.) then
          pcol(3) = 1
          pstat(m) = 2
        endif
        if (px2(m,1) >= a) then
          pcol(1) = 1
          pstat(m) = 2
        endif
        if (px2(m,2) <= 0.) then
          pcol(4) = 1
          pstat(m) = 2
        endif
        if (px2(m,2) >= b) then
          pcol(2) = 1
          pstat(m) = 2
        endif
        if (px2(m,3) <= 0.) then
          pstat(m) = 3
        endif
        if (px2(m,3) >= lw) then
          pstat(m) = 3
        endif

	! Particles interact with walls if boundary was crossed
        if (pstat(m) == 2) then

	  if (fields == 1 .OR. fields == 3) then
            ! Interpolate additional fields needed
	    pB0 = interB2particle(Hx0,Hy0,Hz0,Hii,Hip,Hri,pwall,sB,geo,px1(m,3))
	    pE1 = interE2particle(Ex1,Ey1,Ez1,Eii,Eip,Eri,pwall)

          elseif (fields == 2) then      
	    call pFIELDSanalytic(geo,d4x,px1(m,:),wave,sB,n-2,pB0)
                      
          endif
	  
	  ! Pack all info about the particle in one variable
          particle = particle_t(pcol,px1(m,:),pv1(m,:),pv2(m,:),pB0,pB1,pB2,pE1,pE2)

	  if (fields == 1 .OR. fields == 3) then               
            call particlewall(px2t,pv2t,s,d4x,x3,geo,particle,sey,fields,wave,n,sB,SEYvsE,atype, &
		  Hx0,Hy0,Hz0,Hx1,Hy1,Hz1,Ex1,Ey1,Ez1,Hx2,Hy2,Hz2,Ex2,Ey2,Ez2)
          elseif (fields == 2) then
            call particlewall(px2t,pv2t,s,d4x,x3,geo,particle,sey,fields,wave,n,sB,SEYvsE,atype)
          endif

	  call activate_secondary(px2,pv2,px2t,pv2t,pstat,Np,Nt,pact,full,s,m,lost)

	  deallocate (px2t, pv2t)

	endif	!!!!!!!!!! end of pstat==2 if

      endif	!!!!!!!!!!! end of pstat==1 if

      if (full .AND. sum(pstat, MASK = pstat .EQ. 3)==0 .AND. sum(pstat, MASK = pstat .EQ. 4)==0) then
	write(*,*) 'Condition if pas comprise. Voir multipac.f90 ligne 498' !pact = Nt
	exit
      endif

      if (pstat(m) == 1 .AND. (px2(m,1) < 0 .OR. px2(m,1) > a)) then
        write(*,*) 'Error. Terminating.'
        call exit
      endif
      
      do l = 1,3      
        if ( isNaN(pv2(m,l)) ) then
	  write(*,*) 'pv1(m,:) = ', pv1(m,:)
	  write(*,*) 'pv2(m,:) = ', pv2(m,:)
          write(*,*) 'NaN. Terminating'
          call exit
        endif
      enddo

    enddo	!!!!!!!!!!! end of m=1,pact do loop

    ! Delete escape/absorbed particles by shifting arrays upwards
    do m = 1,pact
      if (pstat(m) == 3 .OR. pstat(m) == 4) then
                
        if (pstat(m) == 3) then
          lost(1) = lost(1) + 1
        else
          lost(2) = lost(2) + 1
        endif
        
        if (m < Nt) then        
          do i = m,Nt-1
            px2(i,:) = px2(i+1,:)
            pv2(i,:) = pv2(i+1,:)
            pstat(i) = pstat(i+1)
          enddo
	endif
        px2(Nt,:) = 0.
        pv2(Nt,:) = 0.
        pstat(Nt) = 0

      endif

    enddo

    pact = sum(pstat, MASK = pstat .EQ. 1)  
    
    !!!!!! Other diagnoses !!!!!!
    ! Store number of particles active vs. time
    if (mod(n,50) == 0 .OR. pact == Nt) then
      pactc = pactc + 1
      call add21Darray(pactts,1)		! Add one more elemet to pactts
      call add21DarrayI(pactps,1)		! Add one more elemet to pactps
      pactts(pactc) = time+dt
      pactps(pactc) = pact
    endif
        
    ! Store average energy along each dimension for each power run
    pvm(1) = mean(0.5*me*(pv2(:,1)*pv2(:,1))/(-e))
    pvm(2) = mean(0.5*me*(pv2(:,2)*pv2(:,2))/(-e))
    pvm(3) = mean(0.5*me*(pv2(:,3)*pv2(:,3))/(-e))
    write(24,'(F20.16,1X,F20.16,1X,F20.16)') pvm(1), pvm(2), pvm(3)

    ! Determine if MP will develop or not2                                                                                                   fields

    if (n >= pstart+50) then     ! Because we need pacts
      complete = MPdetector(pact,pactts,pactps,Np,Nt,pstart,sB,n,time,pstat)
    endif

    ! MULH time results. 4 columns: Simulation number, Power, time, active particles
    if (mod(n,50)==0 .OR. complete /= 0) then
      write(4,12) simno, P_I, pactts(pactc), pactps(pactc)
    endif

  endif		!!!!!!!!!!!!!! end of n>pstart if

  ! Save particle positions to files
  if (psave /= 0 .AND. (mod(n,psave) == 0 .OR. complete/=0)) then
    write(itn,'(I8)') n		! Convert n into a string
    itn = adjustl(itn)
    open(unit=7,file=trim(dirname)//'/px/'//'px'//trim(itn)//'.txt',status='new')
    do i = 1,Nt
      if (pstat(i) == 1) write(7,*) (px2(i,j), j=1,3)
    enddo
    close(unit=7)
    open(unit=8,file=trim(dirname)//'/pv/'//'pv'//trim(itn)//'.txt',status='new')
    do i = 1,Nt
      if (pstat(i) == 1) write(8,*) (pv2(i,j), j=1,3)
    enddo
    close(unit=8)
    it = n
  endif
  
enddo		!!!!!!!!!!!! end of do while complete

if (atype == 8 .OR. atype == 9) then
  close(unit=16)

  if (atype == 9) then
    close(unit=17)
    close(unit=18)
    close(unit=19)
    close(unit=20)
  else
    close(unit=21)
  endif

endif

!open(unit=10,file=trim(dirname)//'/SEYvsE.txt',status='unknown')
!do i = 1,size(SEYvsE,3)
!  write(10,13) SEYvsE(1,1,i),SEYvsE(1,2,i),SEYvsE(2,1,i), SEYvsE(2,2,i),SEYvsE(3,1,i),SEYvsE(3,2,i),SEYvsE(4,1,i),SEYvsE(4,2,i)
!enddo
 !close(unit=10)

!!!!!! Save setup !!!!!!
open(unit=13,file=trim(dirname)//'/setup.txt',status='unknown')
if (psave /= 0) then
  write(13,10) it,pstart,psave,NOC,NOC_PML,dt,dx,dy,dz,a,b,f_I,P_I,E1,Emax,delta_b,deltamax
else
  write(13,11) pstart,psave,NOC,NOC_PML,dt,dx,dy,dz,a,b,f_I,P_I,E1,Emax,delta_b,deltamax
endif
 close(unit=13)

!!!!!!!!!!!!!!!!!!!!!! DEALLOCATE ARRAYS !!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate (x,y,z,x3,px1,px2,pv1,pv2,pstat,pphases,phi0s,pactts,pactps)

if (fields == 1) deallocate (ExPML,EyPML,EzPML,Ex2,Ey2,Ez2,Ex1,Ey1,Ez1,HxPML,HyPML,HzPML,Hx2,Hy2,Hz2,Hx1,Hy1,Hz1,Hx0,Hy0,Hz0,&
			CE,Dh,CEsigz,DHsigz,Jx,Jy,Jz)
if (fields == 3) deallocate (Exn,Eyn,Ezn,Ex2,Ey2,Ez2,Ex1,Ey1,Ez1,Hxn,Hyn,Hzn,Hx2,Hy2,Hz2,Hx1,Hy1,Hz1,Hx0,Hy0,Hz0)

1 format(1X,A16,1X,I1,A2)
2 format(1X,A17,1X,I2,A2)
3 format(1X,A17,1X,I3,A2)
4 format(1X,A20,1X,F4.2,A3)
5 format(3X,A4,1X,I3,A2)
6 format(3X,A6,1X,I3,A2)
7 format(3X,A11,1X,F4.2)
8 format(3X,A5,1X,A10,1X,A8,1X,ES9.3E2,A27,1X,I3)
9 format(I2.1,1X,I6,ES17.15E2,1X,I3)
10 format(I8,1X,I3,1X,I3,1X,I3,1X,I3,1X,I3,1X,I2,1X,ES12.6E2,1X,F8.6,1X,F8.6,1X,F8.6,1X,F8.6,1X, &
F8.6,1X,ES10.4E2,1X,I7,1X,I3,1X,I3,1X,F5.3,1X,F5.3)
11 format(I3,1X,I3,1X,I3,1X,I3,1X,I3,1X,I2,1X,ES12.6E2,1X,F8.6,1X,F8.6,1X,F8.6,1X,F8.6,1X, &
F8.6,1X,ES10.4E2,1X,I7,1X,I3,1X,I3,1X,F5.3,1X,F5.3)
12 format(I2,1X,I7,1X,ES24.17E3,1X,I5)
13 format(F10.6,1X,F10.6,1X,F10.6,1X,F10.6,1X,F10.6,1X,F10.6,1X,F10.6,1X,F10.6)

end subroutine multipactor

end module multipac

