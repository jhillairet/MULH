module particle_launch

contains
subroutine particlelaunch(atype,geo,d4x,pact,phases,pstart,pstat,phi0s,Np,wave,px,pv,method,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!				PARTICLE LAUNCH
!
! Launch (activate) particles at different phases so that there isn't a bias towards a resonant
! or unresonant phase.
! method is used so that different ways of launching can be programmed in the future. At the moment
! particles are just launched at random phases.
! An option to read launch phases from a file is included in case identical runs are desired for
! debugging and such.
!
! INPUTS:
!   - atype: integer indicating what type analysis is being performed
!   - geo: user defined type array with info regarding geometry of waveguide and PML
!   - d4x: real array with the discretization steps of the 4 dimensions
!   - pact: integer number of active particles
!   - phases: real array, phases at which particles are launched (2nd column: whether particle at that phase has been launched)
!   - pstart: integer, when particles begin to launch
!   - pstat: integer array indicating the status of the particles
!   - phi0s: real array with phase and fraction of the field felt at which electron was launched
!   - Np: integer, number of seed electrons
!   - wave: real array with information about RF EM wave
!   - px: real 2D array with particle locations
!   - pv: real 2D array with particle velocities
!   - method: integer, method of activating electrons
!   - n: integer, time iteration counter
!
! OUTPUTS:
!   - pact: integer number of active particles
!   - phases: real array, phases at which particles are launched (2nd column: whether particle at that phase has been launched)
!   - pstat: integer array indicating the status of the particles
!   - phi0s: real array with phase and fraction of the field felt at which electron was launched
!   - px: real 2D array with particle locations
!   - pv: real 2D array with particle velocities
!


use precision_def
use def_types
use constants
use ecuyer_taus_rng
implicit none
! Subroutine arguments
type(geo_t), intent(in) :: geo
real(long), intent(in) :: d4x(4),wave(4)
integer, allocatable, intent(inout) :: pstat(:)
real(long), allocatable, intent(inout) :: px(:,:),pv(:,:),phases(:,:),phi0s(:,:)
integer, intent(in) :: atype, pstart, Np, method, n
integer, intent(inout) :: pact
! Unpacked variables
real(long) :: b, dt, w, beta
! Dummy variables
real(long) :: twopi, eps, pfact, phi_0, temp_px(3), temp_pv(3)
integer :: i, m, swaps, length, temp_stat


!!!!!!!!! Unpack variables !!!!!!!!!!!!!
b = geo%b
dt = d4x(4)
w = wave(1)
beta = wave(2)

!!!!!!!!!!!! For now this only supports method == 3 !!!!!!!!!!!
twopi = 2*pi

if (method == 3) then	! Activates particles at random phases
  
  if (phases(1,1) == 0.) then
    if (atype == 8 .OR. atype == 9) then
      do i = 1,Np
        read(16,*) phases(i,1)
      enddo
    else
      do i = 1,Np
       phases(i,1) = taus88() * twopi
      enddo
    endif
  endif
  
  eps = pi/(Np-1)	! To account for discreteness of space/time, launch if located within +/-eps of the required phase

  do m = 1,Np
    if (pstat(m) == 0) then
      pfact = sin(pi*px(m,2)/b);
      phi_0 = modulo(w*(n+1)*dt-beta*px(m,3),twopi)	! Calculate current phi_0 for each particle

      do i = 1,Np
	if (phases(i,2)==0 .AND. ((phi_0 >= phases(i,1)-eps) .AND. (phi_0 <= phases(i,1)+eps))) then

	  
	  pstat(m) = 1
	  pact = pact + 1
	  phi0s(m,1) = phi_0
	  phi0s(m,2) = pfact*sin(phi_0)
	  phases(i,2) = 1
	  exit

	endif
      enddo

    endif
  enddo

  ! Activate particles that were launched at one of the random phases chosen
  do m = 1,Np
    if (n>=(pstart*3) .AND. pact/=Np .AND. pstat(m)==0) then
      pfact = sin(pi*px(m,2)/b)
      phi_0 = mod(w*(n+1)*dt-beta*px(m,3),twopi)		! Calculate current phi_0 for each particle
            
      if (pfact*sin(phi_0) > 0.6) then
        phi0s(m,1) = phi_0				! Store phase at which particle was launched
        phi0s(m,2) = pfact*sin(phi_0)		! Fraction of E0 felt by particle when launched
        pstat(m) = 1
        pact = pact + 1
      endif
    endif
  enddo

  ! Organize px, pv, pstat so that active particles are at top of array (bubble sort)
  length = Np
  swaps = 1
  do while (swaps /= 0)
    swaps = 0
        
    do i = 1,length-1
      if (pstat(i) < pstat(i+1)) then
        temp_px = px(i,:)
        temp_pv = pv(i,:)
        temp_stat = pstat(i)
        px(i,:) = px(i+1,:)
        pv(i,:) = pv(i+1,:)
        pstat(i) = pstat(i+1)
        px(i+1,:) = temp_px
        pv(i+1,:) = temp_pv
        pstat(i+1) = temp_stat
        swaps = swaps + 1
      endif
    enddo
        
    length = length - 1
  enddo

endif

end subroutine particlelaunch

end module particle_launch
