module secondary_launch

contains
subroutine activate_secondary(px2,pv2,px2t,pv2t,pstat,Np,Nt,pact,full,s,m,lost)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!				ACTIVATE SECONDARY ELECTRONS
!
! This subroutine activates secondary particles and reorganizes other arrays according to the number
! of electrons produced by subroutine particlewall.
!
! INPUTS:
!   - px2: real 2D array, particle locations
!   - pv2: real 2D array, particle velocities
!   - px2t: real 2D array, particle location of secondary electrons
!   - pv2t: real 2D array, particle velocity of secondary electrons
!   - pstat: 1D integer array, particle status
!   - Np: integer, number of seed particles
!   - Nt: integer, total number of particles simulated
!   - pact: integer, number of active particles
!   - full: boolean, indicate if all Nt particles have been activated
!   - s: integer, number of secondary electrons
!   - m: integer, particle loop counter
!   - lost: integer 1D array counting particles that escape or are absorbed
!
! OUTPUTS:
!   - px2: real 2D array, particle locations
!   - pv2: real 2D array, particle velocities
!   - pstat: 1D integer array, particle status
!   - pact: integer, number of active particles
!   - full: boolean, indicate if all Nt particles have been activated
!   - lost: integer 1D array counting particles that escape or are absorbed
!
!

use precision_def
use reshape_array
implicit none
! Subroutine Arguments
real(long), allocatable, intent(inout) :: px2(:,:),pv2(:,:)
real(long), allocatable, intent(in) :: px2t(:,:),pv2t(:,:)
integer, allocatable, intent(inout) :: pstat(:)
integer, intent(in) :: Np, Nt, s, m
integer, intent(inout) :: pact, lost(2)
logical, intent(inout) :: full
! Dummy arguments
integer :: i,j,k,l,o
integer, allocatable :: s3(:)

! Activate reflected/secondary electrons
! Put first secondary back where the primary was
! Put subsequent secondaries where electrons lost or absorbed in this same time iteration were
! If not electrons were lost/absorbed then stack subsequent secondaries at the bottom of active particles
if (s > 0) then
                            
  px2(m,:) = px2t(1,:)
  pv2(m,:) = pv2t(1,:)
  pstat(m) = 1
                            
  if (s > 1) then

    o = 0	! Counter used to store indices of electrons lost/absorbed in this time iteration
    ! Store indices of electrons lost/absorbed in this time iteration in s3
    do j = 1, m-1
      if (pstat(j) == 3 .OR. pstat(j) == 4) then
	if (allocated(s3)) then
	  call add21DarrayI(s3,1)		! Add one more element to s3
	else
	  allocate(s3(1))
	endif
	o = o + 1
	s3(o) = j
      endif
    enddo

    l = 1
    ! Put new secondaries in as many s3 locations as possible
    i = 0
    do while (i < o .AND. l < s)

      l = l + 1
      i = i + 1
      if (pstat(s3(i)) == 3) then
        lost(1) = lost(1) + 1
      elseif (pstat(s3(i)) == 4) then
        lost(2) = lost(2) + 1
      endif
      pv2(s3(i),:) = pv2t(l,:)
      px2(s3(i),:) = px2t(l,:)
      pstat(s3(i)) = 1

    enddo

    ! Put remaining secondaries at the bottom of the list
    if (l < s) then

      do k = 1,s-l

	i = k
	if (pact+k > Nt) then
	  full = .TRUE.
	  exit
	endif
        l = l+1
        pv2(pact+k,:) = pv2t(l,:)
        px2(pact+k,:) = px2t(l,:)
	pstat(pact+k) = 1

      enddo

      ! Counters keeping track of secondaries stacked at the bottom
      pact = pact + i
    endif

  endif
                        
elseif (s == 0) then
  pstat(m) = 4
endif

end subroutine activate_secondary

end module secondary_launch
