module MP_detector

contains
integer function MPdetector(pact,pactts,pactps,Np,Nt,pstart,sB,n,time,pstat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!			MULTIPACTOR DETECTOR
!
! Detect whether the multipactor developed based on the number of particles over
! time. At the moment it is based on reaching a certain threshold, but ideally
! it should truly analyze the trend over time and see if increase in number of
! electrons is exponential or not.
!
! INPUTS:
!   - pact: integer, number of active particles
!   - pactts: real 1D array with time at which number of active particles have been saved
!   - pactps: real 1D array with saved number of active particles
!   - Np: integer, number of seed electrons
!   - Nt: integer, total number of simulated electrons
!   - pstart: integer, n iteration when particle dynamics began.
!   - sB: 1D real array with info about static B fields
!   - n: integer, time loop counter
!   - time: real, current simulation time
!   - pstat: 1D integer array, particle status
!
! OUTPUTS:
!   - MPdetector: integer, =1 breakdown, =2 no breakdown, =0 continue analysis.
!
!

use precision_def
implicit none
! Function arguments
integer, intent(in) :: pact, Np, Nt, pstart, n
integer, allocatable, intent(in) :: pactps(:)
real(long), allocatable, intent(in) :: pactts(:)
real(long), intent(in) :: sB(5), time
integer, allocatable, intent(in) :: pstat(:)
! Dummy variables
integer :: nop, p, changed, i

if (sB(1) == 0. .AND. sB(2) == 0. .AND. sB(3) == 0.) then

  ! Find index for time 0.5e-7s ago
  nop = size(pactts)
  p = nop
  do while (time - pactts(p) < 0.5e-7 .AND. p > 0)
    p = p-1
  enddo

  ! Determine whether number of particles has changed over that period
  do i = p,nop
    if (pact /= pactps(i)) then
      changed = 1
      exit
    else
      changed = 0
    endif
  enddo
    
  ! Determine whether MP will develop based on following criteria
  if (pact >= Nt) then
    write(*,*) '=> BREAKDOWN (1)'
    MPdetector = 1
  elseif ((n >= pstart*2) .AND. pact <= Np*0.25) then
    write(*,*) '=> NO BREAKDOWN (2)'
    MPdetector = 2
  elseif (time > 1e-6 .AND. pact >= Np*3) then
    write(*,*) '=> BREAKDOWN (3)'
    MPdetector = 1
  elseif (time >= 1e-6 .AND. pact < Np*3) then
    write(*,*) '=> NO BREAKDOWN (4)'
    MPdetector = 2
  elseif (time > 1e-7 .AND. pact < Np*3 .AND. changed == 0) then
    write(*,*) '=> NO BREAKDOWN (5)'
    MPdetector = 2
  elseif (time > 1e-7 .AND. pact >= Np*3 .AND. changed == 0) then
    write(*,*) '=> BREAKDOWN (6)'
    MPdetector = 2
  else
    MPdetector = 0
  endif
  
else
    
  ! Find index for time 0.5e-7s ago
  nop = size(pactts)
  p = nop
  do while (time - pactts(p) < 0.5e-7)
    p = p-1
  enddo
    
  ! Determine whether number of particles has changed over that period
  do i = p,nop
    if (pact /= pactps(i)) then
      changed = 1
    else
      changed = 0	
    endif
  enddo
    
  ! Determine whether MP will develop based on following criteria
  if (pact >= Nt) then
    write(*,*) '=> BREAKDOWN (1)'
    MPdetector = 1
  elseif ((n >= pstart*2) .AND. pact <= Np*0.25) then
    write(*,*) '=> NO BREAKDOWN (2)'
    MPdetector = 2
  elseif (time > 1e-6 .AND. pact >= Np*3) then
    write(*,*) '=> BREAKDOWN (3)'
    MPdetector = 1
  elseif (time >= 1e-6 .AND. pact < Np*3) then
    write(*,*) '=> NO BREAKDOWN (4)'
    MPdetector = 2
  elseif (time > 1e-7 .AND. pact < Np*3 .AND. changed == 0) then
    write(*,*) '=> NO BREAKDOWN (5)'
    MPdetector = 2
  elseif (time > 1e-7 .AND. pact >= Np*3 .AND. changed == 0) then
    write(*,*) '=> BREAKDOWN (6)'
    MPdetector = 2
  else
    MPdetector = 0
  endif
    
endif


return
end function MPdetector

end module MP_detector        

