module mod_Vaughan

contains
subroutine modVaughan(E_0,E_0p,E1,Emax,deltamax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Modify Vaughan's Formula
!
! This routine modifies Vaughan's formula so that given E_1, Emax and deltamax
! E_0 is matched to the experimental data. This is useful because previous studies
! have found E1 to be one of the most important parameters when determining the 
! multipactor threshold.
!
! INPUTS:
!   - E1: integer, first crossover energy
!   - Emax: integer, energy at which yield is max
!   - deltamax: real, maximum secondary yield
!
! OUTPUTS:
!   - E_0: integer, Vaughan E_0 parameter
!   - E_0p: integer, E_0 plus parameter
!
!


use precision_def
use def_types
implicit none
! Subroutine arguments
integer, intent(out) :: E_0, E_0p
integer, intent(in) :: E1, Emax
real(long), intent(inout) :: deltamax
! Dummy arguments
integer :: penf1, found, from_below, from_above, E_0s
real(long) :: xi, k, delta, deltas

if (E1==35 .AND. Emax==165 .AND. abs(deltamax-2.3) < 1e-6) then
  ! This is FEST Cu
  E_0 = 21
  E_0p = 20
else
    
  ! Initial values
  penf1 = E1
  E_0 = 13             ! Cutoff energy in eV

  found = 0               ! Flag indicating E_0 & E_0p were found
  from_below = 0
  from_above = 0
  do while (found == 0)

    ! Calculate yield for E1   
    xi = (real(penf1 - E_0))/(Emax - E_0)
    if (xi <= 3.6) then
      if (xi <= 1) then
        k = 0.56
      elseif (xi > 1 .OR. xi <= 3.6) then
        k = 0.25
      endif
        delta = deltamax * ( (xi * exp(1-xi))**k )
    elseif (xi > 3.6) then
      delta = deltamax * 1.125 / (xi**0.35)
    endif

    ! Compare yield to 1 to decide whether to reduce or to increase E_0
    if (delta < 1) then

      if (from_above == 0) then
        from_below = 1
        E_0 = E_0 - 1
      else
        found = 1
        E_0s = E_0 + 1
      endif


    elseif (delta > 1) then

      if (from_below == 0) then
        from_above = 1
        E_0 = E_0 + 1
      else
        found = 1
        E_0s = E_0 - 1
      endif

    endif

    ! Calculate delta for one more reduction/increment of E_0 and decide which one gives a yield closer to 1 for E1
    if (found == 1) then
      xi = real((penf1 - E_0s))/(Emax - E_0s)
      if (xi <= 3.6) then
        if (xi <= 1) then
          k = 0.56
        elseif (xi > 1 .OR. xi <= 3.6) then
          k = 0.25
        endif
        deltas = deltamax * ( (xi * exp(1-xi))**k )
      elseif (xi > 3.6) then
        deltas = deltamax * 1.125 / (xi**0.35)
      endif

      if (abs(deltas-1) < abs(delta-1)) then
        E_0 = E_0s
      endif
    endif  

  enddo

  if (E_0 <= 10) then
    E_0p = 10
  else
    E_0p = E_0-1
  endif
endif

end subroutine modVaughan

end module mod_Vaughan
