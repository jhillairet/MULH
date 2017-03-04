module particle_load

contains
subroutine particleload(atype,d4x,geo,seed,px_i,fields,px,pv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!			PARTICLE LOAD
!
! Load the seed electrons onto the computational domain. Electrons are placed
! at random locations in the center of the waveguide, unless atype=3, when electrons
! are placed at the plasma-end of the waveguide to asses effect of suprious modes.
! Electrons are initialized with random velocities from a Maxwellian distribution.
! px_i is used in case other locations are programmed in the future.
! 
! This subroutine also includes the option of using pre-calculated particle locations
! and velocities in case identical results want to be reproduced.
!
! INPUTS:
!   - atype: integer type of analysis being performed
!   - d4x: real array with the discretization steps of the 4 dimensions
!   - geo: user defined type array with info regarding geometry of waveguide and PML
!   - seed: user defined type array w/ info about seed electrons
!   - px_i: integer flag indicating how to intialize electrons
!   - fields: integer, which fields are being used
!
! OUTPUTS:
!   - px: real array with particle locations
!   - pv: real array with particle velocities
!
!

use precision_def
use ecuyer_taus_rng
use constants
use def_types
implicit none
!!!!!!!!! Variable type declaration !!!!!!!!!!!
! subroutine arguments
real(long), intent(in) :: d4x(:)
integer, intent(in) :: atype, px_i, fields
real(long), intent(out) :: px(:,:),pv(:,:)
type(geo_t) :: geo
type(seed_t) :: seed
! Unpacked variables
real(long) :: a,b,lw,dx,dy,dz,vra
integer :: Np,Nt,vth,fmax,nsamplev
! local dummy variables
real(long) :: xb0,xb1,yb0,yb1,zb0,zb1,q,r,s,rad,theta,phi,vthr
integer :: i1,i2,i3,val(8),i,m,j

!!!!!!! Unpack info from input arguments !!!!!!!
a = geo%a
b = geo%b
lw = geo%lw
dx = d4x(1)
dy = d4x(2)
dz = d4x(3)
Np = seed%Np
Nt = seed%Nt
vth = seed%vth
fmax = seed%fmax
nsamplev = seed%nsamplev
vra = seed%vra

!!!!!! Initialize particle position !!!!!!!!
if (atype == 8 .OR. atype == 9) then
  do i = 1,Np
    read(14,*) (px(i,j), j=1,3)
    read(15,*) (pv(i,j), j=1,3)
  enddo
else

  if (px_i == 3) then

    if (fields == 3) then

      zb0 = lw*0.5
      zb1 = lw

    else

      zb0 = lw*0.25
      zb1 = lw*0.5

    endif
    
    xb0 = a*0.25
    xb1 = a*0.5
      
    yb0 = b*0.25
    yb1 = b*0.5
  
    ! Use Ecuyer's Tausworthe RNG. Use seconds in current time to produce i1, i2, i3 seeds
    call date_and_time(VALUES=val)
    i1 = val(8)
    i2 = val(8) + val(7)
    i3 = val(8) - val(7)
    call init_seeds(i1, i2, i3)

    do i = 1,Np
      q = taus88()
      r = taus88()
      s = taus88()
        
      px(i,1) = xb0 + q*xb1
      px(i,2) = yb0 + r*yb1
      px(i,3) = zb0 + s*zb1
    enddo

  endif

  !!!!!!!! Initialize particle velocity !!!!!!!

  if (vth /= 0) then

    vthr = sqrt(2*vth*(-eme))
    
    do i = 1,Np
      r = taus88()
      theta = 2*pi*r
      phi = pi*r
      Rad = vthr * sqrt(-log(taus88()))

      pv(i,1) = Rad * cos(theta) * sin(phi)
      pv(i,2) = Rad * sin(theta) * sin(phi)
      pv(i,3) = Rad * cos(phi)
    enddo

    do i = Np+1,Nt
      pv(i,1) = 0.
      pv(i,2) = 0.
      pv(i,3) = 0.
    enddo

  endif
endif


end subroutine particleload

end module particle_load
