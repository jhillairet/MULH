module seem

contains
subroutine see(atype,pvf,s,pcol,pvf1,sey,SEYvsE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!			SECONDARY ELECTRON EMISSION MODEL
!
! Numerics of the Secondary Electron Emission model. seec is used to decide between a model
! that spun as a hybrid of Cheng's, Furman & Pivi's and analysing FEST3D results and the
! Furman & Pivi model.
!
! INPUTS
!   - atype: integer type of analysis being performed
!   - pcol: integer array indicating which wall electron collided with.
!   - pvf1: real 1D array with velocity of primary electron at time of impact.
!   - sey: user defined array with info about SEY model
!   - SEYvsE: array used for storing SEY at each impacting energy.
!
! OUTPUTS:
!   - pvf: real 2D array with velocity of outgoing electrons at time of ejection
!   - s: integer, number of outgoing electrons
!
!

use precision_def
use def_types
use constants
use random
use ecuyer_taus_rng
use reflect
use incomplete_gamma
use invincbeta
use reshape_array
implicit none
! Subroutine arguments
real(long), intent(out), allocatable :: pvf(:,:)
integer, intent(out) :: s
real(long), intent(in) :: pvf1(3)
type(sey_t), intent(in) :: sey
integer, intent(in) :: atype, pcol(4)
real(long), allocatable, intent(inout) :: SEYvsE(:,:,:)
! Unpacked variables
real(long) :: Eom, deltamax, p_n, delta_b, E_0, E_0p, a_lara, ks, kse, z_lara
integer :: Emax, ReRr, seec
! Dummy variables
real(long) :: pvmf1,penf1,theta,ke,delta,Rr,Ee1,Ee2,Re,C1,a0,a1,a2,a3,e0,f,Pe,Pr, &
		r,ans,qans,y,sint,alpha,Es,lnbeta,inbeta,xi,rt,delta0,P_1e_inf,Phat_1e,Ehat_e,W,p, &
		sigma_e,e1,e2,P_1r_inf,E_r,q,r1,r2,deltahat_ts,Ehat_ts,S_fp,t1,t2,t3,t4, &
		delta_e,delta_r,Ehat,deltahat,delta_ts,delta_tsp,rs,pu,pl,a_e,a_r,u,x_0
real(long), allocatable :: yn(:),penf(:),vs(:),phis(:),thetas(:),P_n_ts(:)
integer :: g, gi(1), Eb, ierr, i, j, ifault, k, h, M
logical, save :: first=.TRUE.



Eom = sey%Eom      ! Average energy of Maxwellian distribution of secondary electrons emitted in joules
pvmf1 = absolute(pvf1(1:3))             ! Magnitude of impact velocity
penf1 = 0.5*me*pvmf1*pvmf1/(-e)     ! Particle Energy in eV
seec = sey%seec      ! SEE model, =1 Vaughan, =2 Furman & Pivi
p_n = sey%p_n

!write(*,*) 'penf1 = ',penf1

gi = maxloc(pcol(1:4))
g = gi(1)
! Angle of incidence
if (g<3) then
  theta = acos( pvf1(2-mod(g,2)) / pvmf1)
elseif (g>=3) then
  theta = acos( - pvf1(2-mod(g,2)) / pvmf1)
endif

if (seec == 1) then	!************* Vaughan model **************!
    
  Emax = sey%Emax     ! Emax(delta=max,theta=0) in eV
  deltamax = sey%deltamax      ! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
  ReRr = sey%ReRr
  delta_b = sey%delta_b
  E_0 = sey%E_0
  E_0p = sey%E_0p
       
  if (penf1 <= E_0p) then

    delta = delta_b

  elseif (penf1 > E_0p .AND. penf1 < E_0) then
        
    delta = -delta_b*(penf1 - E_0)
        
  else
        
    ! Unpack user inputs regarding particle-wall interaction and SEE
    kse = sey%kse           ! Roughness factor for energy, =0 for rough, =1 for dull, =2 for smooth and anything in between
    ks = sey%ks            ! Roughness factor for angle, =0 for rough, =1 for dull, =2 for smooth and anything in between

    ! Determine secondary electron yields (Vaughan's theory (1993))
    Emax = Emax * (1 + kse*theta*theta/(2*pi))
    deltamax = deltamax * (1 + ks*theta*theta/(2*pi))

    xi = real((penf1 - E_0))/(Emax - E_0)

    if (xi <= 3.6) then

      if (xi <= 1) then
        ke = 0.56
      elseif (xi > 1 .OR. xi <= 3.6) then
        ke = 0.25
      endif

        delta = deltamax * ( (xi * exp(1-xi))**ke )

    elseif (xi > 3.6) then

        delta = deltamax * 1.125 / (xi**0.35)

    endif

  endif

  !if (allocated(SEYvsE)) then
   ! call add23Darray(SEYvsE,(/0,0,1/))
    !k = size(SEYvsE,3)
    !SEYvsE(4,1,k) = penf1
    !SEYvsE(4,2,k) = delta
  !else
   ! allocate (SEYvsE(4,2,1))
   ! SEYvsE(4,1,1) = penf1
    !SEYvsE(4,2,1) = delta
  !endif

  if (delta > 10) then
    write(*,*) 'high delta'
    call exit
  endif

  !if (first) then
  !  s = random_Poisson(real(delta,4),.TRUE.)        ! Number of secondary electrons produced
  !  first = .FALSE.
  !else
  !  s = random_Poisson(real(delta,4),.FALSE.)
  !endif

  if (isNaN(real(delta,8))) then
    write(*,*) 'theta = ',theta
    write(*,*) 'penf1 = ',penf1
    write(*,*) 'E_0 = ',E_0
    write(*,*) 'E_0p = ',E_0p
    write(*,*) 'delta = ',delta
    write(*,*) 'E_0p = ',E_0p
    write(*,*) 'E_0 = ',E_0
  endif

  if (atype == 9) then
    read(17,'(F17.15,1X,I1)') delta0, s
    if (delta0-delta>1.5e-3) then
      write(*,*) 'discrepancy in delta'
      write(*,*) 'penf1 = ',penf1
      write(*,*) 'theta = ',theta
      write(*,*) 'delta0 = ',delta0
      write(*,*) 'delta = ',delta
      write(*,*) 'diff = ',delta0-delta
      call exit
    endif
  else
    s = poisdev(real(delta,8),atype)
  endif
    
  if (s == 1) then

    allocate (pvf(1,3))	! Only one electron comes off the wall
    
    if (ReRr == 1) then
      ! Determine Re and Rr from de Lara's (2006) empiricial fit
      a_lara = sey%a_lara
      z_lara = sey%z_lara
      Eb = 300 + 175*z_lara
      Rr = a_lara * (1 - (3e-5)*penf1) * (penf1**0.56) * exp(-((penf1/Eb)**0.70))

      Ee1 = 50./sqrt(z_lara)
      Ee2 = 0.25 * z_lara * z_lara
      Re = 0.93/(1+penf1/Ee1) + 0.07/(1+penf1/Ee2)

      ! Take into account angle dependence of Re and Rr
      C1 = 0.89 * Rr / (Re + Rr)
      Rr = (Rr**cos(theta)) * (C1**(1-cos(theta)))
      Re = (Re**cos(theta)) * ((0.89-C1)**(1-cos(theta)))
        
    elseif (ReRr == 2) then	! Empirical fit from LHC
      a0 = 20.69989
      a1 = -7.07605
      a2 = 0.483547
      a3 = 0
      e0 = 56.914686
      f = exp(a0 + a1*log(penf1+e0) + a2*((log(penf1+e0))**2) + a3*((log(penf1+e0))**3))
        
      ! Based on publications saying that 7% and 3% of secondaries are elastically and inelastically reflected respectively
      Re = f*0.3
      Rr = f-Re
    endif

!open(unit=63,file='ReRr.txt',status='unknown',position='append')
!write(63,*) Re,Rr
!close(unit=63) 

    ! Determine probability of events
    Pe = Re            ! Prob that electron is elastically reflected
    Pr = Rr            ! Prob that electron is inelastically reflected
    Pe = 0.03
    Pr = 0.07
    
    ! Now generate a random number uniformly distributed between 0 and 1
    if (atype == 8) then
      read(21,*) r
    elseif (atype == 9) then
      read(18,*) r
    else
      r = taus88()
    endif
    
    if (r < Pe) then

      ! Reflect electron elastically
      pvf = reflect_electron(pcol,penf1*(-e),pvf1,1)

      !if (allocated(SEYvsE)) then
	!call add23Darray(SEYvsE,(/0,0,1/))
	!k = size(SEYvsE,3)
        !SEYvsE(1,1,k) = penf1
        !SEYvsE(1,2,k) = delta
      !else
       ! allocate (SEYvsE(4,2,1))
        !SEYvsE(1,1,1) = penf1
        !SEYvsE(1,2,1) = delta
      !endif

    elseif (r >= Pe .AND. r < (Pe+Pr)) then

      ! Electron is reflected inelastically
      !pvf = reflect_electron(pcol,penf1*(-e),pvf1,2,atype)

      allocate (thetas(1))
      allocate (phis(1))
      allocate (penf(1))
      allocate (vs(1))

      penf = taus88() * penf1
      vs = sqrt(2*penf*(-eme))
                        
      ! Determine which is the normal component of the outgoing velocity
      if (pcol(1)==1 .OR. pcol(3)==1) then
        i = 1              ! Normal component
        j = 2              ! Parallel component(not z)
      elseif (pcol(2)==1 .OR. pcol(4)==1) then
        i = 2              ! Normal component
        j = 1              ! Parallel component (not z)
      endif
                
      rt = taus88()
      phis = rt * 2 * pi        ! Calculate emission azimuthal angle of secondaries

      rt = taus88()
      thetas = asin(2*rt - 1)   ! Calculate emission angle of secondaries with respect to the normal
            
      pvf(1,i) = vs(1) * cos(thetas(1))   ! Normal component of secondary velocity

      pvf(1,j) = sin(phis(1)) * vs(1) * sin(thetas(1))
      pvf(1,3) = cos(phis(1)) * vs(1) * sin(thetas(1))

      if (pcol(1) == 1 .OR. pcol(2) == 1) then
        pvf(1,i) = - abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
        pvf(1,i) = abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      endif
         
      deallocate (thetas)
      deallocate (phis)
      deallocate (penf)
      deallocate (vs)

      !if (allocated(SEYvsE)) then
	!call add23Darray(SEYvsE,(/0,0,1/))
	!k = size(SEYvsE,3)
        !SEYvsE(2,1,k) = penf1
        !SEYvsE(2,2,k) = delta
      !else
       ! allocate (SEYvsE(4,2,1))
        !SEYvsE(2,1,1) = penf1
        !SEYvsE(2,2,1) = delta
      !endif

    else
         
      !if (allocated(SEYvsE)) then
	!call add23Darray(SEYvsE,(/0,0,1/))
	!k = size(SEYvsE,3)
        !SEYvsE(3,1,k) = penf1
        !SEYvsE(3,2,k) = delta
      !else
       ! allocate (SEYvsE(4,2,1))
        !SEYvsE(3,1,1) = penf1
        !SEYvsE(3,2,1) = delta
      !endif

      allocate (yn(1))
      allocate (thetas(1))
      allocate (phis(1))
      allocate (penf(1))
      allocate (vs(1))   

      ! Initialize variables
      pvf = 0.       ! velocity of 1 secondary electron

      ! Scalar parameters needed to obtain velocity components
      call gratio(real(p_n*s,8), real(penf1/Eom,8), ans, qans, 0)
      if (atype == 8) then
        read(21,*) rt
        ans = ans*rt
      elseif (atype == 9) then
        read(20,*) rt
        ans = ans*rt
      else
        ans = ans*taus88()
      endif
      call gaminv(real(p_n*s,8), yn(1), real(0,8), ans, 1-ans, ierr)

      if (ierr < 0) then
	write(*,*) 'error in calculating inverse incomplete gamma function. Terminating.'
	call exit
      endif

      yn = sqrt( yn ) ! magnitude of normalized velocity of secondary electrons

      penf = Eom * yn * yn
      vs = sqrt(2*penf*(-eme))
                        
      ! Determine which is the normal component of the outgoing velocity
      if (pcol(1)==1 .OR. pcol(3)==1) then
        i = 1              ! Normal component
        j = 2              ! Parallel component(not z)
      elseif (pcol(2)==1 .OR. pcol(4)==1) then
        i = 2              ! Normal component
        j = 1              ! Parallel component (not z)
      endif
                
      if (atype == 8) then
        read(21,*) rt
      elseif (atype == 9) then
        read(20,*) rt
      else
        rt = taus88()
      endif
      phis = rt * 2 * pi        ! Calculate emission azimuthal angle of secondaries

      if (atype == 8) then
        read(21,*) rt
      elseif (atype == 9) then
        read(20,*) rt
      else
        rt = taus88()
      endif
      thetas = asin(2*rt - 1)   ! Calculate emission angle of secondaries with respect to the normal
            
      pvf(1,i) = vs(1) * cos(thetas(1))   ! Normal component of secondary velocity

      pvf(1,j) = sin(phis(1)) * vs(1) * sin(thetas(1))
      pvf(1,3) = cos(phis(1)) * vs(1) * sin(thetas(1))

      if (pcol(1) == 1 .OR. pcol(2) == 1) then
        pvf(1,i) = - abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
        pvf(1,i) = abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      endif
         
      deallocate (yn)
      deallocate (thetas)
      deallocate (phis)
      deallocate (penf)
      deallocate (vs)
   
    endif
  elseif (s > 1) then

      !if (allocated(SEYvsE)) then
	!call add23Darray(SEYvsE,(/0,0,1/))
	!k = size(SEYvsE,3)
        !SEYvsE(3,1,k) = penf1
        !SEYvsE(3,2,k) = delta
      !else
       ! allocate(SEYvsE(4,2,1))
        !SEYvsE(3,1,1) = penf1
        !SEYvsE(3,2,1) = delta
      !endif

    allocate (pvf(s,3))
    allocate (yn(s))
    allocate (thetas(s))
    allocate (phis(s))
    allocate (penf(s))
    allocate (vs(s))
        
    ! Initialize variables
    pvf = 0.       ! m x 3 array of 3D velocities of m secondary electrons
    yn = 0.        ! magnitude of velocity of secondary electrons
    thetas = 0.     ! Emission angle of secondaries with respect to the normal
    phis = 0.       ! Azimuthal emission angle of secondaries
    penf = 0.
    vs = 0.

    call gratio(real(p_n*s,8), real(penf1/Eom,8), ans, qans, 0)
    if (atype == 8) then
      read(21,*) rt
    elseif (atype == 9) then
      read(20,*) rt
    else
      rt = taus88()
    endif
    ans = ans*rt
    call gaminv(real(p_n*s,8), y, real(0,8), ans, 1-ans, ierr)

    Y = sqrt( Y )    ! v^tilde

    !! Scalar parameters needed to obtain velocity components
    sint = 1.
    do k = 1,s-1
      lnbeta = betaln(real(p_n*(s-k),4),real(p_n,4))

      if (atype == 8) then
        read(21,*) rt
      elseif (atype == 9) then
        read(20,*) rt
      else
        rt = taus88()
      endif
      inbeta = betain( real(rt,8), real(p_n*(s-k),8), real(p_n,8), lnbeta, ifault )
      if (ifault /= 0) then
	write(*,*) 'error calculating incomplete beta function. Terminating.'
	call exit
      endif

      alpha = xinbta( real(p_n*(s-k),8), real(p_n,8), lnbeta, inbeta, ifault )
      if (ifault /= 0) then
	write(*,*) 'error calculating inverse incomplete beta function. Terminating.'
	call exit
      endif

      alpha = asin( sqrt( alpha ) )	! Alpha angles used to calculate the magnitude of velocity

      yn(k) = Y * sint * cos(alpha)	! magnitude of outgoing velocities

      sint = sint * sin(alpha)		! Spherical coordinates factor                
    enddo
    yn(s) = Y * sint
                        
    ! Determine which is the normal component of the outgoing velocity
    if (pcol(1)==1 .OR. pcol(3)==1) then
      i = 1              ! Normal component
      j = 2              ! Parallel component(not z)
    elseif (pcol(2)==1 .OR. pcol(4)==1) then
      i = 2              ! Normal component
      j = 1              ! Parallel component (not z)
    endif

    do k = 1,s
      penf(k) = Eom * (yn(k))**2
      vs(k) = sqrt(2*penf(k)*(-eme))

      if (atype == 8) then
        read(21,*) rt
      elseif (atype == 9) then
        read(20,*) rt
      else
        rt = taus88()
      endif
      phis(k) = rt * 2 * pi        ! Calculate emission azimuthal angle of secondaries

      if (atype == 8) then
        read(21,*) rt
      elseif (atype == 9) then
        read(20,*) rt
      else
        rt = taus88()
      endif
      thetas(k) = asin(2*rt - 1)   ! Calculate emission angle of secondaries with respect to the normal
            
      pvf(k,i) = vs(k) * cos(thetas(k))   ! Normal component of secondary velocity
      pvf(k,j) = sin(phis(k)) * vs(k) * sin(thetas(k))
      pvf(k,3) = cos(phis(k)) * vs(k) * sin(thetas(k))
    enddo

    if (pcol(1) == 1 .OR. pcol(2) == 1) then
      do k = 1,s
	pvf(k,i) = - abs(pvf(k,i))   ! Normal velocity pointing into the inside of the waveguide
      enddo
    elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
      do k = 1,s
        pvf(k,i) = abs(pvf(k,i))   ! Normal velocity pointing into the inside of the waveguide
      enddo
    endif

    deallocate (yn)
    deallocate (thetas)
    deallocate (phis)
    deallocate (penf)
    deallocate (vs)
            
  elseif (s == 0) then
    allocate (pvf(1,3))
    pvf = 0.
  endif

elseif (seec == 2) then	

  !************** Furman & Pivi SEE model *****************!

  ! Fitting parameters
  P_1e_inf = 0.022
  Phat_1e = 0.496
  Ehat_e = 0
  W = 60.86
  p = 1
  sigma_e = 2
  e1 = 0.26
  e2 = 2
  P_1r_inf = 0.22
  E_r = 0.041
  r = 0.104
  q = 0.5
  r1 = 0.26
  r2 = 2
  deltahat_ts = 2.3-P_1e_inf-P_1r_inf
  Ehat_ts = 165!276.8
  S_fp = 1.54
  t1 = 0.66
  t2 = 0.8
  t3 = 0.7
  t4 = 1
    
  M = 10

  delta_e = P_1e_inf + (Phat_1e - P_1e_inf)*exp(- ((abs(penf1-Ehat_e)/W)**p)/p)
  delta_e = delta_e * (1 + e1*(1 - (cos(theta))**e2))
    
  delta_r = P_1r_inf * (1 - exp(- (penf1/E_r)**r))
  delta_r = delta_r * (1 + r1*(1 - (cos(theta))**r2))
    
  deltahat = deltahat_ts * (1 + t1*(1 - (cos(theta))**t2))
  Ehat = Ehat_ts * (1 + t3*(1 - (cos(theta))**t4))
    
  delta_ts = deltahat * ( S_fp * (penf1/Ehat)) / (S_fp - 1 + (penf1/Ehat)**S_fp)
    
  delta_tsp = delta_ts/(1. - delta_e - delta_r)
    
  allocate(P_n_ts(M))
  P_n_ts = 0.
  rs = taus88()
  pu = 0.
  pl = 0.

  do h = 1,M
    P_n_ts(h) = (delta_tsp**(h-1)) * exp(-delta_tsp) / factorial(h-1)
        
    if (h == 2) then
      P_n_ts(h) = (1 - delta_e - delta_r)*P_n_ts(h) + delta_e + delta_r
    else
      P_n_ts(h) = (1 - delta_e - delta_r)*P_n_ts(h)
    endif
        
    if (h > 1) pl = pl + P_n_ts(h-1)
    pu = pu + P_n_ts(h)
        
    if (rs>=pl .AND. rs<pu) then
      s = h - 1
    elseif (rs >= pu) then
      s = h
    endif
  enddo

  if (s == 1) then

    allocate (thetas(1))
    allocate (phis(1))
    allocate (penf(1))
    allocate (vs(1)) 
        
    a_e = delta_e/P_n_ts(2)
    a_r = delta_r/P_n_ts(2)
        
    u = taus88()
        
    if (u < a_e) then
      penf = -1.
      do while (penf(1) < 0)
        penf = penf1 - sigma_e * abs(random_normal())
      enddo
            
    elseif (u >= a_e .AND. u < a_e+a_r) then
            
      penf = penf1 * (taus88()**(1/(1+q)))
            
    else
      
      allocate (yn(1))
      ! Scalar parameters needed to obtain velocity components
      call gratio(real(p_n*s,8), real(penf1/Eom,8), ans, qans, 0)
      ans = ans*taus88()

      call gaminv(real(p_n*s,8), yn(1), real(0,8), ans, 1-ans, ierr)

      if (ierr < 0) then
	write(*,*) 'error in calculating inverse incomplete gamma function. Terminating.'
	call exit
      endif

      penf = Eom * yn
            
    endif
        
    vs = sqrt(2*penf*(-e)/me)
        
  elseif (s > 1) then
        
    allocate (yn(s))
    allocate (thetas(s))
    allocate (phis(s))
    allocate (penf(s))
    allocate (vs(s))
        
    x_0 = penf1/Eom

    call gratio(real(p_n*s,8), real(x_0,8), ans, qans, 0)
    rt = taus88()
    ans = ans*rt
    call gaminv(real(p_n*s,8), y, real(0,8), ans, 1-ans, ierr)

    Y = sqrt( Y )    ! Y = sqrt(P_0)
        
    SINt = 1
    do k = 1,s-1
      lnbeta = betaln(real(p_n*(s-k),4),real(p_n,4))
      rt = taus88()
      inbeta = betain( real(rt,8), real(p_n*(s-k),8), real(p_n,8), lnbeta, ifault )
      if (ifault /= 0) then
	write(*,*) 'error calculating incomplete beta function. Terminating.'
	call exit
      endif

      alpha = xinbta( real(p_n*(s-k),8), real(p_n,8), lnbeta, inbeta, ifault )
      if (ifault /= 0) then
	write(*,*) 'error calculating inverse incomplete beta function. Terminating.'
	call exit
      endif

      alpha = asin( sqrt( alpha ) )	! Alpha angles used to calculate the magnitude of velocity

      yn(k) = Y * sint * cos(alpha)	! magnitude of outgoing velocities

      sint = sint * sin(alpha)		! Spherical coordinates factor                
    enddo
    yn(s) = Y * sint
        
    do k = 1,s
      penf(k) = Eom * (yn(k))**2
      vs(k) = sqrt(2*penf(k)*(-e)/me)
    enddo
           
  endif

  if (s > 0) then
        
    allocate (pvf(s,3))       ! m x 3 array of 3D velocities of m secondary electrons

    ! Determine which is the normal component of the outgoing velocity
    if (pcol(1)==1 .OR. pcol(3)==1) then
      i = 1              ! Normal component
      j = 2              ! Parallel component(not z)
    elseif (pcol(2)==1 .OR. pcol(4)==1) then
      i = 2              ! Normal component
      j = 1              ! Parallel component (not z)
    endif

    do k = 1,s
      phis(k) = taus88() * 2 * pi            ! Calculate azimuthal emission angle of secondaries
      thetas(k) = asin(2*taus88() - 1)       ! Calculate emission angle of secondaries with respect to the normal

      pvf(k,i) = vs(k) * cos(thetas(k))   ! Normal component of secondary velocity

      pvf(k,j) = sin(phis(k)) * vs(k) * sin(thetas(k))
      pvf(k,3) = cos(phis(k)) * vs(k) * sin(thetas(k))
    enddo
        
    if (pcol(1) == 1 .OR. pcol(2) == 1) then
      pvf(:,i) = - abs(pvf(:,i))   ! Normal velocity pointing into the inside of the waveguide
    elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
      pvf(:,i) = abs(pvf(:,i))   ! Normal velocity pointing into the inside of the waveguide
    endif
    
  else 
    allocate (pvf(1,3))     
    pvf = 0.   
  endif

elseif (seec == 3) then

  ! Cheng model with de Lara's approach to calculating Re and Rr

  ! Determine Re and Rr from de Lara's paper
  a_lara = sey%a_lara
  z_lara = sey%z_lara
  Eb = 300 + 175*z_lara
  Rr = a_lara * (1 - (3e-5)*penf1) * (penf1**0.56) * exp(-((penf1/Eb)**0.70))

  Ee1 = 50./sqrt(z_lara)
  Ee2 = 0.25 * z_lara * z_lara
  Re = 0.93/(1+penf1/Ee1) + 0.07/(1+penf1/Ee2)

  ! Take into account angle dependence of Re and Rr
  C1 = 0.89 * Rr / (Re + Rr)
  Rr = (Rr**cos(theta)) * (C1**(1-cos(theta)))
  Re = (Re**cos(theta)) * ((0.89-C1)**(1-cos(theta)))

  ! Determine probability of events
  Pe = Re            ! Prob that electron is elastically reflected
  Pr = Rr            ! Prob that electron is inelastically reflected

  r = taus88()
  if (r < Pe) then

    s = 1
    allocate (pvf(1,3))	! Only one electron comes off the wall
    pvf = 0.       ! velocity of 1 secondary electron
    ! Reflect electron elastically
    pvf = reflect_electron(pcol,penf1*(-e),pvf1,1)

  elseif (r >= Pe .AND. r < (Pe+Pr)) then

    s = 1
    allocate (pvf(1,3))	! Only one electron comes off the wall
    pvf = 0.       ! velocity of 1 secondary electron
    ! Electron is reflected inelastically
    pvf = reflect_electron(pcol,penf1*(-e),pvf1,2,atype)

  else

    Emax = sey%Emax     ! Emax(delta=max,theta=0) in eV
    deltamax = sey%deltamax      ! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
    delta_b = sey%delta_b
    E_0 = sey%E_0
    E_0p = sey%E_0p
       
    if (penf1 <= E_0p) then

      delta = delta_b

    elseif (penf1 > E_0p .AND. penf1 < E_0) then
     
      delta = -delta_b*(penf1 - E_0)
        
    else

        
      ! Unpack user inputs regarding particle-wall interaction and SEE
      kse = sey%kse           ! Roughness factor for energy, =0 for rough, =1 for dull, =2 for smooth and anything in between
      ks = sey%ks            ! Roughness factor for angle, =0 for rough, =1 for dull, =2 for smooth and anything in between

      ! Determine secondary electron yields (Vaughan's theory (1993))
      Emax = Emax * (1 + kse*theta*theta/(2*pi))
      deltamax = deltamax * (1 + ks*theta*theta/(2*pi))

      xi = real((penf1 - E_0))/(Emax - E_0)

      if (xi <= 3.6) then

        if (xi <= 1) then
          ke = 0.56
        elseif (xi > 1 .OR. xi <= 3.6) then
          ke = 0.25
        endif

        delta = deltamax * ( (xi * exp(1-xi))**ke )

      elseif (xi > 3.6) then

        delta = deltamax * 1.125 / (xi**0.35)

      endif

    endif

    if (delta > 10) then
      write(*,*) 'high delta'
      call exit
    endif

    if (isNaN(real(delta,8))) then
      write(*,*) 'theta = ',theta
      write(*,*) 'penf1 = ',penf1
      write(*,*) 'E_0 = ',E_0
      write(*,*) 'E_0p = ',E_0p
      write(*,*) 'delta = ',delta
      write(*,*) 'E_0p = ',E_0p
      write(*,*) 'E_0 = ',E_0
    endif

    s = poisdev(real(delta,8),atype)

    if (s == 1) then

      allocate (pvf(1,3))	! Only one electron comes off the wall
      allocate (yn(1))
      allocate (thetas(1))
      allocate (phis(1))
      allocate (penf(1))
      allocate (vs(1))   

      ! Initialize variables
      pvf = 0.       ! velocity of 1 secondary electron

      ! Scalar parameters needed to obtain velocity components
      call gratio(real(p_n*s,8), real(penf1/Eom,8), ans, qans, 0)
      ans = ans*taus88()
      call gaminv(real(p_n*s,8), yn(1), real(0,8), ans, 1-ans, ierr)

      if (ierr < 0) then
	write(*,*) 'error in calculating inverse incomplete gamma function. Terminating.'
	call exit
      endif

      yn = sqrt( yn ) ! magnitude of velocity of secondary electrons
                        
      ! Determine which is the normal component of the outgoing velocity
      if (pcol(1)==1 .OR. pcol(3)==1) then
        i = 1              ! Normal component
        j = 2              ! Parallel component(not z)
      elseif (pcol(2)==1 .OR. pcol(4)==1) then
        i = 2              ! Normal component
        j = 1              ! Parallel component (not z)
      endif

      penf = Eom * yn * yn
      vs = sqrt(2*penf*(-eme))
      rt = taus88()
      phis = rt * 2 * pi        ! Calculate emission azimuthal angle of secondaries

      rt = taus88()
      thetas = asin(2*rt - 1)   ! Calculate emission angle of secondaries with respect to the normal
            
      pvf(1,i) = vs(1) * cos(thetas(1))   ! Normal component of secondary velocity

      pvf(1,j) = sin(phis(1)) * vs(1) * sin(thetas(1))
      pvf(1,3) = cos(phis(1)) * vs(1) * sin(thetas(1))

      if (pcol(1) == 1 .OR. pcol(2) == 1) then
        pvf(1,i) = - abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
        pvf(1,i) = abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      endif
         
      deallocate (yn)
      deallocate (thetas)
      deallocate (phis)
      deallocate (penf)
      deallocate (vs)

    elseif (s > 1) then

      allocate (pvf(s,3))
      allocate (yn(s))
      allocate (thetas(s))
      allocate (phis(s))
      allocate (penf(s))
      allocate (vs(s))
        
      ! Initialize variables
      pvf = 0.       ! m x 3 array of 3D velocities of m secondary electrons
      yn = 0.        ! magnitude of velocity of secondary electrons
      thetas = 0.     ! Emission angle of secondaries with respect to the normal
      phis = 0.       ! Azimuthal emission angle of secondaries
      penf = 0.
      vs = 0.

      call gratio(real(p_n*s,8), real(penf1/Eom,8), ans, qans, 0)
      rt = taus88()
      ans = ans*rt
      call gaminv(real(p_n*s,8), y, real(0,8), ans, 1-ans, ierr)

      Y = sqrt( Y )    ! v^tilde

      !! Scalar parameters needed to obtain velocity components
      sint = 1.
      do k = 1,s-1
        lnbeta = betaln(real(p_n*(s-k),4),real(p_n,4))

        rt = taus88()
        inbeta = betain( real(rt,8), real(p_n*(s-k),8), real(p_n,8), lnbeta, ifault )
        if (ifault /= 0) then
	  write(*,*) 'error calculating incomplete beta function. Terminating.'
	  call exit
        endif

        alpha = xinbta( real(p_n*(s-k),8), real(p_n,8), lnbeta, inbeta, ifault )
        if (ifault /= 0) then
	  write(*,*) 'error calculating inverse incomplete beta function. Terminating.'
	  call exit
        endif

        alpha = asin( sqrt( alpha ) )	! Alpha angles used to calculate the magnitude of velocity

        yn(k) = Y * sint * cos(alpha)	! magnitude of outgoing velocities

        sint = sint * sin(alpha)		! Spherical coordinates factor                
      enddo
      yn(s) = Y * sint
                        
      ! Determine which is the normal component of the outgoing velocity
      if (pcol(1)==1 .OR. pcol(3)==1) then
        i = 1              ! Normal component
        j = 2              ! Parallel component(not z)
      elseif (pcol(2)==1 .OR. pcol(4)==1) then
        i = 2              ! Normal component
        j = 1              ! Parallel component (not z)
      endif

      do k = 1,s
        penf(k) = Eom * (yn(k))**2
        vs(k) = sqrt(2*penf(k)*(-eme))

        rt = taus88()
        phis(k) = rt * 2 * pi        ! Calculate emission azimuthal angle of secondaries

        rt = taus88()
        thetas(k) = asin(2*rt - 1)   ! Calculate emission angle of secondaries with respect to the normal
              
        pvf(k,i) = vs(k) * cos(thetas(k))   ! Normal component of secondary velocity
        pvf(k,j) = sin(phis(k)) * vs(k) * sin(thetas(k))
        pvf(k,3) = cos(phis(k)) * vs(k) * sin(thetas(k))
      enddo

      if (pcol(1) == 1 .OR. pcol(2) == 1) then
        do k = 1,s
	  pvf(k,i) = - abs(pvf(k,i))   ! Normal velocity pointing into the inside of the waveguide
        enddo
      elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
        do k = 1,s
          pvf(k,i) = abs(pvf(k,i))   ! Normal velocity pointing into the inside of the waveguide
        enddo
      endif

      deallocate (yn)
      deallocate (thetas)
      deallocate (phis)
      deallocate (penf)
      deallocate (vs)

    elseif (s == 0) then
      allocate (pvf(1,3))
      pvf = 0.
    endif

  endif

elseif (seec == 4) then	! New model as of Jul12 2013. FEST3Dish?

  Emax = sey%Emax     ! Emax(delta=max,theta=0) in eV
  deltamax = sey%deltamax      ! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
  E_0 = sey%E_0
  E_0p = sey%E_0p

  if (penf1 <= E_0*4) then	! Reflect electron

      ! Reflect electron inelastically
!      pvf = reflect_electron(pcol,penf1*(-e),pvf1,2)

      allocate (thetas(1))
      allocate (phis(1))
      allocate (penf(1))
      allocate (vs(1))

      penf = taus88() * penf1
      vs = sqrt(2*penf*(-eme))
                        
      ! Determine which is the normal component of the outgoing velocity
      if (pcol(1)==1 .OR. pcol(3)==1) then
        i = 1              ! Normal component
        j = 2              ! Parallel component(not z)
      elseif (pcol(2)==1 .OR. pcol(4)==1) then
        i = 2              ! Normal component
        j = 1              ! Parallel component (not z)
      endif
                
      rt = taus88()
      phis = rt * 2 * pi        ! Calculate emission azimuthal angle of secondaries

      rt = taus88()
      thetas = asin(2*rt - 1)   ! Calculate emission angle of secondaries with respect to the normal
            
      pvf(1,i) = vs(1) * cos(thetas(1))   ! Normal component of secondary velocity

      pvf(1,j) = sin(phis(1)) * vs(1) * sin(thetas(1))
      pvf(1,3) = cos(phis(1)) * vs(1) * sin(thetas(1))

      if (pcol(1) == 1 .OR. pcol(2) == 1) then
        pvf(1,i) = - abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
        pvf(1,i) = abs(pvf(1,i))   ! Normal velocity pointing into the inside of the waveguide
      endif
         
      deallocate (thetas)
      deallocate (phis)
      deallocate (penf)
      deallocate (vs)
        
  else		! True secondary(ies)
        
    ! First calculate yield according to modified Vaughan

    ! Unpack user inputs regarding particle-wall interaction and SEE
    kse = sey%kse           ! Roughness factor for energy, =0 for rough, =1 for dull, =2 for smooth and anything in between
    ks = sey%ks            ! Roughness factor for angle, =0 for rough, =1 for dull, =2 for smooth and anything in between

    ! Determine secondary electron yields (Vaughan's theory (1993))
    Emax = Emax * (1 + kse*theta*theta/(2*pi))
    deltamax = deltamax * (1 + ks*theta*theta/(2*pi))

    xi = real((penf1 - E_0))/(Emax - E_0)

    if (xi <= 3.6) then

      if (xi <= 1) then
        ke = 0.56
      elseif (xi > 1 .OR. xi <= 3.6) then
        ke = 0.25
      endif

        delta = deltamax * ( (xi * exp(1-xi))**ke )

    elseif (xi > 3.6) then

        delta = deltamax * 1.125 / (xi**0.35)

    endif

    if (delta > 10) then
      write(*,*) 'high delta'
      call exit
    endif

    if (isNaN(real(delta,8))) then
      write(*,*) 'theta = ',theta
      write(*,*) 'penf1 = ',penf1
      write(*,*) 'E_0 = ',E_0
      write(*,*) 'E_0p = ',E_0p
      write(*,*) 'delta = ',delta
      write(*,*) 'E_0p = ',E_0p
      write(*,*) 'E_0 = ',E_0
    endif

    s = poisdev(real(delta,8),atype)

    if (s > 0) then

      allocate (pvf(s,3))
      allocate (yn(s))
      allocate (thetas(s))
      allocate (phis(s))
      allocate (penf(s))
      allocate (vs(s))
        
      ! Initialize variables
      pvf = 0.       ! m x 3 array of 3D velocities of m secondary electrons
      yn = 0.        ! magnitude of velocity of secondary electrons
      thetas = 0.     ! Emission angle of secondaries with respect to the normal
      phis = 0.       ! Azimuthal emission angle of secondaries
      penf = 0.
      vs = 0.

      ! Calculate cummulative energy of secondary electrons from a Gaussian distribution
!      Y = -1.
!      do while (Y < 0.)
!        Y = penf1 - Eom * abs(random_normal())
!      enddo
!      Y = sqrt( Y/Eom )    ! v^tilde

      call gratio(real(p_n*s,8), real(penf1/Eom,8), ans, qans, 0)
      rt = taus88()
      ans = ans*rt
      call gaminv(real(p_n*s,8), y, real(0,8), ans, 1-ans, ierr)

      Y = sqrt( Y )    ! v^tilde

      !! Scalar parameters needed to obtain velocity components
      sint = 1.
      do k = 1,s-1
        lnbeta = betaln(real(p_n*(s-k),4),real(p_n,4))

        rt = taus88()
        inbeta = betain( real(rt,8), real(p_n*(s-k),8), real(p_n,8), lnbeta, ifault )
        if (ifault /= 0) then
	  write(*,*) 'error calculating incomplete beta function. Terminating.'
	  call exit
        endif

        alpha = xinbta( real(p_n*(s-k),8), real(p_n,8), lnbeta, inbeta, ifault )
        if (ifault /= 0) then
	  write(*,*) 'error calculating inverse incomplete beta function. Terminating.'
	  call exit
        endif

        alpha = asin( sqrt( alpha ) )	! Alpha angles used to calculate the magnitude of velocity

        yn(k) = Y * sint * cos(alpha)	! magnitude of outgoing velocities

        sint = sint * sin(alpha)		! Spherical coordinates factor                
      enddo
      yn(s) = Y * sint
                        
      ! Determine which is the normal component of the outgoing velocity
      if (pcol(1)==1 .OR. pcol(3)==1) then
        i = 1              ! Normal component
        j = 2              ! Parallel component(not z)
      elseif (pcol(2)==1 .OR. pcol(4)==1) then
        i = 2              ! Normal component
        j = 1              ! Parallel component (not z)
      endif

      do k = 1,s
        penf(k) = Eom * (yn(k))**2
        vs(k) = sqrt(2*penf(k)*(-eme))

        rt = taus88()
        phis(k) = rt * 2 * pi        ! Calculate emission azimuthal angle of secondaries

        rt = taus88()
        thetas(k) = asin(2*rt - 1)   ! Calculate emission angle of secondaries with respect to the normal
            
        pvf(k,i) = vs(k) * cos(thetas(k))   ! Normal component of secondary velocity
        pvf(k,j) = sin(phis(k)) * vs(k) * sin(thetas(k))
        pvf(k,3) = cos(phis(k)) * vs(k) * sin(thetas(k))
      enddo

!	write(*,*) 's = ',s
!	write(*,*) 'penf1 = ', penf1
!	write(*,*) 'penf(:) = ', penf
!	write(*,*) 'vs(:) = ', vs
!	stop

      if (pcol(1) == 1 .OR. pcol(2) == 1) then
        do k = 1,s
	  pvf(k,i) = - abs(pvf(k,i))   ! Normal velocity pointing into the inside of the waveguide
        enddo
      elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
        do k = 1,s
          pvf(k,i) = abs(pvf(k,i))   ! Normal velocity pointing into the inside of the waveguide
        enddo
      endif

      deallocate (yn)
      deallocate (thetas)
      deallocate (phis)
      deallocate (penf)
      deallocate (vs)
            
    elseif (s == 0) then
      allocate (pvf(1,3))
      pvf = 0.
    endif

  endif
    
endif

!!!!!! Diagnostics and checks !!!!!!

if (s > 0) then
	
  Es = 0.
  do k = 1,s
    Es = Es + (0.5*me*(absolute(pvf(k,:)))**2)/(-e)
  enddo
   
  if (Es/penf1 > 1.0001) then
    write(*,*) 's = ',s
    write(*,*) 'energy imbalance'
    write(*,*) 'Es/penf1 = ',Es/penf1
    write(*,*) 'Es = ',Es
    write(*,*) 'penf1 = ',penf1
    write(*,*) 'r = ',r
    write(*,*) 'Pe = ',Pe
    write(*,*) 'Pr = ',Pr
    write(*,*) 'Pe+Pr = ',Pe+Pr
    write(*,*) 'pvf = ',pvf
    call exit
  endif

endif

end subroutine see


! Function to compute the absolute value/magnitude of a vector
real(long) function absolute(v)

use precision_def
implicit none
real(long) :: v(3)

absolute = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3) ) 

return
end function absolute


! Function to compute the factorial of a number
integer function factorial(n)

implicit none
integer, intent(in) :: n
integer :: i, ans

ans = 1
do i = 1,n
  ans = ans * i
enddo
factorial = ans

end function factorial

end module seem
