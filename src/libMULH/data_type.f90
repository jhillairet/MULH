module data_type

use precision_def
use def_types
implicit none

! variables related to geometry
real(long) :: a,b,lw,d4x(4),CLfactor,dx,dy,dz,rLH,r0,dt
integer :: NOC_lambda,NOC_PML,NOC(3),ndt
real(long), allocatable :: x(:),y(:),z(:),x3(:)

! variables related to electromagnetic wave/fields
real(long) :: eta,lambda,beta,sB(5),d,E0,R_max,T,gamm,w,wave(4),f_I
integer :: fields,ramp,m_PML
real(long), allocatable :: Exn(:,:,:,:),Eyn(:,:,:,:),Ezn(:,:,:,:),Hxn(:,:,:,:),Hyn(:,:,:,:),Hzn(:,:,:,:), &
			Ex2(:,:,:),Ey2(:,:,:),Ez2(:,:,:),Hx2(:,:,:),Hy2(:,:,:),Hz2(:,:,:),Ex1(:,:,:), &
			Ey1(:,:,:),Ez1(:,:,:),Hx1(:,:,:),Hy1(:,:,:),Hz1(:,:,:),Hx0(:,:,:),Hy0(:,:,:),Hz0(:,:,:), &
			CE(:,:,:,:,:),DH(:,:,:,:,:),CEsigz(:,:,:,:),DHsigz(:,:,:,:),ExPML(:,:,:,:), &
			EyPML(:,:,:,:),EzPML(:,:,:,:),HxPML(:,:,:,:),HyPML(:,:,:,:),HzPML(:,:,:,:), &
			Jx(:,:,:),Jy(:,:,:),Jz(:,:,:)
real(long), allocatable :: Ext(:,:,:)

! variables related to particles
real(long) :: vra,pB0(3),pB1(3),pB2(3),pE1(3),pE2(3),part(28),Eri(3,3),Hri(3,3),pvm(3)
integer :: Np,psave,vth,px_i,fmax,nsamplev,launch,Ns,pstart,Nt,pcol(4),pwall(4), &
		ps_c,lost(2),launched,pactc,ps_c0,pact,Eii(3,3),Eip(3,3), &
		Hii(3,3),Hip(3,3),s
real(long), allocatable :: px1(:,:),px2(:,:),pv1(:,:),pv2(:,:),pphases(:,:),phi0s(:,:), px2t(:,:), pv2t(:,:), &
			pactts(:)
integer, allocatable :: pstat(:), pactps(:)

! variables related to SEE model
real(long) :: delta_b,deltamax,p_n,ks,kse,Eom,a_lara,z_lara
integer :: seec,ReRr,E1,Emax,E_0,E_0p
real(long), allocatable :: SEYvsE(:,:,:)
 character(15) :: mat

! Flags, counters and others
 character(50) :: makedir
 character*8 :: itn
 character*10 :: runtime
integer :: i,j,k,l,m,n,it
logical :: equals,full
real (long) :: time

! User defined types, mostly for packing info into single array for ease of transport  
type(sey_t) :: sey
type(geo_t) :: geo
type(seed_t) :: seed
type(particle_t) :: particle

! Testing & debuggin
real(long) :: ans, qans, yn(1),lnbeta,inbeta,alpha,betalog
real(long), allocatable :: As(:,:,:)
integer :: pact0, ranpoi(5000), ll
end module
