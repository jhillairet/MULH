module MULHs_dec

!***************************************************
!	VARIABLE DECLARATIONS FOR MULHs.f90
!
!

use precision_def
use def_types
implicit none

! variables related to geometry
real(long) :: a,b,lw,lwu,CLfactor
real(long), allocatable :: is(:), CLfactors(:)
integer :: NOC_lambda,NOC_PML,NOC(3),ndt,NOC_lambdau
integer, allocatable :: NOC_lambdas(:)

! Variables used for controlling power sweep
real(long) :: Prec
integer :: simno, Pl0, Pu0, Pl, Pu, P, complete, atype, ios
logical :: found, pro, multio, existe
 character*10 :: runtime
 character*8 :: rundate

! Electromagnetic wave/fields
real(long) :: eta,lambda,sBx,sBy,sBz,E0,R_max,T,gamm,w,wave(3),f_I,sBxu,sByu,sBzu
integer :: fields,ramp,m_PML

! Particles
real(long) :: vra,pvms(3),pvm(3)
integer :: Np,vth,px_i,fmax,nsamplev,launch,Ns
!real(long), allocatable :: pvm(:,:)

! variables related to SEE model
real(long) :: delta_b,deltamax,p_n,ks,kse,Eom,a_lara
integer :: seec,ReRr,E1,Emax,z_lara,material

! Saving outputs
 character(42) :: dirname
 character(34) :: counterdir
 character(50) :: makedir
 character*10 :: folder
 character*6 :: cdate, datecount
integer :: datecountn, psave
type(output_s) :: outputs

! Dummy variables
integer :: i,j,k,m,sw,ip,il,iu
 character(20) :: machine
 character(15) :: mat


end module MULHs_dec
