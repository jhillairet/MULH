program MULHs

!********************************************************************************
!	MUltipactor in Lower Hybrid antenna waveguides
!			MULH
!
! This code does a power sweep to calculate at which point the multipactor
! starts to develop inside a single rectangular waveguide.
! EM fields can be resolved analytically for the TE10 mode of a rectangular
! waveguide, or the FDTD method with a TE10 drive could be used, or other
! fields could be imported (e.g. from COMSOL).
! The particle wall interaction is modeled using G. Chen, L. Liu (2011)
! ejection-collection algorithm but the actual method for deciding between
! absorption, elastic or inelastic collision or secondaries has been slightly
! modified to match FEST3D results. The option to use the Furman & Pivi model
! is also made available.
!
! The Ecuyer Taus random number generator from Alan Miller's library is used.
! It has a period of 2^88 and the poisson deviates are generated using this
! generator and the waiting time method for means < 15.
!
! AUTHOR: Manaure Francisquez
! Thayer School of Engineering
! Dartmouth College
! 14 Engineering Dr.
! Hanover, NH 03755
! USA
!
! contact email: mfrancisquez@gmail.com
!
! Spring 2011 - Spring 2012
!
!
!

use MULHs_dec
use multipac
use constants
use reshape_array
implicit none
character(42) :: output_path
character(42) :: project_path = '../../Work_MULH'
 call date_and_time(DATE=rundate,TIME=runtime)
!****************************************************************************************************************************
! 			USER INPUTS
!
!
! Analysis type (atype)
! 0 = single run at a user specified power, not a power sweep. Can create a video by changing psave
! 1 = single power sweep
! 2 = varying length analysis. Requires same inputs as 1 but also to specify lower and upper lw (in multiples of 10)
! 3 = convergence study. Specify lower and upper NOC_lambda.
! 4 = sensitivity to toroidal magnetic field (Bx)
! 5 = sensitivity to poloidal magnetic field (By)
! 6 = sensitivity to radial magnetic field (Bz)
! 7 = threshold for different machines using different materials at different frequencies
! 8 = single sweep w/ controlled stochasticity. Provide px2, pv2, pphases and a stock of random deviates (need enough of them)
! 9 = single sweep w/ controlled stochasticity. Provide px2, pv2, pphases, srand, drand, irand, secrand
! 10 = time convergence study. Specify lower and upper CLfactor will be 0.99.
! 11 = sensitivity to a_lara unknown factor that has to be between 7e-3 and 10e-3
! 12 = stochastic oscillation of results
! 
! NOTE = 2-7, 10, 11 and 12 can take a very long time so it is recommended to run them remotely (in a server)
!

open(unit=30, file=trim(project_path)//'/config.mulh',status='old')
read(30,*) atype

!!!!!!!!! Geometry !!!!!!!!!
read(30,*) b		! Height of the waveguide in mm
read(30,*) a		! Width of the waveguide in mm
read(30,*) lw		! Depth/Length of the waveguide in mm. If changed make sure there are still enough cells per wavelength
read(30,*) lwu		! Upper waveguide length (if atype==2)
read(30,*) NOC_lambda	! Number of Cells per wavelength
read(30,*) NOC_lambdau	! Upper number of Cells per wavelength (for atype==3)
read(30,*) CLfactor	! Factor to make sure Courant-Levy condition is obeyed. Decrease CLfactor for smaller time step
read(30,*) Prec		! Threshold precision (dB)

eta = mu_0 * c
!!!!!!!!! Input Wave/Field solver !!!!!!!!!
read(30,*) f_I	! Frequency of input wave in Hz
lambda = c/f_I	! Wavelength
read(30,*) Pl0	! Lower power limit (watts)
read(30,*) Pu0	! Upper power limit (watts)

read(30,*) fields	! Field solver. =1 FDTD, =2 analytic TE10 mode, =3 exported from other solver(need NOC_PML=1 w/ fields=3)
read(30,*) ramp		! Increase fields slowly over ramp periods, integer
read(30,*) NOC_PML	! Number Of Cells in PML, integer
read(30,*) R_max	! Reflection error for normally incident wave (as a fraction)
read(30,*) m_PML	! PML grading order, integer

!!!!!!!!!! Static DC magnetic field (poloidal + toroidal). Gauss format, e.g. 2/10000. Minimum field allowed = 1G!!!!!!!!!
read(30,*) sBx		! Toroidal magnetic field at plasma center (T)
read(30,*) sBxu		! Upper toroidal magnetic field at plasma center (T) (for atype==4)
read(30,*) sBy		! Poloidal magnetic field in waveguide (T)
read(30,*) sByu		! Upper poloidal magnetic field in waveguide (T)
read(30,*) sBz		! Radial magnetic field in waveguide (T)
read(30,*) sBzu		! Upper radial magnetic field in waveguide (T)

!!!!!!!!!! Particles !!!!!!!!!
read(30,*) Np		! Number of primary particles (has to be even, preferably multiples of 16)
read(30,*) vth		! Initial energy of seed electrons (eV), integer
read(30,*) vra		! Ratio of vth_perpendicular to vth_parallel (sqrt(2) for isotropic)
read(30,*) px_i		! Position of seed e 1=Side walls 2=4 planes parallel to side walls 3= Randomly scattered in centered region
read(30,*) fmax		! vth*fMax is the largest velocity represented for the sample array fSample
read(30,*) nsamplev	! # of sample fraction values for creating a Maxwellian velocity distribution function
read(30,*) launch	! Launch method 1=each at a different phase [0,360] 2=range of phases every 5 or 10 degrees [0,360] 3=random

!!!!!!!!!! Particle-Wall interaction !!!!!!!!!
read(30,*) seec		! SEE model, =1 Modified Vaughan, =2 Furman & Pivi, =3 Cheng w/ de Lara Re/Rr, =4 FEST3Dish?
read(30,*) ReRr		! Contributions from elastically (Re) and rediffused (Rr) electrons. =1 de Lara, =2 CERN LHC report
read(30,*) E1		! First crossover
read(30,*) Emax		! Emax(delta=max,theta=0) in eV
read(30,*) delta_b	! Yield below E_0
read(30,*) deltamax	! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
read(30,*) p_n		! p_n phenomelogical parameter in Furman & Pivi, ==2 in Cheng
read(30,*) kse		! Energy Roughness factor, [0,2], 0=rough,2=smooth
read(30,*) ks		! Angle Roughness factor, [0,2], 0=rough,2=smooth
read(30,*) Eom		! Av energy of Maxwellian distribution of secondary electrons emitted (eV)
read(30,*) Ns		! Number of secondary electrons simulated (in multiples of Np)
read(30,*) a_lara	! Material dependent coeff for elastic contribution in de Lara's paper
read(30,*) z_lara	! Atomic number of coating material in Lara's fit
read(30,*) mat		! Material

! Save outputs
read(30,*) psave	! Save particle position and velocity every psave iteration, =0 for not saving
read(30,*) output_path	! Where batch results will be stored

close(unit=30)
!
!
!		END OF USER INPUTS
!
!*************************************************************************************************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create arrays containing the different values used in a given study
!
b = b/1000.              ! Height of the waveguide in meters
a = a/1000.               ! Width of the waveguide in meters
lw = lw/1000.             ! Depth/Length of the waveguide in meters
lwu = lwu/1000.             ! Upper Depth/Length of the waveguide in meters
if (atype == 0 .OR. atype == 1 .OR. atype == 8 .OR. atype == 9) then
  j = 1
elseif (atype == 2 .OR. atype == 4 .OR. atype == 5 .OR. atype == 6) then

  if (atype == 2) then
    iu = anint(lwu*1000)
    il = anint(lw*1000)
    write(*,*) '***** VARYING LENGTH ANALYSIS *****'
  elseif (atype == 4) then
    iu = anint(sBxu*10000)
    il = anint(sBx*10000)
    write(*,*) '***** SENSITIVITY TO TOROIDAL B ANALYSIS *****'
  elseif (atype == 5) then
    iu = anint(sByu*10000)
    il = anint(sBy*10000)
    write(*,*) '***** SENSITIVITY TO POLOIDAL B ANALYSIS *****'
  elseif (atype == 6) then
    iu = anint(sBzu*10000)
    il = anint(sBz*10000)
    write(*,*) '***** SENSITIVITY TO RADIAL B ANALYSIS *****'
  endif
  
  if (il == 0) il = 1

  ! Calculate the order of magnitude of lower parameter
  i = 0
  k = il
  do
    k = k/10
    if (k >= 1) then
      i = i + 1
    else
      exit
    endif
  enddo

  if (iu/il >= 10) then

    ! Calculate the order of magnitude of upper parameter
    ip = 0
    k = iu
    do
      k = k/10
      if (k >= 1) then
        ip = ip + 1
      else
        exit
      endif
    enddo

    ! Determine the length of the array containing all the parameters
    j = 10 - il/(10**i) + 9*(ip-i-1) + iu/(10**ip)
    allocate(is(j))

    ! Populate the array containing all parameters
    m = il
    ip = 0
    do k=1,j

      ip = ip + 1
      is(k) = m + (ip-1)*(10**i)
      if (is(k)/(10**i) >= 9) then
	i = i + 1
	ip = 0
	m = 10**i
      endif

    enddo

  else

    ! Determine the length of the array containing all the lengths
    j = (iu - il)*2/(10**i) + 1
    allocate(is(j))

    ! Populate the array containing all the lengths
    do k=1,j
      is(k) = il + (k-1)*5*(10**(i-1))
    enddo

  endif

elseif (atype == 3) then

  i = NOC_lambda
  j = 1
  do while (i < NOC_lambdau)
    i = i * 2
    j = j + 1
  enddo

  allocate (NOC_lambdas(j))

  NOC_lambdas(1) = NOC_lambda
  do i = 2,j-1
    NOC_lambdas(i) = NOC_lambdas(i-1)*2
  enddo
  NOC_lambdas(j) = NOC_lambdau

  write(*,*) '***** CONVERGENCE STUDY *****'

elseif (atype == 7) then

  write(*,*) '***** ANALYSIS OF DIFFERENT MACHINES *****'
  ! Import specs of different machines
  call file_row_count('../data/LHs.txt',i)
  open(unit=22,file='../data/LHs.txt',status='old')  

  call file_row_count('../data/machines.txt',j)
  open(unit=23,file='../data/machines.txt',status='old')  

  if (i /= j) then
    write(*,*) 'Number of records in LHs.txt and machines.txt do not coincide. Please revise.'
    call exit
  endif

elseif (atype == 10) then

  j = (0.95 - CLfactor)*20 + 2
  allocate(CLfactors(j))

  do i = 1,j-1
    CLfactors(i) = CLfactor + (i-1)*0.05
  enddo
  CLfactors(j) = 0.99

  write(*,*) '**** CONVERGENCE STUDY OF TIME STEP SIZE *****'

elseif (atype == 11) then

  j = 16
  allocate (is(j))
  is = (/7e-3,7.2e-3,7.4e-3,7.6e-3,7.8e-3,8e-3,8.2e-3,8.4e-3,8.6e-3,8.8e-3,9e-3,9.2e-3,9.4e-3,9.6e-3,9.8e-3,10e-3/)

  write(*,*) '**** SENSITIVITY TO a_Lara FACTOR STUDY *****'

elseif (atype == 12) then

  j = 100
 
endif


do sw = 1,j
  
  if (atype == 2) then
    lw = is(sw)/1000.
  elseif (atype == 3) then
    NOC_lambda = NOC_lambdas(sw)
  elseif (atype == 4) then
    sBx = is(sw)/10000.
  elseif (atype == 5) then
    sBy = is(sw)/10000.
  elseif (atype == 6) then
    sBz = is(sw)/10000.
  elseif (atype == 7) then
    read(22,*) b,a,f_I,material
    read(23,*) machine
    if (material == 0) then 	! Copper
      E1 = 35                		! First crossover
      Emax = 165             		! Emax(delta=max,theta=0) in eV
      delta_b = 0.5          		! Yield below E_0
      deltamax = 2.3         		! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
      z_lara = 29			! Atomic number
      mat = 'copper'
    elseif (material == 1) then ! Stainless steel
      E1 = 30                		! First crossover
      Emax = 300             		! Emax(delta=max,theta=0) in eV
      delta_b = 0.5          		! Yield below E_0
      deltamax = 2         		! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
      z_lara = 26			! Atomic number
      mat = 'stainless steel'
    elseif (material == 2) then ! Carbon
      E1 = 67                		! First crossover
      Emax = 250             		! Emax(delta=max,theta=0) in eV
      delta_b = 0.5          		! Yield below E_0
      deltamax = 0.88         		! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
      z_lara = 6			! Atomic number
      mat = 'carbon'
    elseif (material == 3) then	! Titanium
      E1 = 30                		! First crossover
      Emax = 240             		! Emax(delta=max,theta=0) in eV
      delta_b = 0.5          		! Yield below E_0
      deltamax = 1.9        		! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
      z_lara = 22			! Atomic number
      mat = 'titanium'
    elseif (material == 4) then	! Gold
      E1 = 100                		! In the case of gold this is actually E_0 to be used in Vaughan's formula. No E1.
      Emax = 800             		! Emax(delta=max,theta=0) in eV
      delta_b = 0.5          		! Yield below E_0
      deltamax = 0.53         		! Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
      z_lara = 79			! Atomic number
      mat = 'gold'
    endif
  elseif (atype == 10) then
    CLfactor = CLfactors(sw)
  elseif (atype == 11) then
    a_lara = is(sw)
  endif

  !**************************************************
  !
  ! PRINT SCENARIO TO SCREEN
  !
  !
  write(*,8)'Began simulation on:',trim(rundate(7:8)),'/',trim(rundate(5:6)),'/',trim(rundate(3:4)), &
	'@',trim(runtime(1:2)),':',trim(runtime(3:4)),':',trim(runtime(5:6))
  write(*,*)'=>   SCENARIO   <='
  if (atype == 7) write(*,*) 'Machine: ',adjustl(machine)
  write(*,1)'Waveguide width:',a*1e3,'mm'
  write(*,2)'Waveguide height:',b*1e3,'mm'
  write(*,3)'Waveguide length:',INT(anINT(lw*1000.)),'mm'
  if (atype == 3) write(*,'(A29,I3)') 'No. of cells per wavelength: ',NOC_lambda
  if (atype == 10) write(*,*) 'Courant-Levy factor: ',CLfactor
  write(*,4)'Operation frequency:',REAL(f_I/1e9),'GHz'
  if (atype == 11) write(*,*)'a_lara: ',a_lara
  if (sBx /= 0.) write(*,*) 'Toroidal magnetic field at plasma center: ', sBx
  if (sBy /= 0.) write(*,*) 'Poloidal magnetic field: ', sBy
  if (sBz /= 0.) write(*,*) 'Radial magnetic field: ', sBz
  write(*,*)'_____________________________________'
  write(*,*)'SEY of ',trim(mat)
  write(*,5)'E1 =',E1,'eV'
  write(*,6)'Emax =',Emax,'eV'
  write(*,7)'Delta_max =',deltamax
  write(*,*)'_____________________________________'
  !
  !
  !***************************************************


  !!!!!!!!!!!!!!!! Counters, flags & others !!!!!!!!!!!!!!!!!
  simno = 0		! Simulation number (within one power sweep)
  pro = .TRUE.		! Boolean indicating whether to proceed or not (due to error)
  found = .FALSE.	! Flag indicating whether threshold has been found
  multio = .FALSE.	! Indicate if multipactor occurred at all during power sweep
  existe = .FALSE.	! Folder w/ counter.txt for this study already exists

  Pl = Pl0		! Lower Power limit
  Pu = Pu0		! Upper Power limit
  P = Pl		! Current Power delivered by waveguide (watts)

  pvms = 0.		! Average energy of all particles in time

  ! Pack outputs. Indicate which outputs to save.
  outputs%psave = psave

  !!!!!!!!!!!!!!!!! DIRECTORY WHERE RESULTS ARE STORED !!!!!!!!!!!!!!!!!!!!!!
  if (sw == 1) then

    counterdir = '../data/results'
    if (atype == 2) then
      inquire(file='../data/results/length/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/length'
    elseif (atype == 3) then
      inquire(file='../data/results/convergence/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/convergence'
    elseif (atype == 4) then
      inquire(file='../data/results/bt/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/bt'
    elseif (atype == 5) then
      inquire(file='../data/results/bp/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/bp'
    elseif (atype == 6) then
      inquire(file='../data/results/br/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/br'
    elseif (atype == 7) then
      inquire(file='../data/results/freq/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/freq'
    elseif (atype == 10) then
      inquire(file='../data/results/timeconvergence/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/timeconvergence'
    elseif (atype == 11) then
      inquire(file='../data/results/a_lara/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/a_lara'
    elseif (atype == 12) then
      inquire(file='../data/results/oscillations/counter.txt', exist = existe)
      if (existe) counterdir = '../data/results/oscillations'
    endif

    open(unit=2,file=trim(counterdir)//'/counter.txt',status='old')	! Open counter with info of last run
    read(2,'(A6)') datecount	! Read date of last run as a character
    read(2,*) datecountn	! Read in number of last run in the day it was run
    call date_and_time(DATE=folder)
    cdate = folder(3:8)	! Store current date with format yymmdd

    ! Compare dates and if equal add to counter, else change date and restart counter
    if (cdate==datecount) then
      datecountn = datecountn + 1
    else
      datecount = cdate
      datecountn = 1
    endif

    rewind(unit=2)
    write(2,'(A6)') datecount	! Replace old date with new date
    write(2,*) datecountn 	! Update counter
    close(unit=2)

  else

    if (.NOT. existe) then
      ! Copy counter.txt into the directory used by atype
      call system('cp -n ../data/results/counter.txt '//trim(makedir)//'/counter.txt')
      existe = .TRUE.
    endif

    open(unit=2,file=trim(makedir)//'/counter.txt',status='old')	! Open counter with info of last run
    read(2,'(A6)') datecount	! Read date of last run as a character
    read(2,*) datecountn	! Read in number of last run in the day it was run
    call date_and_time(DATE=folder)
    cdate = folder(3:8)	! Store current date with format yymmdd

    ! Compare dates and if equal add to counter, else change date and restart counter
    if (cdate==datecount) then
      datecountn = datecountn + 1
    else
      datecount = cdate
      datecountn = 1
    endif

    rewind(unit=2)
    write(2,'(A6)') datecount	! Replace old date with new date
    write(2,*) datecountn 		! Update counter
    close(unit=2)

  endif

  write(cdate,'(I6)') datecountn	! Convert datecountn into a string named cdate
  cdate = adjustl(cdate)
  folder = datecount(1:6)//"_"//cdate	! Name of the folder where results of this run will be stored
  if (atype == 0 .OR. atype == 1) then	! Single simulation or single sweep can go directly to date-named folder
    dirname = '../data/results/'//trim(folder)
  elseif (atype == 2) then	! Create folder for results from varying length analysis
    makedir = '../data/results/length'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/length/'//trim(folder)
  elseif (atype == 3) then	! Create folder for results from convergence study
    makedir = '../data/results/convergence'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/convergence/'//trim(folder)
  elseif (atype == 4) then	! Create folder for results from sensitivity to toroidal B analysis
    makedir = '../data/results/bt'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/bt/'//trim(folder)
  elseif (atype == 5) then	! Create folder for results from sensitivity to poloidal B analysis
    makedir = '../data/results/bp'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/bp/'//trim(folder)
  elseif (atype == 6) then	! Create folder for results from sensitivity to radial B analysis
    makedir = '../data/results/br'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/br/'//trim(folder)
  elseif (atype == 7) then	! Create folder for results from analysis of different machines
    makedir = '../data/results/freq'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/freq/'//trim(folder)
  elseif (atype == 10) then
    makedir = '../data/results/timeconvergence'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/timeconvergence/'//trim(folder)
  elseif (atype == 11) then
    makedir = '../data/results/a_lara'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/a_lara/'//trim(folder)
  elseif (atype == 12) then
    makedir = '../data/results/oscillations'
    call system('mkdir -p '//trim(makedir))
    dirname = '../data/results/oscillations/'//trim(folder)
  endif
  call system('mkdir -p '//trim(dirname))
 
  open(unit=1,file=trim(dirname)//'/MULH_power_results.txt',status='new')

  open(unit=4,file=trim(dirname)//'/MULH_time_results.txt',status='unknown')

  open(unit=24,file=trim(dirname)//'/energies.txt',status='new')

!write(1,16)'Began simulation on:',trim(rundate(7:8)),'/',trim(rundate(5:6)),'/',trim(rundate(3:4)), &
!	'@',trim(runtime(1:2)),':',trim(runtime(3:4)),':',trim(runtime(5:6))
!write(1,*)'************* SCENARIO ****************'
  if (atype == 7) write(1,*) 'Machine: ',trim(machine)
  write(1,9)'Waveguide width:',a*1e3,'mm'
  write(1,10)'Waveguide height:',b*1e3,'mm'
  write(1,11)'Waveguide length:',INT(anINT(lw*1000.)),'mm'
  if (atype == 3) write(1,'(A29,I3)') 'No. of cells per wavelength: ',NOC_lambda
  if (atype == 10) write(1,'(A21,F4.2)') 'Courant-Levy factor: ',CLfactor
  write(1,12)'Operation frequency:',REAL(f_I/1e9),'GHz'
  if (atype == 11) write(1,'(A8,ES10.3E3)') 'a_lara: ',a_lara
  if (sBx /= 0.) write(1,'(A41,1X,ES8.1E3)') 'Toroidal magnetic field at plasma center:', sBx
  if (sBy /= 0.) write(1,'(A24,1X,ES8.1E3)') 'Poloidal magnetic field:', sBy
  if (sBz /= 0.) write(1,'(A22,1X,ES8.1E3)') 'Radial magnetic field:', sBz
  write(1,'(A37)')'_____________________________________'
  write(1,'(A6,1X,A15)')'SEY of',adjustl(mat)
  write(1,13)'E1 =',E1,'eV'
  write(1,14)'Emax =',Emax,'eV'
  write(1,15)'Delta_max =',deltamax
  write(1,'(A37)')'_____________________________________'
  write(1,'(A9,3X,A9,3X,A8)')'Power (W)','Breakdown','Energies'

!!!!!!!!!!!!!!!!! Print scenario to a file !!!!!!!!!!!!!!!!!!
! Order has to match the read order in MULH.f90 ! Check if changed.
  open (unit=11,file=trim(dirname)//'/scenario.txt',status='unknown')
  write(11,*) a
  write(11,*) b
  write(11,*) lw
  write(11,*) NOC_lambda
  write(11,*) CLfactor
  write(11,*) ramp
  write(11,*) NOC_PML
  write(11,*) R_max
  write(11,*) m_PML
  write(11,*) f_I
  write(11,*) fields
  write(11,*) sBx
  write(11,*) sBy
  write(11,*) sBz
  write(11,*) Np
  write(11,*) vth
  write(11,*) vra
  write(11,*) px_i
  write(11,*) fmax
  write(11,*) nsamplev
  write(11,*) launch
  write(11,*) seec
  write(11,*) mat
  write(11,*) ReRr
  write(11,*) E1
  write(11,*) Emax
  write(11,*) delta_b
  write(11,*) deltamax
  write(11,*) p_n
  write(11,*) kse
  write(11,*) ks
  write(11,*) Eom
  write(11,*) Ns
  write(11,*) a_lara
  write(11,*) z_lara
  close(unit=11)

  !*****************************************************************************
  !
  !	RUN SIMULATION
  !
  !
  do while (.NOT. found)

    simno = simno + 1
    i = 0
    ios = 0

    ! Run analysis
    call multipactor(atype,dirname,P,outputs,pro,simno,complete)

    if (.NOT. pro) then
      write(*,*) 'MULHs pro error. Check code/simulation. Terminating.'
      exit
    endif

    ! Calculate the mean energy of all particles over time
    rewind(unit=24)
    do
      if (ios /= 0) then
	exit
      else
        read(24,'(F20.16,1X,F20.16,1X,F20.16)',iostat=ios) pvm(1), pvm(2), pvm(3)
	pvms(1) = pvms(1) + pvm(1)
	pvms(2) = pvms(2) + pvm(2)
	pvms(3) = pvms(3) + pvm(3)
        i = i + 1
      endif
    enddo
    rewind(unit=24)
    pvms = pvms/i


    if (complete == 2) then

      ! MULH power results. Power BD/NoBD | Mean Energies
      write(1,17) P,'NoBD | Mean energies (eV):',pvms(1),pvms(2),pvms(3)

      if (atype == 0) exit	! Exit power sweep when only one simulation is desired

      if (10*log10(real(Pu)/P) <= Prec) then
        found = .TRUE.
      else

        Pl = P
        if (P*2 < Pu) then
	  P = P*2
        else
	  P = P + (Pu-P)*0.5
        endif
      endif

    elseif (complete == 1) then

      ! MULH power results. Power BD/NoBD | Mean Energies
      write(1,17) P,'BD   | Mean energies (eV):',pvms(1),pvms(2),pvms(3)

      if (atype == 0) exit	! Exit power sweep when only one simulation is desired

      multio = .TRUE.

      if (10*log10(real(P)/Pl) <= Prec) then
        found = .TRUE.
      else
        Pu = P
        P = Pl + (P-Pl)*0.5;
      endif

    endif

  enddo

  close(unit=1)
  close(unit=4)
  close(unit=24)

  ! get ride of energies.txt file
  call system('rm '//trim(dirname)//'/energies.txt')


  write(*,*)'_____________________________________'
  open(unit=31, file=trim(project_path)//'/'//trim(output_path),status='old', position='append')
  if (.NOT. multio) then
    write(*,*) ' **** NO BREAKDOWN FOUND **** '
    write(31,*) ' **** NO BREAKDOWN FOUND **** ', sBx, sBy, sBz
  elseif (complete == 1) then
    write(*,*) 'Multipactor Threshold (W): ',Pl
    write(31,*) Pl, sBx, sBy, sBz
  elseif (complete == 2) then
    write(*,*) 'Multipactor Threshold (W): ',P
    write(31,*) P, sBx, sBy, sBz
  endif
  close(unit = 31)
  call date_and_time(DATE=rundate,TIME=runtime)
  write(*,16) '===> Terminee le ',trim(rundate(7:8)),'/',trim(rundate(5:6)),'/',trim(rundate(3:4)), &
	  '@',trim(runtime(1:2)),':',trim(runtime(3:4)),':',trim(runtime(5:6))
  write(*,*) 'Results written to ',trim(dirname)


enddo


  !**************************************************************************************
  !
  !				LIST OF FILES OPENED/USED
  !
  ! (CI) = Closed Immediately so unit number can be reused
  !
  ! unit=11 scenario.txt: temporarily saving the scenario for current power sweep
  ! unit=2 counter.txt: date and run number of last power sweep. Replace with new one (CI)
  ! unit=1 MULH_power_results: storing the results from each Power run within a power sweep
  ! unit=22 LHs.txt: import specs of different machines (CI)
  ! unit=23 machines.txt: import names of different machines (CI)
  ! unit=11 scenario.txt: loading the scenario into subroutine multipactor
  ! unit=14 px2.txt: loading initial electron locations previously generated (for comparison) (CI)
  ! unit=15 pv2.txt: loading initial electron velocities previously generated (for comparison) (CI)
  ! unit=16 phases0.txt: loading initial electron launch phases previously generated (for comparison)
  ! unit=17 srand.txt: random variable s (number of secondaries) previously generated (for comparison)
  ! unit=18 drand.txt: random variable r (decision) previously generated (for comparison)
  ! unit=19 irand.txt: random numbers used in inelastic collisions previously generated (for comparison)
  ! unit=20 secrand.txt: random numbers used in calculating velocity of secondaries previously generated (for comparison)
  ! unit=21 rstock.txt: stock of random numbers previously generated (for comparison)
  ! unit=3 pxn.txt: storing particle locations during iteration n (CI)
  ! unit=3 pvn.txt: storing particle velocities during iteration n (CI)
  ! unit=4 MULH_time_results.txt: storing number of particles in time for each simulation in the power sweep
  ! unit=7 pxn.txt: storing particle locations during iteration n (CI)
  ! unit=8 pvn.txt: storing particle velocities during iteration n (CI)
  ! unit=10 SEYvsE.txt: storing each impact energy and each SEY (delta)
  ! unit=13 setup.txt: storing the setup used in this simulation inorder to plot particles in time
  ! unit=24 mean energy along each direction in time
  ! unit=25 enorm from PICCOLO (COMSOL) and comsolfields.m
  ! unit=30 config.mulh: load all constants MULH needs in order to work
  ! unit=31 date_hour.txt: save BD power 


1 format(1X,A16,1X,F5.2,A2)
2 format(1X,A17,1X,F5.1,A2)
3 format(1X,A17,1X,I6,A2)
4 format(1X,A20,1X,F4.2,A3)
5 format(3X,A4,1X,I3,A2)
6 format(3X,A6,1X,I3,A2)
7 format(3X,A11,1X,F4.2)
8 format(1X,A20,1X,A2,A1,A2,A1,A2,1X,A1,1X,A2,A1,A2,A1,A2)
9 format(A16,1X,F5.2,A2)
10 format(A17,1X,F5.1,A2)
11 format(A17,1X,I6,A2)
12 format(A20,1X,F4.2,A3)
13 format(A4,1X,I3,A2)
14 format(A6,1X,I3,A2)
15 format(A11,1X,F4.2)
16 format(A20,1X,A2,A1,A2,A1,A2,1X,A1,1X,A2,A1,A2,A1,A2)
17 format(I7,1X,A26,1X,F7.3,1X,F6.3,1X,F6.3)

end program MULHs
