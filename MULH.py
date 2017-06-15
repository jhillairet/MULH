##/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 11:03:24 2017

@author: AP252436
"""

import numpy as np
import os as os
import shutil

atype = 1

#==============================================================================
# ######### Geometry #########
#==============================================================================
b = 70      # Height of the waveguide in mm
a = 8       # Width of the waveguide in mm
lw = 150    # Depth/Length of the waveguide in mm. If changed make sure there are still enough cells per wavelength
lwu = 160   # Upper waveguide length (if atype==2)

NOC_lambda = 20     # Number of Cells per wavelength
NOC_lambdau = 50    # Upper number of Cells per wavelength (for atype==3)
CLfactor = 0.95     # Factor to make sure Courant-Levy condition is obeyed. Decrease CLfactor for smaller time step
Prec = 0.1          # Threshold precision (dB)

#==============================================================================
# ######### Input Wave/Field solver #########
#==============================================================================
f_I = 3.7e9     # Frequency of input wave in Hz
Pl0 = 1.0e4     # Lower power limit (watts)
Pu0 = 1.0e6    # Upper power limit (watts)

fields = 1      # Field solver. =1 FDTD, =2 analytic TE10 mode, =3 exported from other solver(need NOC_PML=1 w/ fields=3)
ramp = 1        # Increase fields slowly over ramp periods, integer
NOC_PML = 1     # Number Of Cells in PML, integer
R_max = 1e-10   # Reflection error for normally incident wave (as a fraction)
m_PML = 3       # PML grading order, integer

#==============================================================================
# ########## Static DC magnetic field (poloidal + toroidal). Gauss format, e.g. 2/10000. Minimum field allowed = 1G#########
#==============================================================================
sBx = 0.        # Toroidal magnetic field at plasma center (T)
sBxu = 0.       # Upper toroidal magnetic field at plasma center (T) (for atype==4)
sBy = 0.        # Poloidal magnetic field in waveguide (T)
sByu = 0.       # Upper poloidal magnetic field in waveguide (T)
sBz = 0.        # Radial magnetic field in waveguide (T)
sBzu = 0.       # Upper radial magnetic field in waveguide (T)

#==============================================================================
# ########## Particles #########
#==============================================================================
Np = 100        # Number of primary particles (has to be even, preferably multiples of 16)
vth = 5         # Initial energy of seed electrons (eV), integer
vra = np.sqrt(2.)  # Ratio of vth_perpendicular to vth_parallel (sqrt(2) for isotropic)
px_i = 3        # Position of seed e 1=Side walls 2=4 planes parallel to side walls 3= Randomly scattered in centered region
fmax = 4        # vth*fMax is the largest velocity represented for the sample array fSample
nsamplev = 250  # # of sample fraction values for creating a Maxwellian velocity distribution function
launch = 3      # Launch method 1=each at a different phase [0,360] 2=range of phases every 5 or 10 degrees [0,360] 3=random

#==============================================================================
# ########## Particle-Wall interaction #########
#==============================================================================
seec = 1        # SEE model, =1 Modified Vaughan, =2 Furman & Pivi, =3 Cheng w/ de Lara Re/Rr, =4 FEST3Dish?
ReRr = 1        # Contributions from elastically (Re) and rediffused (Rr) electrons. =1 de Lara, =2 CERN LHC report
E1 = 35         # First crossover
Emax = 165      # Emax(delta=max,theta=0) in eV
delta_b = 0.5   # Yield below E_0
deltamax = 2.3  # Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
p_n = 1         # p_n phenomelogical parameter in Furman & Pivi, ==2 in Cheng
kse = 1         # Energy Roughness factor, [0,2], 0=rough,2=smooth
ks = 1          # Angle Roughness factor, [0,2], 0=rough,2=smooth
Eom = 1         # Av energy of Maxwellian distribution of secondary electrons emitted (eV)
Ns = 4          # Number of secondary electrons simulated (in multiples of Np)
a_lara = 7.5e-3 # Material dependent coeff for elastic contribution in de Lara's paper
z_lara = 29     # Atomic number of coating material in Lara's fit
mat = 'copper'  # Material

#==============================================================================
# # Save outputs
#==============================================================================
psave = 0	# Save particle position and velocity every psave iteration, =0 for not saving

#==============================================================================
# #
# #
# #		END OF USER INPUTS
# #
# #*************************************************************************************************************************
#==============================================================================
os.mkdir('tempconfig')
np.savetxt('tempconfig/temp1_config.txt', (atype,b,a,lw,lwu,NOC_lambda,NOC_lambdau), fmt = '%i')
np.savetxt('tempconfig/temp2_config.txt',(CLfactor,Prec), fmt= '%f')
np.savetxt('tempconfig/temp3_config.txt', [f_I], fmt = '%e')
np.savetxt('tempconfig/temp4_config.txt',(Pl0, Pu0, fields, ramp, NOC_PML), fmt= '%i')
np.savetxt('tempconfig/temp5_config.txt',[R_max], fmt= '%.15f')
np.savetxt('tempconfig/temp6_config.txt', [m_PML], fmt = '%i')
np.savetxt('tempconfig/temp7_config.txt',(sBx, sBxu, sBy, sByu, sBz, sBzu), fmt= '%.15f')
np.savetxt('tempconfig/temp8_config.txt',(Np, vth), fmt= '%i')
np.savetxt('tempconfig/temp9_config.txt',[vra], fmt= '%.15f')
np.savetxt('tempconfig/temp10_config.txt',(px_i, fmax, nsamplev, launch, seec, ReRr, E1, Emax), fmt= '%i')
np.savetxt('tempconfig/temp11_config.txt',(delta_b, deltamax), fmt= '%.15f')
np.savetxt('tempconfig/temp12_config.txt',(p_n, kse, ks, Eom, Ns), fmt= '%i')
np.savetxt('tempconfig/temp13_config.txt',[a_lara], fmt= '%.15f')
np.savetxt('tempconfig/temp14_config.txt',[z_lara], fmt= '%i')
np.savetxt('tempconfig/temp15_config.txt',[mat], fmt= '%s', )
np.savetxt('tempconfig/temp16_config.txt',[psave], fmt= '%i')
#filenames = ['temp1_config.txt', 'temp2_config.txt','temp3_config.txt', 'temp4_config.txt','temp5_config.txt', 'temp6_config.txt','temp7_config.txt', 'temp8_config.txt','temp9_config.txt', 'temp10_config.txt','temp11_config.txt', 'temp12_config.txt','temp13_config.txt','temp14_config.txt','temp15_config.txt']
filenames = [0]*17
for i in range(1,17):
    filenames[i] = 'tempconfig/temp' + str(i) + '_config.txt'
with open('config.mulh', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            outfile.write(infile.read())
shutil.rmtree('tempconfig')