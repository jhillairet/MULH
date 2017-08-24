##/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 11:03:24 2017

@author: AP252436
"""

import numpy as np
import os as os
#import shutil
import subprocess as sub
import fileinput

#==============================================================================
# Messages colors
#==============================================================================
def printc(message, color='red'):
    if color in ('red','r','warning'):
        escape_code = '\x1b[41m'
    if color in ('blue', 'b','message'):
        escape_code = '\x1b[34m'
    if color in ('green', 'g','results'):
        escape_code = '\x1b[32m'
    if color in ('magenta', 'm','error'):
        escape_code = '\x1b[35m'
    if color in ('cyan', 'c','info'):
        escape_code = '\x1b[36m'
    normal_code = '\x1b[0m'
    print(escape_code+message+normal_code)


#==============================================================================
# Class MULHcl
#==============================================================================
class MULHcl(object):
    """
    MULH simulation object.
    
    IMPORTANT :
        A config.mulh file should be present in the MULH/ directory.
    
        Results will be stored in the MULH/data/results_python folder.
        """
    MULH_PATH = '/Home/AP252436/MULH/'
    BIN_PATH = MULH_PATH + 'bin/'
    counter = 0 #Counter in order to compile MULH only once
    error_counter = 0 # Counts every time a computation has to be made once again

    if not os.path.exists(MULH_PATH):
        raise OSError('Incorrect MULH folder absolute path')  
    if not os.path.exists(BIN_PATH):
        raise OSError('Incorrect bin relative path -- check that there is a bin directory in the MULH folder')
    


    def __init__(self, project_path, 
                 config_path, 
                 output_path):
        """
        Constructor
        
        Arguments:
             - project_path: absolute project path
            (- config_file: file called by MULHs.f90)
            (- results_file: name of the file where results are stored)
        """
        
        self.project_path = project_path                
#        if not os.path.exists(project_path):
#            raise OSError('Incorrect project directory absolute path')
        
        self.config_path = config_path
#        if not os.path.isfile(config_path):
#            print(config_path)
#            raise OSError('Incorrect relative config filename')
        
        self.output_path = output_path
#        if not os.path.isfile(output_path):
#            print(output_path)
#            raise OSError('Incorrect relative output filename')


            
            
        
    def create_config_file(self, config_path):
        """
        Set config.mulh
        
        Argument:
             - project_path: where config.mulh should be stored
            
        Returns:
             - config.mulh: a config file in project_path
        
        Resources:
             - MULH_PATH/names.mulh: a text file containing all variable names
             - MULH_PATH/types.mulh: a text file containing all formats
        
        """
        config = np.zeros(46, dtype = 'object')

        #****************************************************************************************************************************
        # 			USER INPUTS
        #
        #
        # Analysis type (atype)
        # 0 = single run at a user specified power, not a power sweep. Can create a video by changing psave
        # 1 = single power sweep
        # 2 = varying length analysis. Requires same inputs as 1 but also to specify lower and upper lw (in multiples of 10)
        # 3 = convergence study. Specify lower and upper NOC_lambda.
        # 4 = sensitivity to toroidal magnetic field (Bx)
        # 5 = sensitivity to poloidal magnetic field (By)
        # 6 = sensitivity to radial magnetic field (Bz)
        # 7 = threshold for different machines using different materials at different frequencies
        # 8 = single sweep w/ controlled stochasticity. Provide px2, pv2, pphases and a stock of random deviates (need enough of them)
        # 9 = single sweep w/ controlled stochasticity. Provide px2, pv2, pphases, srand, drand, irand, secrand
        # 10 = time convergence study. Specify lower and upper CLfactor will be 0.99.
        # 11 = sensitivity to a_lara unknown factor that has to be between 7e-3 and 10e-3
        # 12 = stochastic oscillation of results
        # 
        # NOTE = 2-7, 10, 11 and 12 can take a very long time so it is recommended to run them remotely (in a server)
        
        config[0] = 1       # atype
        #==============================================================================
        # ######### Geometry #########
        #==============================================================================
        config[1] = 70      # b: Height of the waveguide in mm
        config[2] = 2       # a: Width of the waveguide in mm
        config[3] = 70      # lw: Depth/Length of the waveguide in mm. If changed make sure there are still enough cells per wavelength
        config[4] = 160     # lwu: Upper waveguide length (if atype==2)
        
        config[5] = 20      # NOC_lambda: Number of Cells per wavelength
        config[6] = 50      # NOC_lambdau: Upper number of Cells per wavelength (for atype==3)
        config[7] = 0.95    # CLfactor: Factor to make sure Courant-Levy condition is obeyed. Decrease CLfactor for smaller time step
        config[8] = 0.1     # Prec: Threshold precision (dB)
        
        #==============================================================================
        # ######### Input Wave/Field solver #########
        #==============================================================================
        config[9] = 3.7e9   # f_I: Frequency of input wave in Hz
        config[10] = 1.0e4  # Pl0: Lower power limit (watts)
        config[11] = 1.0e6  # Pu0: Upper power limit (watts)
        
        config[12] = 1      # fields: Field solver. =1 FDTD, =2 analytic TE10 mode, =3 exported from other solver(need NOC_PML=1 w/ fields=3)
        config[13] = 1      # ramp: Increase fields slowly over ramp periods, integer
        config[14] = 1      # NOC_PML: Number Of Cells in PML, integer
        config[15] = 1e-10  # R_max: Reflection error for normally incident wave (as a fraction)
        config[16] = 3      # m_PML: PML grading order, integer
        
        #==============================================================================
        # ########## Static DC magnetic field (poloidal + toroidal). Gauss format, e.g. 2/10000. Minimum field allowed = 1G#########
        #==============================================================================
        config[17] = 0.     # sBx: Toroidal magnetic field at plasma center (T)
        config[18] = 0.     # sBxu: Upper toroidal magnetic field at plasma center (T) (for atype==4)
        config[19] = 0.     # sBy: Poloidal magnetic field in waveguide (T)
        config[20] = 0.     # sByu: Upper poloidal magnetic field in waveguide (T)
        config[21] = 0.     # sBz: Radial magnetic field in waveguide (T)
        config[22] = 0.     # sBzu: Upper radial magnetic field in waveguide (T)
        
        #==============================================================================
        # ########## Particles #########
        #==============================================================================
        config[23] = 100            # Np: Number of primary particles (has to be even, preferably multiples of 16)
        config[24] = 5              # vth: Initial energy of seed electrons (eV), integer
        config[25] = np.sqrt(2.)    # vra: Ratio of vth_perpendicular to vth_parallel (sqrt(2) for isotropic)
        config[26] = 3              # px_i: Position of seed e 1=Side walls 2=4 planes parallel to side walls 3= Randomly scattered in centered region
        config[27] = 4              # fmax: vth*fMax is the largest velocity represented for the sample array fSample
        config[28] = 250            # nsamplev: # of sample fraction values for creating a Maxwellian velocity distribution function
        config[29] = 3              # launch: Launch method 1=each at a different phase [0,360] 2=range of phases every 5 or 10 degrees [0,360] 3=random
        
        #==============================================================================
        # ########## Particle-Wall interaction #########
        #==============================================================================
        config[30] = 1          # seec: SEE model, =1 Modified Vaughan, =2 Furman & Pivi, =3 Cheng w/ de Lara Re/Rr, =4 FEST3Dish?
        config[31] = 1          # ReRr: Contributions from elastically (Re) and rediffused (Rr) electrons. =1 de Lara, =2 CERN LHC report
        config[32] = 35         # E1: First crossover
        config[33] = 165        # Emax: Emax(delta=max,theta=0) in eV
        config[34] = 0.5        # delta_b: Yield below E_0
        config[35] = 2.3        # deltamax: Maximum secondary electron yield (at Emax) for normal incidence (theta=0)
        config[36] = 1          # p_n: p_n phenomelogical parameter in Furman & Pivi, ==2 in Cheng
        config[37] = 1          # kse: Energy Roughness factor, [0,2], 0=rough,2=smooth
        config[38] = 1          # ks: Angle Roughness factor, [0,2], 0=rough,2=smooth
        config[39] = 1          # Eom: Av energy of Maxwellian distribution of secondary electrons emitted (eV)
        config[40] = 4          # Ns: Number of secondary electrons simulated (in multiples of Np)
        config[41] = 7.5e-3     # a_lara: Material dependent coeff for elastic contribution in de Lara's paper
        config[42] = 29         # z_lara: Atomic number of coating material in Lara's fit
        config[43] = 'copper'   # mat: Material
        
        #==============================================================================
        # ########## Save outputs ##########
        #==============================================================================
        config[44] = 0	        # psave: Save particle position and velocity every psave iteration, =0 for not saving"
        #==============================================================================
        # -------------------------END OF USER INPUTS----------------------------------
        #==============================================================================
        
        
        #==============================================================================
        # Store datas in config.mulh
        #==============================================================================
        # Read-only files were datatype and name of variables are written
        types = np.loadtxt(self.MULH_PATH + 'types.mulh', dtype=bytes).astype(str) # Trick to avoid the
        names = np.loadtxt(self.MULH_PATH + 'names.mulh', dtype=bytes).astype(str) # 'b' in front of strings

        # Write 
        outfile = open(config_path, 'w')
        for i in range(45):
            # Columns <100 : User inputs with 'types.mulh' format ; Columns >100 : Variables names
            outfile.write(str(types[i]%config[i]) + names[i] + '\n')
        outfile.close()
        
        printc(  '\n' + '_____________________________________' 
                  + '\n' + '        config.mulh created!         '
                  + '\n' + '_____________________________________'
                  + '\n',color='message')




    def set_config_parameter(self, config_path, param, value):
        """
        Changes the inital value of param by value
        
        Arguments:
            - param: the name of the parameter that you want to be changed. It must be in the list 'MULH_PATH/names.mulh'
            - value: the new desired value of param
        
        Returns:
            None
        """
        config_path = self.config_path
        types = np.loadtxt(self.MULH_PATH + 'types.mulh', dtype=bytes).astype(str) # Trick to avoid the
        names = np.loadtxt(self.MULH_PATH + 'names.mulh', dtype=bytes).astype(str) # 'b' in front of strings
        with fileinput.FileInput(config_path, inplace=True, backup='.bak') as file:
            for line in file:
                if(param in line and param == str(line.split()[1])):
                    print(str(types[file.filelineno()-1]%value + names[file.filelineno()-1])) # Replace the param with the desired value and the good format
                else:
                    print(line, end='')
    
    
    def get_config_parameter(self, param):
        """
        Check if 'param' exists or has not been deleted by error
        
        Arguments:
            - param: the parameter. It must be in the list 'MULH_PATH/names.mulh'
        
        Returns:
            - value of param
        """
        config_path = self.config_path
        
        value = []
        with fileinput.FileInput(config_path) as file:
            for line in file:
                if(param in line and param == str(line.split()[1])):
                    value = (line.split()[0])
                    print('Parameter: ' + param + ' = ' + value)
        return value
         
        
    def fortran_compile(self):
        """
        Compile the code with the 'make' command.
        
        Arguments:
            None
            
        Returns:
            Compilation status message
        """
        printc('Compiling...', color='message')
        cmd_compile = 'cd && cd ' + self.MULH_PATH + ' && make'
        try:
            sub.call(cmd_compile, shell=True)
            printc(  '\n' + '_____________________________________' 
                  + '\n' + '        Compilation complete!        '
                  + '\n' + '_____________________________________'
                  + '\n',color='message')

        except sub.CalledProcessError as e:
            printc(  '\n' + '_____________________________________' 
                  + '\n' + '         Compilation failed!         '
                  + '\n' + '             Error ' + e
                  + '\n' + '_____________________________________'
                  + '\n',color='error')
        printc('Complete !', color='message')
        
        
    def run(self):
        """
        Run the MULH modeling.
        
        Arguments:
            None
            
        Returns:
            None
        """
        try:
            env = os.environ
            with sub.Popen(self._get_run_command(), shell=True, env=env, 
                       stdout=sub.PIPE, stderr=sub.PIPE, universal_newlines=True) as p:
                for lines in p.stdout:
                    print(lines, end=' ') # Print MULH messages
#            print('Debug renvoi None MULH.py l 311 ' + p.returncode)
        except OSError as e:
            printc('Error ! ' + e, color='error')
#        print('MULH ligne 315')
        
    
    def _get_run_command(self):
        """
        Define bash command that will run the MULH code
        
        Arguments:
            None
            
        Returns:
            command line
        """
        cmd = 'cd && cd ' + self.BIN_PATH + ' && ./MULH ' + self.config_path + ' ' + self.output_path
        return(cmd)
   
        
        
    def get_results(self):
        """
        Returns MULH run results
        
        Arguments:
            None
        
        Returns:
            power: array of breakdown power
            sBx, sBy, sBz: toroidal, poloidal and radial static magnetic fields
        """
        power, sBx, sBy, sBz = None, None, None, None
        if os.path.isfile(self.output_path):
            power, sBx, sBy, sBz = np.loadtxt(self.output_path, 
                               skiprows=0, 
                               unpack=True)
        self.power = power
        self.sBx = sBx
        self.sBy = sBy
        self.sBz = sBz
        return power, sBx, sBy, sBz

               
        
if __name__ == "__main__":  
    # Absolute path of the project
    project_path = '/Home/AP252436/Work_MULH/'
    
    # Output path
    results_file = 'results/results.txt'
    output_path = os.path.join(project_path, results_file)
    
    # Path to the .mulh file     
    config_file = 'config.mulh'
    config_path = os.path.join(project_path, config_file)
    printc('RESULTS PATH : ' + output_path)
    # Run the MULH Simulation 
    mulh = MULHcl(project_path, config_path, output_path)

    # Create the config file
    mulh.create_config_file(config_path)
    
    
    # Give the parameter name
    parameter_name = 'sBx'

    # Compile only if it is the first execution of the program
    if(mulh.counter == 0):
        mulh.fortran_compile()
    mulh.counter = mulh.counter + 1
    printc("Compteur d'exécutions : " + str(mulh.counter),color='message') # Debug
    try:
        mulh.run()
        power, sBx, sBy, sBz = mulh.get_results()
    except ValueError:
        mulh.error_counter = mulh.error_counter + 1 
        printc("ERREUR !! Compteur d'erreurs : " + str(mulh.error_counter),color='warning')
        mulh.run()
        power, sBx, sBy, sBz = mulh.get_results()

            
    printc('Breakdown power : ' + str(power) + 'W'
            + '\n' + 'Toroidal magnetic field : ' + str(sBx) + 'T'
            + '\n' + 'Poloidal magnetic field : ' + str(sBy) + 'T'
            + '\n' + 'Radial magnetic field : ' + str(sBz) + 'T', color='results')
        
        # Appending the results to a text file
    with open(os.path.join(project_path, 'RESULTS.txt'),'ba') as f_handle:
        np.savetxt(f_handle, [power, sBx])
        
        
        
        
        
        
        
        
