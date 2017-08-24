##/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:58:24 2017

@author: AP252436
"""

import MULH as lh
import os
import time
import tempfile
import numpy as np
from multiprocessing.dummy import Pool
from shutil import copyfile


# Absolute path of the project
project_path = '/Home/AP252436/Work_MULH/Parallel_batch/'
# config files
config_file = 'config.mulh'

B_xs = np.append(0, np.logspace(-3, 0, num=201))
#B_xs = np.linspace(0.12, 0.14, num=41)


NB_RUNS = 1
NB_PROCESSES = 5

def run_and_write_result(mulh):
    param_name = 'sBx'
    param = mulh.get_config_parameter(param_name)
    lh.printc(  '\n' + '_____________________________________' 
                  + '\n' + '        Current parameter :        '
                  + '\n' + '       {} = {}              '.format(param_name, param)
                  + '\n' + '_____________________________________'
                  + '\n',color='info')

    mulh.run()
    
    try:
        power, sBx, sBy, sBz = mulh.get_results()
    except:
        power = np.NaN
        sBx = np.NaN
        sBy = np.NaN
        sBz = np.NaN
    #try:
        #mulh.run()
        ## Append the results to a text file
        #power, sBx, sBy, sBz = mulh.get_results()                        # à voir si c'est les données les plus pertinentes....
        ## Write results
        #with open('/Home/AP252436/Work_MULH/Parallel_batch/Pth_vs_Bx_C3_SEY-MULH_1500electrons_analytic.csv','ba') as f_handle:   
            ## Tricks to save a row vector instead of a column one
            #np.savetxt(f_handle, np.array([float(param), power])[np.newaxis], fmt='%f', delimiter = '\t')
    #except:
        #lh.printc("MULH crashed")
        #with open('/Home/AP252436/Work_MULH/Parallel_batch/Pth_vs_Bx_C3_SEY-MULH_1500electrons_analytic.csv','ba') as f_handle:   
            ## Tricks to save a row vector instead of a column one
            #np.savetxt(f_handle, np.array([float(param), np.NaN])[np.newaxis], fmt='%f', delimiter = '\t')

    # If multipactor simulation success
    if np.isnan(power):
        power = np.NaN        
        sBx = np.NaN
        sBy = np.NaN
        sBz = np.NaN
    elif power==1: # Simulation complete without any multipac detected
        power = 0

    with open('/Home/AP252436/Work_MULH/Parallel_batch/ter_Pth_vs_Bx_C3_SEY-MULH_1500electrons.csv','ba') as f_handle:   
        ## Tricks to save a row vector instead of a column one
        np.savetxt(f_handle, np.array([float(param), power])[np.newaxis], fmt='%f', delimiter = '\t')



# Create the MULH simulations list
project_list = []

lh.printc('Creating the project directories...', color='message')
for B_x in B_xs:
    # Repeat the same simulation NB_RUNS times
    for id in range(NB_RUNS):
        # Create unique directories 
        tmpdir = tempfile.mkdtemp(prefix=os.path.join(project_path, 'tmp/'))+'/'
        tmp_path = os.path.relpath(tmpdir)
        output_path = os.path.relpath(tempfile.mkdtemp(
                        prefix=os.path.join(project_path, 'results/')))+'/'
        output_path = os.path.relpath(output_path)

        
        # Copy the config file to a unique dir
        if not os.path.isfile(os.path.join(project_path, config_file)):
            print(os.path.join(project_path, config_file))
            raise OSError('Incorrect (relative) original configuration filename')
        
        copyfile(os.path.join(project_path, config_file), 
                 os.path.join(project_path, tmp_path, config_file))


        # Create the MULH project
        MULH_project = lh.MULHcl(project_path, 
                                 config_path = os.path.join(tmp_path, config_file), 
                                 output_path = os.path.join(output_path, 'results.txt')) 

        # Modify the parameter under study in the config file   
        MULH_project.set_config_parameter(config_path = os.path.join(tmp_path, config_file), 
                                          param = 'sBx', value = B_x)
        MULH_project._get_run_command()
        with open(os.path.join(output_path, 'results.txt'),'w') as f_handle:
            f_handle.close()

        # Append the simulation to the list of simulation to perform
        project_list.append(MULH_project)
        
lh.printc('Complete !',color='message')


# Compile the Fortran code one time
project_list[0].fortran_compile()


# Run the MULH Simulations in parallel (in independant Threads)       
with Pool(processes = NB_PROCESSES) as pool:
    rs = pool.map_async(run_and_write_result, project_list, chunksize = 1)
    pool.close()
    while(True):
       if (rs.ready()):
           break
       lh.printc('Waiting for {} remaining runs'.format(rs._number_left), color='message')
       time.sleep(60)
    pool.join()
