######################################### Imports #########################################
import sys
import os
from os import path
import datetime
from distance_module import run_distance_calculation
######################################### Run Main #########################################
def run(verbose, outdir, sample, annotation, 
        cpus=1, window=1500, experimental_fimo=False, pre_scan=None): 

######################################## Distance Calculation ######################################### 
    if verbose == True:
        print('--------------Calculating Distances--------------')
        print('Start time: %s' % str(datetime.datetime.now()))  
    ### Experimental Distances ###
    ### If pre_scan is specified then it will be used by default ###
    ### If pre_scan is none then the newly scanned experimental sequences will be run ###
    if pre_scan is not None or experimental_fimo == True:
        if verbose == True:
            print('--------------Experimental Distances--------------')        
        run_distance_calculation(verbose=verbose, outdir=outdir, annotation=annotation, sample=sample, window=window, 
                             cpus=cpus, seq_type='experimental', 
                                 pre_scan=pre_scan) 
    if verbose == True:
        print('--------------Calculating Distances Complete--------------')
        print('Stop time: %s' % str(datetime.datetime.now()))   