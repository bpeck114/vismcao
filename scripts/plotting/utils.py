# IMPORT PACKAGES
# For extracting data from MAOS sims for plotting

import os
import glob
import subprocess
import shutil
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from paarti.utils import maos_utils
import readbin

def get_wfe_metrics_over_field(directory='./', seed=1):
    """
    Function to get various wavefront error
    (WFE) metrics from MAOS PSF sims.

    Inputs:
    ----------------------------------------
    - directory (string): path to psf 
    directory, default is current working 
    directory.
    - seed (int): seed with which the 
    simulation was run, default is 1. 

    Outputs:
    ----------------------------------------
    - open_mean_nm (3darray): a 3darray that
    contains WFE metrics for open-loop MAOS
    results averaged over all the PSF 
    evaluation locations.
    - closed_mean_nm (3darray): a 3darray that
    contains WFE metrics for closed-loop 
    MAOS results averaged over all the PSF
    evaluation locations.
    - open_xx_mean_nm ([N, 3]): a 3darray that
    contains WFE metrics for open-loop MAOS
    results for N number of PSF locations, it 
    will return None if only a single PSF
    location.
    - closed_xx_mean_nm ([N, 3]): a 3darray that
    contains WFE metrics for closed-loop MAOS
    results for N number of PSF locations, it 
    will return None if only a single PSF
    location.    
    """

    # Field-averaged results
    results_file = f'{directory}/Res_{seed}.bin'
    results = readbin.readbin(results_file)
    print("Looking in directory:", directory)

    # Open-loop WFE (nm): Piston removed, TT only, Piston+TT removed
    open_mean_nm = np.sqrt(results[0].mean(axis=0)) * 1.0e9

    # Closed-loop WFE (nm): Piston removed, TT only, Piston+TT removed
    clos_mean_nm = np.sqrt(results[2].mean(axis=0)) * 1.0e9

    # Field-dependent results
    # Determine if we have a field-dependent WFE results file in extra/
    results_xx_file = f'{directory}/Resp_{seed}.bin'
    results_xx_file_old = f'{directory}extra/Resp_{seed}.bin'
    if os.path.exists(results_xx_file):
        results_xx = readbin.readbin(results_xx_file)

        open_xx_mean_nm = np.zeros((results_xx[2].shape[0], 3), dtype=float)
        clos_xx_mean_nm = np.zeros((results_xx[3].shape[0], 3), dtype=float)

        # Loop through PSF positions and get RMS WFE in nm
        for xx in range(open_xx_mean_nm.shape[0]):
            
            # Open-loop WFE (nm): Piston removed, TT only, Piston+TT removed
            open_xx_mean_nm[xx] = np.sqrt(results_xx[2][xx].mean(axis=0)) * 1.0e9
            
            # Closed-loop WFE (nm): Piston removed, TT only, Piston+TT removed
            clos_xx_mean_nm[xx] = np.sqrt(results_xx[3][xx].mean(axis=0)) * 1.0e9

    elif os.path.exists(results_xx_file_old):
        results_xx = readbin.readbin(results_xx_file_old)

        open_xx_mean_nm = np.zeros((results_xx[2].shape[0], 3), dtype=float)
        clos_xx_mean_nm = np.zeros((results_xx[3].shape[0], 3), dtype=float)

        # Loop through PSF positions and get RMS WFE in nm
        for xx in range(open_xx_mean_nm.shape[0]):
            
            # Open-loop WFE (nm): Piston removed, TT only, Piston+TT removed
            open_xx_mean_nm[xx] = np.sqrt(results_xx[2][xx].mean(axis=0)) * 1.0e9
            
            # Closed-loop WFE (nm): Piston removed, TT only, Piston+TT removed
            clos_xx_mean_nm[xx] = np.sqrt(results_xx[3][xx].mean(axis=0)) * 1.0e9

    else:
        results_xx = None
        open_xx_mean_nm = None
        clos_xx_mean_nm = None

    return open_mean_nm, clos_mean_nm, open_xx_mean_nm, clos_xx_mean_nm