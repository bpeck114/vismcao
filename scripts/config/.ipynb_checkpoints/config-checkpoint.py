# FILE config.py

# IMPORT NECESSARY PACKAGES
import numpy as np
import math
from paarti.utils import maos_utils

def calc_actuator_separation(act_count=[4000], sigfigs=3, verbose=True):
    """
    Calculates the separation between actuators (in m) on Keck's 
    primary mirror. 

    Inputs:
    ------------------------------------
    - act_count (array): an array containing actuator counts
    - sigfigs (int): number of significant figures for rounding,
    default is 3. 
    - verbose (bool): if True, prints the results

    Outputs:
    ------------------------------------
    - sides (list) : calculated separation values (dm.dx) in m
    - actuators (list): associated actuator counts    
    """

    # Calculate the side length per suberapture based on actuator count
    sides = [np.round(11 / (2 *((actuator_count/np.pi)**0.5)), sigfigs)
             for actuator_count in act_count]
    if verbose:
        print('--------------------')
        print('Distance between actuators relative to primary mirror (m), dm.dx:')
        print('--------------------')
        for side_value, actuator_value in zip(sides, act_count):
            print(f'{actuator_value:11.0f} {side_value:11.3f}')
    
    return sides, act_count

def calc_act_flux_params(act_count=[4000], lgs_mag=8, tt_mag=8, LGS_WFS='LGSWFS-OCAM2K', LGS_THROUGHPUT=0.36*0.88, LGS_PIXEL_SIZE=1, LGS_THETA_BETA=1.5 *(math.pi/180)/(60*60), LGS_BAND='R', LGS_SIGMA_E=0.5, LGS_PIXPERSA=25, LGS_INTEGRATION=1/1500, LGS_BKGRND=0.1, TT_WFS='TRICK-H', TT_INTEGRATION=1/1500, LBWFS_WFS='LBWFS', LBWFS_INTEGRATION=1/1500):
    """
    Calculates flux parameters for MAOS
    based on actuator count and LGS, TT
    and LBWFS magnitudes.
    Inputs:
    ------------------------------------
    - act_count (array): an array containing actuator counts

    Outputs:
    ------------------------------------
    - siglev (list): powfs.siglev values for LGS, TT and LWBFS
    - bkgrnd (list): powfs.bkgrnd values for LGS, TT and LWBFS
    - nearecon (list): powfs.nearecon values for LGS, TT and LWBFS
    """

    # Generates actuator separation
    sides, actuators = calc_actuator_separation(act_count, verbose=False)

    lgs_flux_values = [maos_utils.keck_nea_photons_any_config(
        wfs = LGS_WFS,
        side=side,
        throughput=LGS_THROUGHPUT,
        ps=LGS_PIXEL_SIZE,
        theta_beta=LGS_THETA_BETA,
        band=LGS_BAND,
        sigma_e=LGS_SIGMA_E,
        pix_per_ap=LGS_PIXPERSA,
        time=LGS_INTEGRATION,
        m=lgs_mag
    ) for side in sides]
    lgs_siglev, lgs_bkgrnd, lgs_nearecon = zip(*[(round(fv[2], 3), LGS_BKGRND, round(fv[1], 3)) for fv in lgs_flux_values])
    # Calculate TT flux parameters
    tt_flux_values = [maos_utils.keck_nea_photons(
        m=tt_mag, wfs=TT_WFS, wfs_int_time=TT_INTEGRATION
    ) for _ in sides]
    tt_siglev, tt_bkgrnd, tt_nearecon = zip(*[(round(fv[2], 3), round(fv[3], 3), round(fv[1], 3)) for fv in tt_flux_values])

    # Calculate LBWFS flux parameters
    truth_flux_values = [maos_utils.keck_nea_photons(
        m=tt_mag, wfs=LBWFS_WFS, wfs_int_time=LBWFS_INTEGRATION
    ) for _ in sides]
    truth_siglev, truth_bkgrnd, truth_nearecon = zip(*[(round(fv[2], 3), round(fv[3], 3), round(fv[1], 3)) for fv in truth_flux_values])

    print('--------------------')
    print('VisMCAO Magnitude-to-Flux Parameters:')
    print('--------------------\n')

    for i, (actuator, side) in enumerate(zip(actuators, sides)):
        print(f'####\n#{lgs_mag}mag LGS ({tt_mag}mag TT)\n####')
        print('#Actuator Count:', actuator)
        print('#dm.dx = [', side, '.168 .168 ]')
        print('#powfs.siglev = [', lgs_siglev[i], tt_siglev[i], truth_siglev[i], ']')
        print('#powfs.bkgrnd = [', lgs_bkgrnd[i], tt_bkgrnd[i], truth_bkgrnd[i], ']')
        print('#powfs.nearecon = [', lgs_nearecon[i], tt_nearecon[i], truth_nearecon[i], ']')
        print('')

    return {
        "siglev": {"LGS": list(lgs_siglev), "TT": list(tt_siglev), "LBWFS": list(truth_siglev)},
        "bkgrnd": {"LGS": list(lgs_bkgrnd), "TT": list(tt_bkgrnd), "LBWFS": list(truth_bkgrnd)},
        "nearecon": {"LGS": list(lgs_nearecon), "TT": list(tt_nearecon), "LBWFS": list(truth_nearecon)}
    }

    


    






    