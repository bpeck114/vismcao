# FILE script_master.py
# For running multiple vismcao simulations

# Importing necessary packages
import numpy as np
import math
import subprocess
import time
import os
from itertools import product
from paarti.utils import maos_utils

import script_vismcao

output_folder           = 'error_budget'           # Name of output folder for simulation study and type of loop-over
loop_over_values        = [2000]
master_vismcao_file     = 'A_vismcao.conf'         # Name of selected master vismcao file for simulations
n_actuator_asm          = 4000                     # Number of actuators on adaptive secondary mirror
h_asm                   = -100                     # Conjugate height of adaptive secondary mirror (meters)
n_downstream_dms        = 2                        # Number of deformable mirrors (1 = SCAO, >2 = MCAO)
side_downstream_dms     = [.168, .168]             # Distance between actuators on downstream deformable mirrors
h_downstream_dms        = [6000, 10000]            # Conjugate height of deformable mirrors (meters)
n_sodium                = 8                        # Number of sodium laser guide stars
w_sodium                = 20                       # Wattage of all sodium laser guide stars
r_sodium                = 30                       # Radius of sodium laser guide stars
b_sodium                = 0.1                      # Sodium laser guide star background (powfs.bkgrnd)
on_axis_1_sodium_beacon = True                     # Overrides the default behavior of single NaLGS is on-axis
n_tiptilt               = 3                        # Number of tip-tilt stars
m_tiptilt               = 8                        # Magnitude of all tip-tilt stars
r_tiptilt               = 30                       # Radius of all tip-tilt stars
field_of_view           = 30                       # Fitted field of view (should match r_sodium)
tomo_iterations         = 40                       # Number of iterations for conjugate gradient (CG) algorithm
integration_time        = 1/1500                   # Integration time for simulation 
simulation_steps        = 1600                     # Number of steps in one MAOS simulation
run_simulation          = False                     # True runs simulations, False prints MAOS command

script_vismcao.main(output_folder=output_folder, loop_over_values=loop_over_values, master_vismcao_file=master_vismcao_file, n_actuator_asm=n_actuator_asm, h_asm=h_asm, n_downstream_dms=n_downstream_dms, side_downstream_dms=side_downstream_dms, h_downstream_dms=h_downstream_dms, n_sodium=n_sodium, w_sodium=w_sodium, r_sodium=r_sodium, b_sodium=b_sodium, n_tiptilt=n_tiptilt, m_tiptilt=m_tiptilt, r_tiptilt=r_tiptilt, field_of_view=field_of_view, tomo_iterations=tomo_iterations, integration_time=integration_time, simulation_steps=simulation_steps, run_simulation=run_simulation)




'''
#2NaLGS + 0RLGS
output_file = "height_study/infW_2S_0R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 2                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [0]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 20                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)


#0NaLGS + 10km 1RLGS
output_file = "height_study/infW_0S_10km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 10                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 20km 1RLGS
output_file = "height_study/infW_0S_20km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 20                                 # 20km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 30km 1RLGS
output_file = "height_study/infW_0S_30km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 30                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 40km 1RLGS
output_file = "height_study/infW_0S_40km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 40                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 50km 1RLGS
output_file = "height_study/infW_0S_50km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 50                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 60km 1RLGS
output_file = "height_study/infW_0S_60km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 60                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 70km 1RLGS
output_file = "height_study/infW_0S_70km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 70                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 80km 1RLGS
output_file = "height_study/infW_0S_80km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 80                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#0NaLGS + 90km 1RLGS
output_file = "height_study/infW_0S_90km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 0                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 90                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 0RLGS
output_file = "height_study/infW_1S_0R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [0]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 20                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 10km 1RLGS
output_file = "height_study/infW_1S_10km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 10                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 20km 1RLGS
output_file = "height_study/infW_1S_20km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 20                                 # 20km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 30km 1RLGS
output_file = "height_study/infW_1S_30km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 30                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 40km 1RLGS
output_file = "height_study/infW_1S_40km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 40                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 50km 1RLGS
output_file = "height_study/infW_1S_50km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 50                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 60km 1RLGS
output_file = "height_study/infW_1S_60km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 60                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 70km 1RLGS
output_file = "height_study/infW_1S_70km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 70                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 80km 1RLGS
output_file = "height_study/infW_1S_80km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 80                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)

#1NaLGS + 90km 1RLGS
output_file = "height_study/infW_1S_90km_R"     # Height Study, 90km RLGS should collapse into NaLGS if MAOS is working properly with "infinite brightness"
master_scao_file = "A_hybrid.conf"
master_mcao_file = "A_mcao_hybrid.conf" 
n_sodium = 1                                    # 0 NaLGS
r_sodium = 0                                    # On-axis NaLGS (if present)
w_sodium = "inf"                                # Infinite Brightness
split_sodium_beacon = False
on_axis_1_sodium_beacon = True                  # On-axis NaLGS (if present)
n_rayleigh = [1]                                # 1 RLGS
r_rayleigh = 0                                  # On-axis
w_rayleigh = "inf"                              # Infinite Brightness
h_rayleigh = 90                                 # 10km Height RLGS
on_axis_1_rayleigh_beacon = True 
m_tiptilt_truth = 8
integration_time = 1/1500 
tomo_iterations = 40
zenith = 30 
single_conjugate = True                         # Single-Conjugate
run_simulation = True

script_hybrid.main(output_file=output_file, master_scao_file=master_scao_file, master_mcao_file=master_mcao_file, n_sodium=n_sodium, r_sodium=r_sodium, w_sodium=w_sodium, split_sodium_beacon=split_sodium_beacon, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon, n_rayleigh=n_rayleigh, r_rayleigh=r_rayleigh, w_rayleigh=w_rayleigh, h_rayleigh=h_rayleigh, on_axis_1_rayleigh_beacon=on_axis_1_rayleigh_beacon, m_tiptilt_truth=m_tiptilt_truth, integration_time=integration_time, tomo_iterations=tomo_iterations, zenith=zenith, single_conjugate=single_conjugate, run_simulation=run_simulation)
'''