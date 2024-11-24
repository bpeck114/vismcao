# FILE vismcao_actuator_script.py
# For running MAOS simulations on NERSC, specifically actuator study

# Importing necessary packages
import numpy as np
import subprocess
import time

command = f"maos -o studies/grid/act_study/3000actuatorsgrid -c master_files/A_keck_mcao_lgs_3000.conf -O"
subprocess.run(command, shell=True, text=True)
