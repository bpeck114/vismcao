# FILE vismcao_actuator_script.py
# For running MAOS simulations on NERSC, specifically actuator study

# Importing necessary packages
import numpy as np
import subprocess
import time

# Set parameters for actuator study
output_file = "studies/grid/act_study"
act_count = np.array([3000])
duration = [] # For recording how long each simulation takes

# Describe what is happening
print("Running actuator study:\n")

# Run MAOS simulations for actuator study
for act in act_count:
    command = f"maos -o {output_file}/{act}actuators -c master_files/A_keck_mcao_lgs_{act}.conf plot.all=1 plot.setup=1 -O"
    
    print("---------------------------------------")
    print("SIM:", act)
    print("COMMAND:", command)
    print("MASTER:",f"A_keck_mcao_lgs_{act}.conf") 
    print("LOCATION:", output_file)
    print("---------------------------------------")

    # Run simulation (comment out the following line for testing on non-NERSC machine)
    start_time = time.time()
    subprocess.run(command, shell=True, text=True)
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Convert time to hours, minutes, seconds and milliseconds
    hours = int(elapsed_time // 3600)
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)
    milliseconds = int((elapsed_time - int(elapsed_time)) * 1000)

    duration.append({
        "sim": act,
        "time_hours": hours,
        "time_minutes" : minutes,
        "time_seconds" : seconds,
        "time_milliseconds": milliseconds
    })

    print("\n")

print("---------------------------------------")
for dur in duration:
    print("SIM ({}) TIME:". format(dur["sim"]), "{}:{}:{}.{}".format(dur["time_hours"], dur["time_minutes"], dur["time_seconds"], dur["time_milliseconds"]))
print("---------------------------------------")
