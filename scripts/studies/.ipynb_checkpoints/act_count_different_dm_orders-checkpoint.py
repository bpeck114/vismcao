# FILE act_count_different_dm_orders.py
# For running "actuator count" KOLA simulations on NERSC

# IMPORT NECESSARY PACKAGES
import numpy as np
import subprocess
import time 
import os
from datetime import datetime
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run actuator count simulations.")
parser.add_argument("--test", action="store_true", help="Run in test mode without executing simulations.")
args = parser.parse_args()

# Set parameters for actuator study
config_folder = "master_files/dm_orders/8mag_lgs/"
default_output_folder = "studies/dm_orders/8mag_lgs/"
test_folder = "studies/dm_orders/8mag_lgs/test/"
dm_orders = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])

# Choose output folder based on test mode
output_folder = test_folder if args.test else default_output_folder
os.makedirs(output_folder, exist_ok=True)

# Initialize log file for this batch
batch_id = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = os.path.join(output_folder, f"act_count_log_{batch_id}.txt")

# Initalize duration list
duration = []

print(f"Beginning dm order study... {'(TEST MODE)' if args.test else ''}\n")

with open(log_file, "w") as log:
    log.write(f"Actuator Count Study Log - Batch ID: {batch_id}\n")
    log.write(f"Output Folder: {output_folder}\n")
    log.write(f"Test Mode: {'Enabled' if args.test else 'Disabled'}\n")
    log.write("---------------------------------------\n\n")

    for dm in dm_orders:
        output_path = os.path.join(output_folder, f"{dm}order")
        config_file = os.path.join(config_folder, f"kola_{dm}order.conf")
        command = f"maos -o {output_path} -c {config_file} -O"

        # Check if the configuration file exists
        if os.path.isfile(config_file):
            print(f"Found configuration file: {config_file}")
            log.write(f"Found configuration file: {config_file}")
        else: 
            error_message = f"Configuration file not found: {config_file}"
            print(error_message)
            log.write(f"ERROR: {error_message}\n\n")
            continue

        log.write(f"SIMULATION: {dm}\n")
        log.write(f"COMMAND: {command}\n")

        if args.test:
            print(f"(TEST MODE) Skipping execution of simulation {dm} order.")
            log.write(f"(TEST MODE) Skipping execution of simulation {dm} order.")
            
            print(command)
            log.write(command)
            continue
            
        # Run the simulation and measure time
        start_time = time.perf_counter()
        try: 
            subprocess.run(command, shell=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            log.write(f"Error: {e}\n\n")
            continue

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        # Convert elapsed time into hours, minutes, seconds, milliseconds
        hours, rem = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(rem, 60)
        milliseconds = int((elapsed_time - int(elapsed_time)) * 1000)

        duration.append({
            "sim": act,
            "time_hours": int(hours),
            "time_minutes": int(minutes),
            "time_seconds": int(seconds),
            "time_milliseconds": milliseconds
        })

        log.write(f"Elapsed Time: {int(hours)}:{int(minutes)}:{int(seconds)}.{milliseconds}\n\n")
        print("\n")
        
    print("---------------------------------------")
    log.write("---------------------------------------\n")
    for dur in duration:
        print("SIM ({}) TIME:".format(dur["sim"]), "{}:{}:{}.{}".format(dur["time_hours"], dur["time_minutes"], dur["time_seconds"], dur["time_milliseconds"]))
        log.write("SIM ({}) TIME: {}:{}:{}.{}\n".format(dur["sim"], dur["time_hours"], dur["time_minutes"], dur["time_seconds"], dur["time_milliseconds"]))
    print("---------------------------------------")
    log.write("---------------------------------------\n")

print(f"Simulation log saved to: {log_file}")
