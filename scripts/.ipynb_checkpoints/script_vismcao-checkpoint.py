# FILE script_vismcao.py
# For running visible-light multi-conjugate adaptive optics simulations with MAOS
# MAOS = Multi-Threaded Adaptive Optics Simulator



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%% IMPORTING NECESSARY PACKAGES  %%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import numpy as np
import math
import subprocess
import time
import os
from itertools import product
from datetime import datetime

# Import internal packages (needs additional installation)
from paarti.utils import maos_utils



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%% MAOS Parameters %%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Name of study (output folder), loop-over values' and master file
# Choose 'asm_actuator_study', 'sodium_wattage_study', 'sodium_count_study', 'sodium_radius_study', 'dm_height_study', 'error_budget', 'other' for output folder
output_folder           = 'actuator_study'            # Name of output folder for simulation study and type of loop-over
loop_over_values        = [2000, 3000, 4000, 5000]
master_vismcao_file     = 'A_keck_mcao_lgs_4000.conf' # Name of selected master vismcao file for simulations

# Adaptive Secondary Mirror
n_actuator_asm          = 4000                        # Number of actuators on adaptive secondary mirror
h_asm                   = -100                        # Conjugate height of adaptive secondary mirror (meters)

# Deformable Mirrors
n_downstream_dms        = 2                           # Number of deformable mirrors (1 = SCAO, >2 = MCAO)
side_downstream_dms     = [.168, .168]                # Distance between actuators on downstream deformable mirrors
h_downstream_dms        = [6000, 10000]               # Conjugate height of deformable mirrors (meters)

# Sodium Laser Guide Stars
n_sodium                = 8                           # Number of sodium laser guide stars
w_sodium                = 60                          # Wattage of all sodium laser guide stars
r_sodium                = 30                          # Radius of sodium laser guide stars
b_sodium                = 0.1                         # Sodium laser guide star background (powfs.bkgrnd)
on_axis_1_sodium_beacon = True                        # Overrides the default behavior of single NaLGS is on-axis (True = 1 NaLGS on-axis)

# Tip-Tilt Natural Guide Stars
n_tiptilt               = 3                           # Number of tip-tilt stars
m_tiptilt               = 8                           # Magnitude of all tip-tilt stars
r_tiptilt               = 30                          # Radius of all tip-tilt stars

# PSF Evaluation + Other
field_of_view           = 30                          # Fitted field of view (should match r_sodium)
tomo_iterations         = 40                          # Number of iterations for conjugate gradient (CG) algorithm 
integration_time        = 1500                        # Integration time for simulation (sim.dt/sim.dtref, set in master file, NOT here)

# Simulation 
simulation_steps       = 1600                        # Number of steps in one MAOS simulation
run_simulation          = False                       # True runs simulations, False prints MAOS command



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%% Loop Parameters %%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#  %%%%%%%% INITIAL FUNCTIONS  %%%%%%%%
def calc_side(actuator_count, sigfigs=3):
    # Calculates the distance (in meters) between actuators on primary mirror
    # Known as dm.dx in MAOS configuration files
    # Needed for keck_nea_photons_any_config in PAARTI
    # Assume 11 meter diameter telescope (Keck parameters)

    # Calculate the side length per suberapture based on actuator count
    side_value = np.round(11 / (2 * ((actuator_count/np.pi)**0.5)), sigfigs)

    return side_value

def calculate_circle_coordinates(n, radius):
    # Calculates (x, y) coordinate points evenly around the circumference of a circle
    # First point always begins on x-axis
    # Meant for calculating guide-star asterisms 
    
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x_coords = np.round(radius * np.cos(angles), 2)
    y_coords = np.round(radius * np.sin(angles), 2)
    return x_coords, y_coords

def calculate_magnitude_from_wattage(sodium_wattage): 
    # Sodium LGS Photometry
    # Calculates the magnitude of a sodium laser guide star given its wattage/power
    # Assumes a 20 W sodium sodium laser guide star is 8.1 magnitude

    assumed_sodium_wattage = 20 # watts
    assumed_sodium_magnitude = 8.1 # magnitude

    magnitude_difference = -2.5 * math.log10(sodium_wattage/assumed_sodium_wattage)
    sodium_magnitude = assumed_sodium_magnitude + magnitude_difference

    return sodium_magnitude

def calculate_nearecon_from_wattage(wattage, height):
    # Rayleigh LGS Photometry
    # Photometry constants from Robo-AO
    
    assumed_height = 10 # kilometers
    assumed_power = 20 # watts
    assumed_nearecon = 104 # milliarcseconds
    scale_height = 8 # kilometers

    # Calculate the multiplier from the wattage
    multiplier = wattage / assumed_wattage
    
    # Calculate the base ten km near econ
    base_ten_km_nearecon = assumed_nearecon * (multiplier ** (-1/2))
    
    # Calculate how noise-equivalent angle increases with altitude
    exponent = np.exp(- ( (assumed_height - height) / scale_height ))
    kilometer_function = (((assumed_height / height) * exponent) ** (1/2))
    nearecon = base_ten_km_nearecon * kilometer_function

    return nearecon

def create_dm_dx_array(n_dms, side_dms, side_asm):
    # Combines distance between actuators for ASM with distance between actautors for downstream deformable mirrors
    # dm.dx value in MAOS configuration files

    if n_dms == len(side_dms): # Checks to make sure that given values match
        return ' '.join(map(str, [side_asm] + side_dms)) # Returns combined array that MAOS accepts
    elif n_dms != len(side_dms):
        raise ValueError("Number of deformable mirrors (n_downstream_dms) does not match number of distances between actuators (side_dms).")

def create_dm_ht_array(n_dms, height_dms, height_asm):
    # Combines the conjugate height of the ASM with the conjugate height of downstream deformable mirrors
    # Almost identical to create_dm_dx_array() above
    
    if n_dms == len(height_dms):
        return ' '.join(map(str, [height_asm] + height_dms))
    elif n_dms != len(side_dms):
        raise ValueError("Number of deformable mirrors (n_downstream_dms) does not match number of heights of deformbale mirrors (height_dms).")


def photometry(n_sodium, n_tiptilt, b_sodium, w_sodium, m_tiptilt, asm_side, integration_time):
    # Takes in number, wattage/magnitude, distance between actuators and integration time
    # Determine signal level, background (unless given), signal-to-noise ratio and noise-equivalent angle
    # For MCAO Keck simulations, specifically, LGS/NGS/Truth WFS
    # DOES NOT HANDLE RAYLEIGH BEACONS!

    if w_sodium == 'inf':
        # Sets really high photon signal level for wavefront sensor 
        sodium_siglev = 1000000 # 1 million photons (might be too much, hasn't been tested...) 
        sodium_nearecon = 1 # Really low noise-equivalent angle (book-keeping essentially)
    else: 
        m_sodium = calculate_magnitude_from_wattage(sodium_wattage=w_sodium)
        print(f"Photometry, m_sodium = {m_sodium}")
        # Setting photometric values for sodium beacon 
        _, nalgs_nearecon_all, nalgs_siglev_full_all, nalgs_bkgrnd_all = maos_utils.keck_nea_photons_any_config(wfs='LGSWFS',
                                                                                                                side = asm_side, # dm.dx value 
                                                                                                                throughput = 0.36 * 0.88,
                                                                                                                ps = 3.0,
                                                                                                                theta_beta = 1.5*(math.pi/180)/(60.0*60.0),
                                                                                                                band = "R",
                                                                                                                sigma_e = 0.5,
                                                                                                                pix_per_ap = 4,
                                                                                                                time = integration_time,
                                                                                                                m = m_sodium)
        # Round the values that PAARTI gives for sodium beacons
        sodium_siglev = np.round(nalgs_siglev_full_all, 2) # powfs.siglev does not affect performance if tomo.precond = 0 
        sodium_nearecon = np.round(nalgs_nearecon_all, 2)
        print(f"Photometry, sodium_siglev = {sodium_siglev}")

    sodium_bkgrnd = b_sodium
            
    # Setting photometric values for tip-tilt stars and truth (low-bandwidth) wavefront sensor
    tt_snr_all, tt_nearecon_all, tt_siglev_all, tt_bkgrnd_all = maos_utils.keck_nea_photons(m=m_tiptilt, wfs='STRAP', wfs_int_time=integration_time)
    truth_snr_all, truth_nearecon_all, truth_siglev_all, truth_bkgrnd_all = maos_utils.keck_nea_photons(m=m_tiptilt, wfs='LBWFS', wfs_int_time=integration_time)
        
    # Round the values that PAARTI gives for tip-tilt and truth
    tt_siglev = np.round(tt_siglev_all, 2)
    print(f"Photometry, tt_siglev = {tt_siglev}")
    tt_bkgrnd = np.round(tt_bkgrnd_all, 2)
    tt_nearecon = np.round(tt_nearecon_all, 2)
        
    truth_siglev = np.round(truth_siglev_all, 2)
    truth_bkgrnd = np.round(truth_bkgrnd_all, 2)
    truth_nearecon = np.round(truth_nearecon_all, 2)

    return sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon

def asterisms(n_sodium, r_sodium, n_tiptilt, r_tiptilt, on_axis_1_sodium_beacon):
    # Takes in number and radius for guide stars
    # Determines postion on sky, evenly spaced around circumference with given radius
    # Defaults to one beacon is placed on-axis for LGS
    if on_axis_1_sodium_beacon: 
        # Places n = 1 sodium beacon on-axis (wavefront-sensor + laser-launch telescope)
            if n_sodium == 1:
                wfs_sodium_x = 0
                wfs_sodium_y = 0
                llt_sodium_x = 0
                llt_sodium_y = 0
            else: 
                # Calculate positions of sodium wave-front sensors and laser launch telescopes 
                wfs_sodium_xx, wfs_sodium_yy = calculate_circle_coordinates(n_sodium, r_sodium)
                llt_sodium_xx, llt_sodium_yy = calculate_circle_coordinates(n_sodium, radius=1) # Multiplied by 6.5 m later
                
                # Set the outputs in a format that MAOS accepts (no commas) 
                wfs_sodium_x = ' '.join(map(str, wfs_sodium_xx))
                wfs_sodium_y = ' '.join(map(str, wfs_sodium_yy))
                llt_sodium_x = ' '.join(map(str, llt_sodium_xx))
                llt_sodium_y = ' '.join(map(str, llt_sodium_yy))
    else: 
        # Calculate positions of sodium wave-front sensors and laser launch telescopes
        wfs_sodium_xx, wfs_sodium_yy = calculate_circle_coordinates(n_sodium, r_sodium)
        llt_sodium_xx, llt_sodium_yy = calculate_circle_coordinates(n_sodium, radius=1) # Multiplied by 6.5 m later
            
        # Set the outputs in a format that MAOS accepts (no commas) 
        wfs_sodium_x = ' '.join(map(str, wfs_sodium_xx))
        wfs_sodium_y = ' '.join(map(str, wfs_sodium_yy))
        llt_sodium_x = ' '.join(map(str, llt_sodium_xx))
        llt_sodium_y = ' '.join(map(str, llt_sodium_yy))

    # Caluclates positions of tip-tilt wavefront sensors
    wfs_tiptilt_xx, wfs_tiptilt_yy = calculate_circle_coordinates(n_tiptilt, r_tiptilt)
            
    # Set the outputs in a format that MAOS accepts (no commas) 
    wfs_tiptilt_x = ' '.join(map(str, wfs_tiptilt_xx))
    wfs_tiptilt_y = ' '.join(map(str, wfs_tiptilt_yy))

    return wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y

def format_time(seconds):
    # Converts a duration of a simulation from seconds to a formatted string
    
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    milliseconds = int((seconds - int(seconds)) * 1000)
    
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}.{milliseconds:03}"



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%% MAIN SCRIPT BEGINS HERE %%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def main(output_folder, loop_over_values, master_vismcao_file, n_actuator_asm, h_asm, n_downstream_dms, side_downstream_dms, h_downstream_dms, n_sodium, w_sodium, r_sodium, b_sodium, n_tiptilt, m_tiptilt, r_tiptilt, field_of_view, tomo_iterations, integration_time, simulation_steps, run_simulation):
    # For running visible-light multi-conjugate laser guide star adaptive optics Keck simulations
    
    initial_cwd = os.getcwd()
    print(f"Simulation starting in {initial_cwd}.")
    print(f"Running '{output_folder}' vismcao study with '{master_vismcao_file}'.\n")

    # Choose which sub-directories to create based on the `loop_over` input
    if output_folder == 'asm_actuator_study':

        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
        
        with open(time_log_file, 'a') as time_file:
            for actuator in loop_over_values:
                # Actuator Study
                side_actuator_asm = calc_side(actuator_count=actuator, sigfigs=3)
                dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

                # Wattage
                sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=w_sodium, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

                # Radius
                wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=n_sodium, r_sodium=r_sodium, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)

                command = f"maos -o {output_folder}/{actuator}actuators -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}'" 

                print(command + "\n")

                # Record the duration of the simulation
                start_time = time.time()
                if run_simulation:
                    subprocess.run(command, shell=True, text=True)
                    
                    end_cwd = os.getcwd()
                    print(f"Finished simulation in {end_cwd}.")
    
                end_time = time.time()
                duration = end_time - start_time

                # Get the current date and time when the simulation ends
                end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Print duration of simulation to a text file
                time_file.write(f"As of {end_datetime}: Simulation for {actuator} actuators completed in {format_time(duration)}.\n")
                print(f"As of {end_datetime}: Simulation for {actuator} actuators completed in {format_time(duration)}.\n")
            
            total_end_time = time.time()
            total_duration = total_end_time - total_start_time
            total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")

            if run_simulation:
                os.chdir("..")
                os.chdir("..")

    elif output_folder == 'sodium_wattage_study':
        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
        
        with open(time_log_file, 'a') as time_file:
            for sodium_wattage in loop_over_values:
                # Wattage Study
                sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=sodium_wattage, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

                # Actuator 
                side_actuator_asm = calc_side(actuator_count=n_actuator_asm, sigfigs=3)
                dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

                # Radius
                wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=n_sodium, r_sodium=r_sodium, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)
            
                command = f"maos -o {output_folder}/{sodium_wattage}sodium_wattage -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}'" 

                print(command + "\n")
             
                # Record the duration of the simulation
                start_time = time.time()
                if run_simulation:
                    subprocess.run(command, shell=True, text=True)
                                        
                    end_cwd = os.getcwd()
                    print(f"Finished simulation in {end_cwd}.")
    
                end_time = time.time()
                duration = end_time - start_time

                # Get the current date and time when the simulation ends
                end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Print duration of simulation to a text file
                time_file.write(f"As of {end_datetime}: Simulation for {sodium_wattage} W sodium LGS completed in {format_time(duration)}.\n")
                print(f"As of {end_datetime}: Simulation for {sodium_wattage} W sodium LGS completed in {format_time(duration)}.\n")
            
            total_end_time = time.time()
            total_duration = total_end_time - total_start_time
            total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")
            
            if run_simulation:
                os.chdir("..")
                os.chdir("..")

    elif output_folder == 'r_sodium':
        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
        
        with open(time_log_file, 'a') as time_file:
            for sodium_radius in loop_over_values:
                # Radius Study
                wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=n_sodium, r_sodium=sodium_radius, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)
            
                # Actuator 
                side_actuator_asm = calc_side(actuator_count=n_actuator_asm, sigfigs=3)
                dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

                # Wattage
                sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=w_sodium, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

                command = f"maos -o {output_folder}/{sodium_radius}sodium_radius -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}'" 

                print(command + "\n")
                             
                # Record the duration of the simulation
                start_time = time.time()
                if run_simulation:
                    subprocess.run(command, shell=True, text=True)
                                                            
                    end_cwd = os.getcwd()
                    print(f"Finished simulation in {end_cwd}.")
    
                end_time = time.time()
                duration = end_time - start_time

                # Get the current date and time when the simulation ends
                end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Print duration of simulation to a text file
                time_file.write(f"As of {end_datetime}: Simulation for {sodium_radius} arcsecond radius of sodium LGS completed in {format_time(duration)}.\n")
                print(f"As of {end_datetime}: Simulation for {sodium_radius} arcsecond radius of sodium LGS completed in {format_time(duration)}.\n")
            
            total_end_time = time.time()
            total_duration = total_end_time - total_start_time
            total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")
            
            if run_simulation:
                os.chdir("..")
                os.chdir("..")

    elif output_folder == 'sodium_count_study':
        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
        
        with open(time_log_file, 'a') as time_file:
            for sodium_number in loop_over_values:
                # Number Study
                wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=sodium_number, r_sodium=r_sodium, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)
            
                # Actuator 
                side_actuator_asm = calc_side(actuator_count=n_actuator_asm, sigfigs=3)
                dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

                # Wattage
                sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=w_sodium, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

                command = f"maos -o {output_folder}/{sodium_number}sodium_number -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}'" 

                print(command + "\n")
                             
                # Record the duration of the simulation
                start_time = time.time()
                if run_simulation:
                    subprocess.run(command, shell=True, text=True)
                    
                    end_cwd = os.getcwd()
                    print(f"Finished simulation in {end_cwd}.")
    
                end_time = time.time()
                duration = end_time - start_time

                # Get the current date and time when the simulation ends
                end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Print duration of simulation to a text file
                time_file.write(f"As of {end_datetime}: Simulation for {sodium_number} sodium LGS completed in {format_time(duration)}.\n")
                print(f"As of {end_datetime}: Simulation for {sodium_number} sodium LGS completed in {format_time(duration)}.\n")
            
            total_end_time = time.time()
            total_duration = total_end_time - total_start_time
            total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")
            
            if run_simulation:
                os.chdir("..")
                os.chdir("..")
    
    elif output_folder == 'dm_two_height_study':
        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
        
        with open(time_log_file, 'a') as time_file:
            for dm_two_height in loop_over_values:
                side_downstream_dms[0] = dm_two_height

                # DM2 Height Study 
                side_actuator_asm = calc_side(actuator_count=n_actuator_asm, sigfigs=3)
                dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

                # Wattage
                sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=w_sodium, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

                # Radius
                wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=n_sodium, r_sodium=r_sodium, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)

                command = f"maos -o {output_folder}/{dm_two_height}dm_two_height -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}'"    

                print(command + "\n")
                                             
                # Record the duration of the simulation
                start_time = time.time()
                if run_simulation:
                    subprocess.run(command, shell=True, text=True)
                                        
                    end_cwd = os.getcwd()
                    print(f"Finished simulation in {end_cwd}.")

                end_time = time.time()
                duration = end_time - start_time

                # Get the current date and time when the simulation ends
                end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Print duration of simulation to a text file
                time_file.write(f"As of {end_datetime}: Simulation for DM2 = {dm_two_height} meters completed in {format_time(duration)}.\n")
                print(f"As of {end_datetime}: Simulation for DM2 = {dm_two_height} meters sodium LGS completed in {format_time(duration)}.\n")
            
            total_end_time = time.time()
            total_duration = total_end_time - total_start_time
            total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")
            
            if run_simulation:
                os.chdir("..")
                os.chdir("..")

    elif output_folder == 'dm_three_height_study':
        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
        
        with open(time_log_file, 'a') as time_file:
            for dm_three_height in loop_over_values:

                side_downstream_dms[1] = dm_two_height

                # DM2 Height Study 
                side_actuator_asm = calc_side(actuator_count=n_actuator_asm, sigfigs=3)
                dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

                # Wattage
                sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=w_sodium, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

                # Radius
                wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=n_sodium, r_sodium=r_sodium, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)

                command = f"maos -o {output_folder}/{dm_three_height}dm_three_height -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}'"   

                print(command + "\n")
                                                             
                # Record the duration of the simulation
                start_time = time.time()
                if run_simulation:
                    subprocess.run(command, shell=True, text=True)
                                        
                    end_cwd = os.getcwd()
                    print(f"Finished {output_folder} simulation in {end_cwd}.")

                end_time = time.time()
                duration = end_time - start_time

                # Get the current date and time when the simulation ends
                end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                # Print duration of simulation to a text file
                time_file.write(f"As of {end_datetime}: Simulation for DM3 = {dm_two_height} meters completed in {format_time(duration)}.\n")
                print(f"As of {end_datetime}: Simulation for DM3 = {dm_two_height} meters sodium LGS completed in {format_time(duration)}.\n")
            
            total_end_time = time.time()
            total_duration = total_end_time - total_start_time
            total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")
            
            if run_simulation:
                os.chdir("..")
                os.chdir("..")

    elif output_folder == 'error_budget':   
        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
        
        with open(time_log_file, 'a') as time_file:
            # Actuator 
            side_actuator_asm = calc_side(actuator_count=n_actuator_asm, sigfigs=3)
            dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

            # Wattage
            sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=w_sodium, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

            # Radius
            wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=n_sodium, r_sodium=r_sodium, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)

            command = f"maos -o {output_folder}/step1 -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}' 'sim.closeloop=0' 'atm.frozenflow=1' 'sim.idealfit=1' 'sim.eptwfs=0' 'tomo.alg=0' 'fit.alg=0' 'surf=[]' 'sim.wspsd=[]'"

            print(command + "\n")
            '''                                                             
            # Record the duration of the simulation
            start_time = time.time()
            if run_simulation:
                subprocess.run(command, shell=True, text=True)
                                    
                end_cwd = os.getcwd()
                print(f"Finished '{output_folder}/step1' simulation in {end_cwd}.")

            end_time = time.time()
            duration = end_time - start_time

            # Get the current date and time when the simulation ends
            end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Print duration of simulation to a text file
            time_file.write(f"As of {end_datetime}: Simulation for {output_folder}/step1 completed in {format_time(duration)}.\n")
            print(f"As of {end_datetime}: Simulation for {output_folder}/step1 completed in {format_time(duration)}.\n")
            '''      
            end_cwd = os.getcwd()
            print(f"Starting '{output_folder}/step2' simulation in {end_cwd}.")

            command = f"maos -o {output_folder}/step2 -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}' 'sim.closeloop=0' 'atm.frozenflow=1' 'powfs.phystep=-1' 'powfs.noisy=0' 'sim.eptwfs=0' 'tomo.alg=0' 'fit.alg=0'"

            print(command + "\n")
                                                                                     
            # Record the duration of the simulation
            start_time = time.time()
            if run_simulation:
                subprocess.run(command, shell=True, text=True)
                                                    
                end_cwd = os.getcwd()
                print(f"Finished '{output_folder}/step2' simulation in {end_cwd}.")

            end_time = time.time()
            duration = end_time - start_time

            # Get the current date and time when the simulation ends
            end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Print duration of simulation to a text file
            time_file.write(f"As of {end_datetime}: Simulation for {output_folder}/step2 completed in {format_time(duration)}.\n")
            print(f"As of {end_datetime}: Simulation for {output_folder}/step2 completed in {format_time(duration)}.\n")
                   
            end_cwd = os.getcwd()
            print(f"Starting '{output_folder}/step3' simulation in {end_cwd}.")

            command = f"maos -o {output_folder}/step3 -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}' 'powfs.phystep=-1' 'powfs.noisy=0'"

            print(command + "\n")
                                                                         
            # Record the duration of the simulation
            start_time = time.time()
            if run_simulation:
                subprocess.run(command, shell=True, text=True)
                                                    
                end_cwd = os.getcwd()
                print(f"Finished '{output_folder}/step3' simulation in {end_cwd}.")

            end_time = time.time()
            duration = end_time - start_time

            # Get the current date and time when the simulation ends
            end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Print duration of simulation to a text file
            time_file.write(f"As of {end_datetime}: Simulation for {output_folder}/step3 completed in {format_time(duration)}.\n")
            print(f"As of {end_datetime}: Simulation for {output_folder}/step3 completed in {format_time(duration)}.\n")
                   
            end_cwd = os.getcwd()
            print(f"Starting '{output_folder}/step4' simulation in {end_cwd}.")
            
            command = f"maos -o {output_folder}/step4 -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}' 'powfs.noisy=0'"

            print(command + "\n")
                                                                         
            # Record the duration of the simulation
            start_time = time.time()
            if run_simulation:
                subprocess.run(command, shell=True, text=True)
                                                    
                end_cwd = os.getcwd()
                print(f"Finished '{output_folder}/step4' simulation in {end_cwd}.")

            end_time = time.time()
            duration = end_time - start_time

            # Get the current date and time when the simulation ends
            end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Print duration of simulation to a text file
            time_file.write(f"As of {end_datetime}: Simulation for {output_folder}/step4 completed in {format_time(duration)}.\n")
            print(f"As of {end_datetime}: Simulation for {output_folder}/step4 completed in {format_time(duration)}.\n")
            
            end_cwd = os.getcwd()
            print(f"Starting '{output_folder}/step5' simulation in {end_cwd}.")
                                                
            command = f"maos -o {output_folder}/step5 -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}' 'powfs.noisy=1'"

            print(command + "\n")
                                                                         
            # Record the duration of the simulation
            start_time = time.time()
            if run_simulation:
                subprocess.run(command, shell=True, text=True)
                                                    
                end_cwd = os.getcwd()
                print(f"Finished '{output_folder}/step5' simulation in {end_cwd}.")
                                                            
            end_time = time.time()
            duration = end_time - start_time

            # Get the current date and time when the simulation ends
            end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Print duration of simulation to a text file
            time_file.write(f"As of {end_datetime}: Simulation for {output_folder}/step5 completed in {format_time(duration)}.\n")
            print(f"As of {end_datetime}: Simulation for {output_folder}/step5 completed in {format_time(duration)}.\n")
            
        total_end_time = time.time()
        total_duration = total_end_time - total_start_time
        total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")
    
    elif output_folder == 'vismcao_single_study': 
        total_start_time = time.time()
        os.makedirs(output_folder, exist_ok=True)
        time_log_file = os.path.join(output_folder, 'maos_sim_time.txt')
                
        with open(time_log_file, 'a') as time_file:
            # Actuator 
            side_actuator_asm = calc_side(actuator_count=n_actuator_asm, sigfigs=3)
            dm_dx = create_dm_dx_array(n_dms=n_downstream_dms, side_dms=side_downstream_dms, side_asm=side_actuator_asm)

            # Wattage
            sodium_siglev, sodium_nearecon, sodium_bkgrnd, tt_siglev, tt_bkgrnd, tt_nearecon, truth_siglev, truth_bkgrnd, truth_nearecon = photometry(n_sodium=n_sodium, n_tiptilt=n_tiptilt, b_sodium=b_sodium, w_sodium=w_sodium, m_tiptilt=m_tiptilt, asm_side=side_actuator_asm, integration_time=integration_time)

            # Radius
            wfs_sodium_x, wfs_sodium_y, llt_sodium_x, llt_sodium_y, wfs_tiptilt_x, wfs_tiptilt_y = asterisms(n_sodium=n_sodium, r_sodium=r_sodium, n_tiptilt=n_tiptilt, r_tiptilt=r_tiptilt, on_axis_1_sodium_beacon=on_axis_1_sodium_beacon)

            command = f"maos -o {output_folder} -c {master_vismcao_file} plot.all=1 plot.setup=1 -O 'dm.dx = [{dm_dx}]' 'powfs.nwfs=[{n_sodium} {n_tiptilt} 1]' 'powfs.siglev = [{sodium_siglev} {tt_siglev} {truth_siglev}]' 'powfs.bkgrnd=[{sodium_bkgrnd} {tt_bkgrnd} {truth_bkgrnd}]' 'powfs.nearecon=[{sodium_nearecon} {tt_nearecon} {truth_nearecon}]' 'wfs.thetax=[{wfs_sodium_x} {wfs_tiptilt_x} 0]' 'wfs.thetay=[{wfs_sodium_y} {wfs_tiptilt_y} 0]' 'powfs0_llt.ox = [{llt_sodium_x}]*6.5' 'powfs0_llt.oy=[{llt_sodium_y}]*6.5' 'fit.fov={field_of_view}' 'tomo.maxit={tomo_iterations}' 'sim.end={simulation_steps}'"

            print(command + "\n")
                                                                         
            # Record the duration of the simulation
            start_time = time.time()
            if run_simulation:
                subprocess.run(command, shell=True, text=True)
                
                end_cwd = os.getcwd()
                print(f"Finished {output_folder} simulation in {end_cwd}.")

            end_time = time.time()
            duration = end_time - start_time

            # Get the current date and time when the simulation ends
            end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # Print duration of simulation to a text file
            time_file.write(f"As of {end_datetime}: Simulation for step5 completed in {format_time(duration)}.\n")
            print(f"As of {end_datetime}: Simulation for step5 completed in {format_time(duration)}.\n")
            
        total_end_time = time.time()
        total_duration = total_end_time - total_start_time
        total_end_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"As of {total_end_datetime}: Total simulation time: {format_time(total_duration)}.\n")
                           
        if run_simulation:
            os.chdir("..")
            os.chdir("..")
                                                
    else: 
        raise ValueError("Not one of the options: 'asm_actuator_study', 'sodium_wattage_study', 'r_sodium', 'sodium_count_study', 'dm_two_height_study', 'dm_three_height_study', 'error_budget', 'vismcao_single_study'. ")
                                            
    final_cwd = os.getcwd()
    print(f"Ending in {final_cwd}.")

if __name__ == "__main__":
    main()