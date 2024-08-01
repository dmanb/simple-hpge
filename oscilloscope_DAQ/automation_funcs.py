#### FUNCTIONS FOR DATA COLLECTION ####

import sys, os
#exec(open("./simple-hpge/oscilloscope_DAQ/automation_funcs.py").read()) ### annoying workaround to include autmation_funcs.py
#from automation_funcs import *
import pyvisa
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import struct 
import os
import datetime

# def addition(a, b):
#     return a + b

## configures power supply to desired format
def configure_power_supply():
    channel = 3
    voltage = 2.8
    current_limit = 50
    
    # Select the channel
    supply.write(f'INST:NSEL {channel}')
    
    # Set voltage
    supply.write(f'VOLT {voltage}')
    
    # Set current limit
    supply.write(f'CURR {current_limit}')
    
    # Turn on the output
    supply.write('OUTP ON')
    
    # Read back current
    current = float(supply.query('MEAS:CURR?'))
    
    # Check if the current is within the desired range
    if .005 <= current <= .010:
        print(f'Current is within range: {current} A')
    else:
        print(f'Current is out of range: {current} A')

#############################################################################

def pulser_v(voltage): ## parameter in units of volts
    return pulser.write(f'C1:BSWV AMP, {voltage}') ## set voltage of pulser


#############################################################################

def set_vert_scale(channel, scale):
    scope.write(f':CH{channel}:SCALE {scale}')

#############################################################################

def set_horizontal_scale(scale):
    scope.write(f':HOR:MAIN:SCALE {scale}')

#############################################################################

## gets measurement from any of the measurement channels on oscilloscope
def get_measurement(measurement_channel):
    # Query the value of the specified measurement channel
    measurement = scope.query(f':MEASU:MEAS{measurement_channel}:VALue?')
    return float(measurement)


#############################################################################

## make sure to adjust which channel you want to use for the trigger
def set_trigger(trigger_level):
    # Construct the SCPI command to set the trigger source
    scope.write('TRIGger:A:EDGE:SOURce CH4') ## makes sure the trigger is set to channel 4
        
    # Construct the SCPI command to set the trigger level
    scope.write(f'TRIGger:A:LEVel:CH4 {trigger_level}') 

#############################################################################


# getting data from waveform
def get_waveform_data(channel):
    # Select the channel
    scope.write(f'DATA:SOURCE CH{channel}')
    scope.write(f'DATA:WIDTH 2')
    
    raw_data = scope.query_binary_values('CURVE?', datatype = 'b', container=np.array)

# Horizontal scale settings to reconstruct time data
    x_increment = float(scope.query('WFMPRE:XINCR?'))
    x_zero = float(scope.query('WFMPRE:XZERO?'))
    x_offset = float(scope.query('WFMPRE:PT_OFF?'))
    # print(x_offset, x_zero)

# Get vertical scale settings to reconstruct voltage data
    y_multiplier = float(scope.query('WFMPRE:YMULT?'))
    y_offset = float(scope.query('WFMPRE:YOFF?'))
    y_zero = float(scope.query('WFMPRE:YZERO?'))

# Convert raw data to numpy array
    data = np.frombuffer(raw_data, dtype=np.int16)  # Adjust dtype if needed

# Reconstruct time and voltage data
    time_data = x_zero + (np.arange(len(data)) - x_offset) * x_increment
    voltage_data = y_zero + (data - y_offset)*y_multiplier

    return time_data, voltage_data


#############################################################################


def adjust_y_divisions(asic_amplitude):
    ## takes in amplitude of the ASIC (a float) and adjusts y divisions of oscilloscope as required
    if 1.6 < asic_amplitude <= 4:
        set_trigger(0.1)
        set_vert_scale(1, 0.5)
        set_vert_scale(4, 0.5)
    elif 0.8 <= asic_amplitude <= 1.6:
        set_vert_scale(1, 0.2)  # 200 mV/div for both ASIC readout, pulser
        set_vert_scale(4, 0.2)
    elif 0.4 <= asic_amplitude < 0.8:
        set_vert_scale(1, 0.1)  # 100 mV/div for ASIC, still 200 for pulser
    elif 0.18 <= asic_amplitude < 0.4:
        set_trigger(0.05)
        set_vert_scale(1, 0.05)  # 50 mV/div on ASIC readout
        set_vert_scale(4, 0.1)   # 100 mV/div for pulser
    elif 0.08 <= asic_amplitude < 0.18:
        set_trigger(0.05)
        set_vert_scale(1, 0.02)  # 20 mV/div
        set_vert_scale(4, 0.05)  # 50 mV/div for pulser
    elif 0.04 <= asic_amplitude < 0.08:
        set_trigger(0.02)
        set_vert_scale(1, 0.01)  # 10 mV/div
        set_vert_scale(4, 0.02)  # 20 mV/div for pulser
    elif 0.015 <= asic_amplitude < 0.04:
        set_trigger(0.01)
        set_vert_scale(1, 0.005)  # 5 mV / div
        set_vert_scale(4, 0.01)   # 10 mV/div for pulser
    elif .008 <= asic_amplitude < 0.015:
        set_trigger(0.004)
        set_vert_scale(1, 0.002)
    elif 0.0 <= asic_amplitude < .008:
        set_vert_scale(1, 0.001)


#############################################################################


## saves ASIC time data, voltage data, pulser voltage data to .csv file as well as other parameters
def save_to_csv(path, time_data, device_voltage_data, pulser_voltage_data, pulser_voltage, asic_voltage, counter = None):
    filename = f'ASIC_A{asic_voltage:.3f}V_P{pulser_voltage:.3f}V_{counter:04d}.csv'.replace('.','p',2)  
    fullpath = os.path.join(path, filename)
    
    with open(fullpath, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        csvwriter.writerow([f'Pulser Voltage: {pulser_voltage:.3f}V', f'ASIC Peak Voltage: {asic_voltage:.3f}V' ])
        csvwriter.writerow(['Time (s)', 'ASIC Voltage (V)','Pulser Voltage (V)'])  # Header row
        for t, v, q in zip(time_data, device_voltage_data, pulser_voltage_data):
            csvwriter.writerow([t, v, q])
   

#############################################################################


## plots ASIC waveform from a .csv
def plot_ASIC_waveform(filename):
    df = pd.read_csv(filename, skiprows = 1)
    
    # Print the column names to debug
    print("Column names:", df.columns)
    
    # Extract the columns for time and voltage
    time_column = df.columns[0]  # Assuming the first column is time
    voltage_column = df.columns[1]  # Assuming the second column is voltage of ASIC
    
    time_data = df[time_column]
    voltage_data = df[voltage_column]
    
    # Plot the waveform
    plt.figure(figsize=(10, 6))
    plt.plot(time_data, voltage_data, label='Waveform')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.title(f'Waveform Data from {filename}')
    plt.legend()
    plt.grid(True)
    plt.show()


#############################################################################


## plots pulser waveform from a .csv
def plot_pulser_waveform(filename):
    # Skip the first 3 rows (which contain metadata)
    df = pd.read_csv(filename, skiprows = 1)
    
    # Print the column names to debug
    print("Column names:", df.columns)
    
    # Extract the columns for time and pulser voltage
    time_column = 'Time (s)'  # Adjust if necessary after inspecting column names
    pulser_voltage_column = 'Pulser Voltage (V)'  # Adjust if necessary after inspecting column names
    
    time_data = df[time_column]
    pulser_voltage_data = df[pulser_voltage_column]
    
    # Plot the waveform
    plt.figure(figsize=(10, 6))
    plt.plot(time_data, pulser_voltage_data, label='Pulser Voltage')
    plt.xlabel('Time (s)')
    plt.ylabel('Pulser Voltage (V)')
    plt.title(f'Pulser Voltage Waveform from {filename}')
    plt.legend()
    plt.grid(True)
    plt.show()


#############################################################################


## creates folder for .csv files to be inserted into. The folder is automatically created in this directory
## in the folder I have called ./Data.
## GIVE THE FOLDER A NAME IN THE get_data FUNCTIONS BELOW
def create_folder(base_dir='../ASIC_data', name = None, run_type = None, n_runs = None, subfolder = None):
    # If a subfolder is specified, use it directly without generating a new folder name
    if subfolder:
        folder_path = os.path.join(base_dir, subfolder)
    else:
        # Generate a folder name if 'name' is not provided and 'subfolder' is not specified
        if name is None:
            name = f'ASIC_data_{int(time.time())}'
        folder_path = os.path.join(base_dir, name)
    
    # Create the folder (or subfolder) path
    os.makedirs(folder_path, exist_ok=True)
    # print(f'Folder created: {folder_path}')
    
    ## create a .csv file with important metadata for main data collection folders
    if not subfolder:
        info_filename = 'run_info.csv'
        info_filepath = os.path.join(folder_path, info_filename)
        current_date = datetime.datetime.now().strftime('%Y-%m-%d')
        current_time = datetime.datetime.now().strftime('%H:%M:%S')

        with open(info_filepath, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(['Oscilloscope: TEKTRONIX MSO44B', 'Pulser: TELEDYNE T3AFG500'])
            csvwriter.writerow([f'Date: {current_date}', f'Time: {current_time}'])
            csvwriter.writerow(['Time units: microseconds', 'Voltage units: Volts'])


            if run_type:
                csvwriter.writerow([f'Run type: {run_type}'])

            if n_runs:
                csvwriter.writerow([f'# of iterations at each voltage: {n_runs}'])
    
    return folder_path


#############################################################################

### MAIN FUNCTION FOR DATA COLLECTION, FUNCTIONS BELOW FOR MORE SPECIFIC / SPECIALIZED DATA COLLECTION

def collect_data(voltage_array, num_runs, folder_title = None, base_directory = '../ASIC_data'):
    ## create folder title for this run with given title, type of run in this case iterates over an array 
    ## at each voltage multiple times. 
    folder_path = create_folder(base_dir = base_directory,
                                name = folder_title, 
                                run_type = f'Array of voltages, each iterated over', 
                                n_runs = num_runs)
    
    ## iterating over each voltage in the voltage array
    for volt in voltage_array:
        
        subfolder_name = f'ASIC_PulserVoltage_{volt:.3f}V'.replace('.','p', 1)
        subfolder_path = create_folder(base_dir = folder_path, subfolder = subfolder_name)
        
        counter = 0      ## to save multiple files at same voltage
        pulser_v(volt)   ## sets the voltage on the pulser
        time.sleep(.25)   ## adds delay for waveform to stabilize before data collection
        
        
        asic_amplitude = get_measurement(3) ## gets ASIC max amplitude so we can adjust scaling
        
        adjust_y_divisions(asic_amplitude)  ## adjusts y scaling for max resolution using asic amplitude
        
        ## collects n waveforms at each voltage
        for n in range(num_runs):
            time.sleep(.001)
            time_wave, voltage_wave               = get_waveform_data(1) ## gets ASIC waveform data out of channel 1
            time_wave_pulser, voltage_wave_pulser = get_waveform_data(4) ## gets pulser voltages 
            save_to_csv(subfolder_path, time_wave, voltage_wave, voltage_wave_pulser, volt, asic_amplitude, counter)
            counter += 1
    
    

#############################################################################

def initialize_scope_and_pulser():
    set_trigger(.015)

    ## pulser parameters
    rise_time = 10e-9 
    fall_time = 10e-9  ## 10 ns rise and fall time

    pulse_width = .1   ## 100 ms pulse width
    pulser_delay = 0.0 ## no delay

    set_horizontal_scale(100e-6) ## sets time axis of oscilloscope to 40 microseconds


    pulser.write(f':PULS:TRAN:LEAD {10e-9}') 
    pulser.write(f':PULS:TRAN:TRA {10e-9}')  ## rise and fall time to 10 nanos
    pulser.write(f':PULS:WIDth {.1}')

    pulser.write(f':PULS:DEL {0.0}') ## 0 delay 



    ## INITIALIZING SCOPE
    set_trigger(.1)
    pulser_v(4.5)
    set_vert_scale(4, .5)
    set_vert_scale(1, .5)


    ## STANDARDIZE SCOPE OFFSET SO ITS IDENTICAL ACROSS RUNS
    scope.write('CH1:POS -4.00')
    scope.write('CH4:POS -4.00')
    scope.write('HORIZONTAL:POS 20')


    ## SET ACQUISITION MODE TO HIGH RESOLUTION
    scope.write(':ACQ:MODE HIRES')

#############################################################################


def close_device_connections():
    scope.close()
    pulser.close()
    supply.close()
    # print('Connections closed')
