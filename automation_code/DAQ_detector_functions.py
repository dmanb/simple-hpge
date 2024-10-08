# CODE FOR AUTOMATION OF DATA COLLECTION Germanium detector 
## LBNL lab 70-0141
## By Damien Bowen , updated: Lisa Schlueter September 2024
## packages
import shutil
import pyvisa
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import struct 
import os
from datetime import datetime
import h5py

## initializing resource manager for the electronics
rm = pyvisa.ResourceManager('@py')

## connecting to the oscilloscope (essential for both detector tests and bench tests)
scope = rm.open_resource('TCPIP0::169.254.7.109')
print(scope.query('*IDN?')) ## prints make and model of oscilloscope to verify connection

# getting waveform data from oscilloscope for a given channel
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

## data to write a waveform into an hdf5 file
def write_waveform_to_hdf5(hdf_file, group_name, time_data, voltage_data, run_num):
    """
    Writes waveform data to HDF5 file under the specified group.

    Parameters:
    hdf_file (h5py.File): The open HDF5 file object.
    group_name (str): The name of the group (e.g., voltage level).
    time_data (list or np.array): The time data points.
    voltage_data (list or np.array): The corresponding voltage data points.
    run_num (int): The current run number to name the dataset.
    """
    # Combine the time and voltage data into a 2D array
    waveform_data = np.vstack((time_data, voltage_data)).T
    
    # Create a dataset for the waveform in the specified group
    dataset_name = f"waveform_{run_num:05d}"
    
    # Write the dataset with compression
    hdf_file[group_name].create_dataset(dataset_name, data=waveform_data, compression="gzip")
    
    """Optional print statement for debugging and tracking."""
    # print(f"Waveform {run_num + 1} saved to group {group_name}.")
    
    
    