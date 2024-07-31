# %%
#### INCLUDING AUTOMATION FUNCTIONS AND PACKAGES NECESSARY FOR DATA COLLECTION ####
from automation_funcs import *
import pyvisa
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import struct 
import os
import datetime


#####################   CONNECT TO INSTRUMENTS   #############################

### resource manager for pyvisa, which talks to instruments

#%% 
rm = pyvisa.ResourceManager('@py')


## assign oscilloscope
scope = rm.open_resource('TCPIP0::169.254.7.109')


## assign pulse generator 
pulser = rm.open_resource('TCPIP0::169.254.7.110')

#%%

# %% 
## power supply is finnicky sometimes
ip_address = '169.254.7.111'
resource_string = f'TCPIP0::{ip_address}::5025::SOCKET'

try:
    # Open the resource with the specified IP address and port
    supply = rm.open_resource(resource_string)

    # Set the read and write termination characters
    supply.write_termination = '\n'
    supply.read_termination = '\n'
    
    # Set the timeout to 5000 ms (5 seconds)
    supply.timeout = 5000

    # Query the device for its ID
    idn = supply.query("*IDN?")
    print("Connected to power supply:", idn)
except pyvisa.VisaIOError as e:
    print("Error connecting to the power supply:", e)
# %%

##############################################################

###########   INITIALIZE OSCILLOSCOPE, PULSE GENERATOR, POWER SUPPLY  ##############

initialize_scope_and_pulser()
configure_power_supply()



