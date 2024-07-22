Pkg.activate("$(@__DIR__)/../")
import Pkg
using TypedTables, Dates
using RadiationDetectorSignals
using Unitful, Formatting, LaTeXStrings, Measures, Measurements
using Measurements: value as mvalue, uncertainty as muncertainty
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP: get_fltpars
using ArraysOfArrays
using Distributed
using RadiationDetectorDSP
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendHDF5IO, LegendDSP, LegendSpecFits, LegendDataTypes
using IntervalSets, TypedTables, StatsBase, PropDicts
using Unitful, Formatting, Printf, Measures, Dates
using Distributed, ProgressMeter, TimerOutputs
using RadiationDetectorDSP
using ArraysOfArrays
using HDF5
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP: get_fltpars
using TypedTables, Dates
using Measurements: value as mvalue, uncertainty as muncertainty
using LegendDataManagement: readlprops, writelprops
using CSV


#### including the util_funcs.jl file 
include("$(@__DIR__)/../bench_test/util_funcs.jl")


### give it your path from where you are now
rel_path = "../../ASIC_data"


### assign the return values of the read_csv function to the variables wvfs and decays
wvfs = read_csv(rel_path)
wvfs
print(typeof(wvfs))

### get decay times from wvfs from csv files 
decay_times = get_decay_times(wvfs)

### getting information on waveforms from simple_dsp function 
table = simple_dsp(wvfs, decay_times)

