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
using DataFrames
using Measurements: value as mvalue, uncertainty as muncertainty
using LegendDataManagement: readlprops, writelprops
using CSV


#### including the util_funcs.jl file 
include("$(@__DIR__)/../bench_test/util_funcs.jl")


### give it your path from where you are now
rel_path = "../../ASIC_data_folder/ASIC_data_test_Jul24"


### assign the return values of the read_csv function to the variables wvfs and decays
wvfs = read_folder(rel_path)
wvfs
print(typeof(wvfs))

wvf1 = wvfs[5]
plot(wvf1)

plotlyjs()
plot(wvf1)

### get decay times from wvfs from csv files 
decay_times = get_decay_times(wvfs)

### getting information on waveforms from simple_dsp function 
table = simple_dsp(wvfs, decay_times)



# columnnames(table)

## getting median decay time from data set for usage in the dsp_trap_rt_optimization function
taus = median(table.tail_τ)


## loading dsp_config locally since it is needed in dsp_trap_rt_optimization function
path_config = "$(@__DIR__)/../config/dsp_config.json"
# get DSP configuration data --> Can be modified in .json filepath = " $(@__DIR__)/../config/dsp_config.json"
dsp_meta = readlprops(path_config)
dsp_config = DSPConfig(dsp_meta.default)
dsp_config.bl_window
dsp_config.tail_window  


## running the dsp_trap_rt_optimization function to get the ENC vs. shaping time plot
enc_trap_grid = dsp_trap_rt_optimization(wvfs, dsp_config, taus,; ft=1.0u"µs")


ENC = Vector{Float64}()
for i in 1:length(dsp_config.e_grid_rt_trap)
    enc_i = flatview(enc_trap_grid)[i,:]
    push!(ENC, std(enc_i))
end 




plot(ustrip.(collect(dsp_config.e_grid_rt_trap)), ENC, xlabel="Rise Time [µs]", ylabel="ENC [e-]", title="ENC vs. RT", legend=:outertopright)


i = 11
enc_i = flatview(enc_trap_grid)[i,:]