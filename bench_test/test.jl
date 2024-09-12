Pkg.activate("$(@__DIR__)/../")
import Pkg
using TypedTables, Dates
using RadiationDetectorSignals
using Unitful, Formatting, LaTeXStrings, Measures, Measurements
using Measurements: value as mvalue, uncertainty as muncertainty
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP
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
using TypedTables, Dates
using DataFrames
using Measurements: value as mvalue, uncertainty as muncertainty
using LegendDataManagement: readlprops, writelprops
using CSV


#### including the util_funcs.jl file 
include("$(@__DIR__)/../bench_test/util_funcs.jl")

## print your current directory to help you understand what the path is to the data
print(@__DIR__)

## use read_data_folder function to get all waveforms from data folder
## read_data_folder function is defined in util_funcs.jl
wvfs_dict = read_data_folder("../ASIC_data/ASIC_BenchTest_08082024/")

## print the keys in the dictionary
println("Available keys in wvfs_dict: ", keys(wvfs_dict))

## get the waveforms from the dictionary
wvfs = wvfs_dict["ASIC_PulserVoltage_4p500V"]


wvf = wvfs[1]

t0 = get_t0(wvf, t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)




typeof(wvfs)

test = ArrayOfRDWaveforms
typeof(test)




## pick the first waveform and plot it (might take a second)
wvf1 = wvfs[1:10]
plotlyjs()
plot(wvf1, ylabel = "Volts (V)")

plot(wvfs[432])
plot!(wvfs[431])


### get decay times from wvfs from csv files 
decay_times = get_decay_times(wvfs)



## code to get rid of last waveform in wvfs, since it was returning an error
pop!(wvfs)
pop!(decay_times)

## test pole zero correction (which occurs already in the simple_dsp function, but can be done separately)

deconv_flt = InvCRFilter.(decay_times)

## pz_wavs = pole zero corrected waveforms


pz_wavs = [deconv_flt[i](wvfs[i]) for i in eachindex(wvfs)]


## plotting the fifth waveform before and after pole zero correction
plot(pz_wavs[5], label="After Pole Zero Correction",ylabel = "Volts (V)")
plot!(wvfs[5], label="Before Pole Zero Correction", ylabel = "Volts (V)")
plot!(legend =:right)
display(plot)



### getting information on waveforms from simple_dsp function 
table = simple_dsp(wvfs, decay_times)



## getting column names of the table to see what information is available
# columnnames(table)

## parameter values for the trap filter gotten from simple_dsp table
trap = table.e_trap



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

## creates ENC noise vs. shaping time graph
plot(ustrip.(collect(dsp_config.e_grid_rt_trap)), ENC, xlabel="Rise Time [µs]", ylabel="ENC [e-]", title="ENC vs. RT", legend=:outertopright)



########### RUNNING DSP FUNCTION STEP BY STEP ####################






   ### config stuff
   path_config = "$(@__DIR__)/../config/dsp_config.json"
   # get DSP configuration data --> Can be modified in .json filepath = " $(@__DIR__)/../config/dsp_config.json"
   dsp_meta = readlprops(path_config)
   config = DSPConfig(dsp_meta.default)
   config.bl_window
   config.tail_window


   bl_window                = config.bl_window
   t0_threshold             = config.t0_threshold
   tail_window              = config.tail_window
   inTraceCut_std_threshold = config.inTraceCut_std_threshold
   sg_flt_degree            = config.sg_flt_degree
   current_window           = config.current_window
   qdrift_int_length        = config.qdrift_int_length
   lq_int_length            = config.lq_int_length

   #### all filter parameters set to default since we do not have access to full library 
   trap_rt, trap_ft = get_fltpars(PropDict(),:trap, config)
   cusp_rt, cusp_ft = get_fltpars(PropDict(), :cusp, config)
   zac_rt, zac_ft   = get_fltpars(PropDict(), :zac, config)
   sg_wl            = get_fltpars(PropDict(), :sg, config)

   # get CUSP and ZAC filter length and flt scale
   flt_length_zac              = config.flt_length_zac
   zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))
   flt_length_cusp             = config.flt_length_cusp
   cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

   # get number of samples the waveform is saturated at low and high of FADC range
   # bit_depth = config.kwargs_pars.fc_bit_depth # of FlashCam FADC
   # sat_low, sat_high = 0, 2^bit_depth - bit_depth
   # sat_stats = saturation.(wvfs, sat_low, sat_high)

   # set tau for CUSP filter to very high number to switch of CR filter
   τ_cusp = 10000000.0u"µs"
   τ_zac = 10000000.0u"µs"


   pop!(wvfs)
   pop!(decay_times)

   ################## ACTUAL WAVEFORM FILTERING AND RECONSTRUCTION, ANALYSIS BEGINS HERE ##################
   bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

   # substract baseline from waveforms
   wvfs = shift_waveform.(wvfs, -bl_stats.mean)
   
   # get wvf maximum
   wvf_max = maximum.(wvfs.signal)
   wvf_min = minimum.(wvfs.signal)
   
   # extract decay times
   tail_stats = tailstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))
   
   # deconvolute waveform 
   # --> wvfs = wvfs_pz
   deconv_flt = InvCRFilter.(decay_times)

   wvfs = [deconv_flt[i](wvfs[i]) for i in eachindex(wvfs)]
   
   # get tail mean, std and slope
   pz_stats = signalstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))
   

   t0 = get_t0(wvfs, t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)

   t0 = get_t0(wvfs[1], t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)