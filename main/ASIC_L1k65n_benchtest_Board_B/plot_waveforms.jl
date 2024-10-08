import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using ArraysOfArrays
using Measures 
using Statistics
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../utils/utils_naming.jl")
# settings 
recompute = false 
Cinj_fF = 500 # capacitance of the ASIC in femto Farrad 
Rf_MOhm = 500 # feedback resistor in Mega Ohm
PulserChargeInj_keV = 300 # pulser charge injected in keV
Temp_K = 77 # temperature in Kelvin 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# read data 
folder = benchtest_filename_raw(pd_data, Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
if !ispath(get(ENV, "ASIC_DATA",".") * folder)
    println("Data folder does not exist: $(get(ENV, "ASIC_DATA",".") * folder)")
end
wvfs_ch1, wvfs_ch2, MetaData = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 2)

# plot 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
 xguidefontsize = 16, xtickfontsize = 12,
 yguidefontsize = 16, ytickfontsize = 12,
 legendfontsize = 14, titlefontsize = 10,
 legendforegroundcolor = :silver)

IdxWvf = 1 
plot(wvfs_ch1[IdxWvf].time, wvfs_ch1[IdxWvf].signal, label = MetaData["Ch1"]["Channel"])
plot!(wvfs_ch2[IdxWvf].time, wvfs_ch2[IdxWvf].signal, label = MetaData["Ch2"]["Channel"])
plot!(xlabel = "Time ($(unit(collect(wvfs_ch1[1].time)[1])))", ylabel = "Voltage (V)", title = "Waveform $(IdxWvf)", legend = :topright)

