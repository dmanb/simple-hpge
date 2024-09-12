import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using Measures 
using Statistics
using ArraysOfArrays
include("$(@__DIR__)/../utils_IO.jl")
include("$(@__DIR__)/../simple_dsp.jl")

# get dsp configuration
dsp_config_path = "$(@__DIR__)/../../config/dsp_config.json"
dsp_config = DSPConfig(readlprops(dsp_config_path).default)

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         legendforegroundcolor = :silver)

# read data 
filename = "GFET_detector_run_test-001.hdf5"

wvfs = read_wvfs_h5(filename; nwvfs = 10)
plot(wvfs.time, wvfs.signal, ylabel = "Volts (V)", label = "ASIC Waveform")

# dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
dsp_par = simple_dsp(wvfs, dsp_config)