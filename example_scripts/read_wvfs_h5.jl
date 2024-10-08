import Pkg
Pkg.activate("$(@__DIR__)/../")
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
using Plots
include("$(@__DIR__)/../utils/utils_IO.jl")
include("$(@__DIR__)/../src/simple_dsp.jl")

# get dsp configuration
config_path = "$(@__DIR__)/../config/"
dsp_config = DSPConfig(readlprops(config_path * "dsp_config.json" ).default)

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         legendforegroundcolor = :silver)

# read data 
folder = "GFET detector runs/"
filename = "GFET_detector_run_test.hdf5"
wvfs = read_wvfs_h5(filename; folder = folder, nwvfs = 5)
plot(wvfs.time, wvfs.signal, ylabel = "Volts (V)")

