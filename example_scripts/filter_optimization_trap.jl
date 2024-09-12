import Pkg
Pkg.activate("$(@__DIR__)/../")
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
include("$(@__DIR__)/../utils/utils_IO.jl")
include("$(@__DIR__)/../src/simple_dsp.jl")

# get dsp configuration
dsp_config_path = "$(@__DIR__)/../config/dsp_config.json"
dsp_config = DSPConfig(readlprops(dsp_config_path).default)

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         legendforegroundcolor = :silver)

# read data 
folder = "ASIC_data_07312024/ASIC_PulserVoltage_0p500V/"
wvfs = read_folder_csv(folder; heading = 2)

# do dsp with default settings
dsp_par = simple_dsp(wvfs, dsp_config)

## running the dsp_trap_rt_optimization function to get the ENC vs. shaping time plot
trap_rt, trap_ft = get_fltpars(PropDict(),:trap, dsp_config) # default rise-time and flattop-time from config
ft = 3.0*u"µs"#trap_ft #
enc_trap_grid = dsp_trap_rt_optimization(wvfs, dsp_config, median(dsp_par.tail_τ),; ft = ft)

ENC = Vector{Float64}()
for i in 1:length(dsp_config.e_grid_rt_trap)
    enc_i = flatview(enc_trap_grid)[i,:]
    push!(ENC, std(enc_i))
end 

## creates ENC noise vs. shaping time graph
plot(ustrip.(collect(dsp_config.e_grid_rt_trap)), ENC, xlabel="Shaping time (µs)", 
    linewidth = 2.5, color = :dodgerblue,
    ylabel="ENC noise (e-)", label = "Flattop time = $ft", legend=:bottomright, yformatter = :plain)
vline!([ustrip(trap_rt)], label = "Default shaping time = $trap_rt", color = :grey, linewidth = 2.5, linestyle = :dashdot)


