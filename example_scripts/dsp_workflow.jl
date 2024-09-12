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
using Plots
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

# dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
dsp_par = simple_dsp(wvfs, dsp_config)


#########  look at results #########
vars = columnnames(dsp_par) # overview of variables in dsp_par

# decay time
p1 = stephist(dsp_par.tail_τ, bins=20, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))

# energy with trapezoidal filter 
p2 = stephist(dsp_par.e_trap, bins=20, xlabel= "Energy (a.u.)", ylabel="Occurrence", label = "Trap filter", fill = true, color = :dodgerblue, ylims = (0, :auto))

# energy with cusp filter
p3 = stephist(dsp_par.e_cusp, bins=20, xlabel= "Energy (a.u.)", ylabel="Occurrence", label = "Cusp filter", fill = true, color = :red2, ylims = (0, :auto))

# overview 
plot(p1, p2, p3, layout = (3,1), size = (800, 1400), margins = 3mm, thickness_scaling = 1.3)
