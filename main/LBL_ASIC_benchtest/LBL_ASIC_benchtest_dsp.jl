import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using Measures 
using Statistics
using Plots
include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../src/simple_dsp.jl")

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         legendforegroundcolor = :silver)

# read data 
wvfs = read_folder_csv(pd_data.datafolder; heading = 2)

# dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
dsp_par = simple_dsp(wvfs, dsp_config)

# create a PropDict save results as json file
dsp_path = "$(@__DIR__)" * pd_data.dspfolder
if !ispath(dsp_path)
    mkpath(dsp_path)
end
dsp_pd = PropDict(Dict(varname => getproperty(dsp_par, varname) for varname in columnnames(dsp_par)))
writelprops(dsp_path * "dsp_pars.json", dsp_pd)
dsp = readlprops(dsp_path * "dsp_pars.json") # test if reading is working

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
