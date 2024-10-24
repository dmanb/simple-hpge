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
include("$(@__DIR__)/../../src/simple_dsp.jl")
recompute = false 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

dspfolder = "$(@__DIR__)$(pd_data.dspfolder)"
dsp_file = dspfolder * "dsp_pars_" * replace(pd_data.datafolder , "/" => "")  * ".json"

if isfile(dsp_file) && !recompute
    dsp_par = readlprops(dsp_file)
    println("Reading DSP results from file: $dsp_file")
else
    # read data 
    folder = pd_data.datafolder
    wvfs = read_folder_csv_oscilloscope(folder; heading = 17, fixtime = true)

    # dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
    dsp_par = simple_dsp(wvfs, dsp_config)

    # create a PropDict save results as json file
    dsp_pd = PropDict(Dict(varname => getproperty(dsp_par, varname) for varname in columnnames(dsp_par)))
    if !ispath(dspfolder)
        mkpath(dspfolder)
    end
    writelprops(dsp_file, dsp_pd)
    println("Writing DSP results to file: $dsp_file")
end 

# #########  look at results #########
 # some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
 xguidefontsize = 16, xtickfontsize = 12,
 yguidefontsize = 16, ytickfontsize = 12,
 legendfontsize = 14, titlefontsize = 10,
 legendforegroundcolor = :silver)

vars = keys(dsp_par) # overview of variables in dsp_par
nbins = round(Int, 0.25*length(dsp_par.t90))

# decay time
p1 = stephist(dsp_par.tail_τ, bins= nbins, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))

# rise time 
p2 = stephist(1e3 .* ustrip.(dsp_par.t90 .- dsp_par.t10), 
    bins= nbins, xlabel= "Rise time (ns)", ylabel="Occurrence", 
    xformatter = x -> @sprintf("%.0f", x),
    label = L"$t_{90} - t_{10}$", 
    fill = true, color = :red2, ylims = (0, :auto))

# baseline std
p3 = stephist(1e3.*filter(x-> x<= quantile(dsp_par.blsigma, 0.95), dsp_par.blsigma), bins=20, xlabel= "Baseline σ (mV)", ylabel="Occurrence", 
label = false, fill = true, color = :orange, ylims = (0, :auto))

# energy with trapezoidal filter 
p4 = stephist(1e3.*dsp_par.e_trap, bins = nbins, 
    xlabel= "Amplitude (mV)", ylabel="Occurrence", 
    label = @sprintf("Trap filter \nσ = %.1f mV",std(1e3.*dsp_par.e_trap)), 
    fill = true, 
    color = :dodgerblue, 
    ylims = (0, :auto))

# energy with cusp filter
p5 = stephist(1e3.*dsp_par.e_cusp, bins = nbins, 
    xlabel= "Amplitude (mV)", ylabel="Occurrence", 
    label = @sprintf("Cusp filter \nσ = %.1f mV",std(1e3.*dsp_par.e_cusp)), fill = true, 
    color = :red2, ylims = (0, :auto))

# overview 
plot(p1, p2, p3, p4, layout = (2,2), size = (1200, 800), 
    left_margins = 2mm , thickness_scaling = 1.4,
    dpi = 300)

# # save figure 
fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
mkpath(fpath)
fname = fpath * "DSPpars.png"# * replace(pd_data.datafolder, "/" => "") 
savefig(fname)

