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
using StatsBase
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/../../utils/utils_IO.jl")

# load dsp parameter from file 
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dspfolder = "$(@__DIR__)$(pd_data.dspfolder)"
dsp_file = dspfolder * "dsp_pars_" * replace(pd_data.datafolder , "/" => "")  * ".json"
if !isfile(dsp_file)
    println("Run dsp.jl first.  DSP pars not found in file : $dsp_file")
else
    dsp_par = readlprops(dsp_file)
    println("Reading DSP results from file: $dsp_file")
end

default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
 xguidefontsize = 16, xtickfontsize = 12,
 yguidefontsize = 16, ytickfontsize = 12,
 legendfontsize = 14, titlefontsize = 10,
 legendforegroundcolor = :silver)
 
# get Standard deviation of amplitudes for different number of waveforms
# waveforms are randomly drawn (without replacement) from the total number of waveforms
nSamples = 5:1:length(dsp_par.e_trap)
stds = [ std(dsp_par.e_trap[sample(1:length(dsp_par.e_trap), i, replace=false)]) for i in nSamples]
p2 = plot(nSamples, 1e3*stds, xlabel = "Number of waveforms", ylabel = "Amplitude σ (mV)", linewidth = 2, color = :dodgerblue, label = "")
hline!([std(1e3.*dsp_par.e_trap)], label = @sprintf("All waveforms: σ = %.3f",std(1e3.*dsp_par.e_trap)), 
        linewidth = 3, linestyle = :solid, color = :orange, size = (500, 320), dpi = 300)

p3 = plot(nSamples, 100*stds./median(dsp_par.e_trap), xlabel = "Number of waveforms", ylabel = "Rel. amplitude σ (%)", linewidth = 2, color = :dodgerblue, label = "")
hline!([100*std(dsp_par.e_trap)./median(dsp_par.e_trap)], label = @sprintf("All waveforms: σ = %.3f%%",100*std(dsp_par.e_trap)./median(dsp_par.e_trap)), 
        linewidth = 3, linestyle = :solid, color = :orange, size = (500, 320), dpi = 300)


nbins = round(Int, 0.25*length(dsp_par.t90))
# energy with trapezoidal filter 
p1 = stephist(1e3.*dsp_par.e_trap, bins = nbins, 
    xlabel= "Amplitude (mV)", ylabel="Occurrence", 
    label = @sprintf("Trap filter \nσ = %.1f mV",std(1e3.*dsp_par.e_trap)), 
    fill = true, 
    color = :dodgerblue, 
    ylims = (0, :auto),)

plot(p1, p2, p3,  layout = (3,1), size = (650, 1000), legend = :topleft, dpi = 300)

# # save figure 
fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
mkpath(fpath)
fname = fpath * "Amplitude_vs_Samples.png"# * replace(pd_data.datafolder, "/" => "") 
savefig(p2, fname)


# # save figure 2 
fname3 = fpath * "RelAmplitude_vs_Samples.png"
savefig(p3, fname3)