
relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
using Plots
using Printf, LaTeXStrings
using PropDicts
using Unitful, Measurements, Measures
using LegendDataManagement: readlprops
using Statistics
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/plot_waveforms.jl")


# get configs (data, dsp )
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
dsp_file = dsp_folder * "DSPpars.json"

plt_path = "$(@__DIR__)" *pd_data.figurefolder * "DSP/"
ispath(plt_path) || mkpath(plt_path)
(csa, board) =  CSAname(pd_data)

# load dsp results 
if isfile(dsp_file)
    dsp_par = readlprops(dsp_file)
    println("Reading DSP results from file: $dsp_file")  
else
    println("DSP results not found") 
end

default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
xguidefontsize = 16, xtickfontsize = 12,
yguidefontsize = 16, ytickfontsize = 12,
legendfontsize = 12, titlefontsize = 10,
legendforegroundcolor = :silver)

#  load quality cuts (qc)
qc_folder = "$(@__DIR__)$(pd_data.dspfolder)"
qc_file = qc_folder * "QualityCuts.json"
cuts_qc = readlprops(qc_file).wvf_keep.all
wvf_qc = collect(1:length(cuts_qc))[cuts_qc] # indices of waveforms that surved qc

# decay time 
p2 = plot_dsp_par(ustrip.(dsp_par.tail_τ), nbins = 200, 
        xlbl = "Decay time", xunit = "$(unit(dsp_par.tail_τ[1]))",
         txt = true)
plot!(title = "CSA $csa Board $board with Ge detector", titlefontsize = 12)
savefig(p2, plt_path * "decaytime.png")

# tail slope 
p3 = plot_dsp_par(1e7 .* ustrip.(dsp_par.tailslope), nbins = 200, 
        xlbl = "Tail slope", xunit = L"$\,\times 10^{-7} \,\textrm{V}^{-1}$",
        lbl = false, txt = true)
plot!(title = "CSA $csa Board $board with Ge detector", titlefontsize = 12)
savefig(p3, plt_path * "tailslope.png")
         

# rise-time
risetime_ns = ustrip.(uconvert.(u"ns", dsp_par.t90 .- dsp_par.t10))
risetime_max = 500
prt = plot_dsp_par(filter(x-> x<=risetime_max, risetime_ns), nbins = 100, 
        xlbl = "Rise time " * L"$t_{90} - t_{10}$", xunit = "ns",
        lbl = false, 
        txt = true)
plot!(title = "CSA $csa Board $board with Ge detector \n" * @sprintf("Outlier removed: %.0f%% have risetime > %.0f ns", 
                100*length(filter(x-> x>risetime_max, risetime_ns))/length(risetime_ns),risetime_max), titlefontsize = 12)
savefig(prt, plt_path * "risetime.png")
    
