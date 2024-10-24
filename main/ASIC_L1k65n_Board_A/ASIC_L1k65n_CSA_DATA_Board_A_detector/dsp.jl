import Pkg
Pkg.activate("$(@__DIR__)/../../../")
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
using JLD2
include("$(@__DIR__)/../../../utils/utils_IO.jl")
include("$(@__DIR__)/../../../src/simple_dsp.jl")
include("$(@__DIR__)/../../../utils/utils_naming.jl")
include("$(@__DIR__)/../../../utils/utils_aux.jl")
include("$(@__DIR__)/../../../src/data_quality.jl")
include("$(@__DIR__)/../../../utils/utils_plot.jl")

# settings 
plotFlag = true
recompute = true  

# get configs (data, dsp )
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# naming; save dsp file as json file to recover later
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
dsp_prelim_file = dsp_folder * "DSPpars_prelim.jl2"
dsp_file = dsp_folder * "DSPpars.json"
qc_folder = "$(@__DIR__)$(pd_data.dspfolder)"
qc_file = qc_folder * "QualityCuts.json"

# 1. Compute or load DSP parameter of all waveforms 
if isfile(dsp_prelim_file) 
    dsp_pd_prelim = jldopen(dsp_prelim_file)["dsp_pd"]
    @info "Reading prelim. DSP results from file: $dsp_file"
else
   @info "Preliminary DSP results file does not exist: $dsp_file \n Run quality cuts (qc first)"
end

# 2. Load quality cuts based on step 10000000.0u
if isfile(qc_file) 
    cuts = readlprops(qc_file)
    @info "Reading quality cuts to file: $qc_file"
else
     @info "Quality cuts file does not exist: $qc_file \n Run quality cuts (qc first)"
end

# apply quality cuts to prelim DSP ans save actual DSP json file 
if isfile(dsp_file) && !recompute
    dsp_qc = readlprops(dsp_file)
    println("Reading DSP results from file: $dsp_file")  
else
    dsp_qc = PropDict(Dict(varname => getproperty(copy(dsp_pd_prelim), varname)[cuts.wvf_keep.all] for varname in keys(copy(dsp_pd_prelim))))
    writelprops(dsp_file, dsp_qc)
    println("Writing DSP results to file: $dsp_file") 
end


(csa, board) =  CSAname(pd_data)
p1 = plot_dsp_par(dsp_qc, :e_trap, nbins = 1000, 
        xlbl = "Energy", xunit = "a.u.", 
        lbl = "Trap. filter", txt = false)
plot!(p1, title = "CSA $csa Board $board with Ge detector", titlefontsize = 12)

fpath = "$(@__DIR__)$(pd_data.figurefolder)Spectrum/"
if !ispath(fpath)
    mkpath(fpath)
end
fname = fpath * "Spectrum_QC_TrapRaw.png"
savefig(fname)
@info "Saving figure to $fname"

