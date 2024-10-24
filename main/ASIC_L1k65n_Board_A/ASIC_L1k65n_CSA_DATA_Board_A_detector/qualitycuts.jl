# quality cuts
relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
using Unitful
using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDSP
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using ArraysOfArrays
using Measures 
using Statistics
using Plots
using Printf, LaTeXStrings
using JLD2
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/src/data_quality.jl")

# settings 
plotFlag = true
recomputeDSP_all = false 
recomputeQC = true 
chunkdata = true 

# get configs (data, dsp )
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# naming; save dsp file as json file to recover later
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
dsp_file = dsp_folder * "DSPpars_prelim.jl2"


# 1. Compute or load DSP parameter of all (non-corrupt) waveforms 
if isfile(dsp_file) && !recomputeDSP_all
    dsp_pd = jldopen(dsp_file)["dsp_pd"]
    println("Reading DSP results from file: $dsp_file")
else
    # read data 
    folder = pd_data.datafolder
    if chunkdata == true 
        filenames = filter!(x -> occursin(".csv", x), readdir(get(ENV, "ASIC_DATA",".") * folder))
        ntot = length(filenames)
        nwvfstep = 500
        niter = ceil(Int, ntot/nwvfstep)
    
        dsp_par = Vector{Table}(undef, niter)
        for i = 1:niter
            if i == niter
                (nwvfmin, nwvfmax) = (1 + nwvfstep * (i - 1), ntot)
            else
                (nwvfmin, nwvfmax) = (1 + nwvfstep * (i - 1), nwvfstep * i)
            end 

            local wvfs, _ = read_folder_csv_oscilloscope(folder; heading = 14, nChannels = 1, nwvfmax = collect(nwvfmin:1:nwvfmax)) # make sure you can read all waveforms 

            # dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
            dsp_par[i] = simple_dsp(wvfs, dsp_config)
        end

        dsp_par = vcat(dsp_par...)
    else
        wvfs, _ = read_folder_csv_oscilloscope(folder; heading = 14, nChannels = 1, nwvfmax = nwvfmax) # make sure you can read all waveforms 
        dsp_par = simple_dsp(wvfs, dsp_config)
    end
  
    # create a PropDict save results as json file
    dsp_pd = PropDict(Dict(varname => getproperty(dsp_par, varname) for varname in columnnames(dsp_par)))
    if !ispath(dsp_folder)
        mkpath(dsp_folder)
    end
    # writelprops(dsp_file, dsp_pd)
    jldsave(dsp_file; dsp_pd)
    println("Writing DSP results to file: $dsp_file")  
end

# 2. Apply and load quality cuts based
# save cuts to json file proper 
# save dsp results with quality cuts to json file proper
qc_folder = "$(@__DIR__)$(pd_data.dspfolder)"
qc_file = qc_folder * "QualityCuts.json"
if isfile(qc_file) && !recomputeQC
    cuts = readlprops(qc_file)
    println("Reading quality cuts to file: $qc_file")
else
    config_qc = readlprops("$(@__DIR__)/config/qc_config.json").default
    cuts = qc_dsp(dsp_pd, config_qc)
    writelprops(qc_file, cuts)
    println("Writing quality cuts to file: $qc_file")
end



