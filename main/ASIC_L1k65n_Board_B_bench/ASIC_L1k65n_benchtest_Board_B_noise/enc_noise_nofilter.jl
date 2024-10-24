# rough estimation of enc noise without filtering
# serve as sanity check 
relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
using Unitful
using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using ArraysOfArrays
using Measures, Measurements
using Measurements: value as mvalue
using Statistics, Distributions
using Plots
using Printf, LaTeXStrings
using JLD2
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

# settings 
recompute = true
set = (Cinj_fF = 3000, Temp_K = 300, noiseAmpFac = 100, Cf_fF = 500)

# load or compute ENC 
pd_data = readlprops("$(@__DIR__)/data_config.json")
fpath = "$(@__DIR__)$(pd_data.noisefolder)"
fname = fpath * "Baseline_noise_simple_Cinj$(set.Cinj_fF)fF_T$(set.Temp_K)K_AmpFac$(set.noiseAmpFac).jld2"

if isfile(fname) & !recompute 
    f = jldopen(fname)
    bl_rms_e = f["bl_rms_e"]
    bl_rms_V = f["bl_rms_V"]
    bl_std_e = f["bl_std_e"]
    bl_std_V = f["bl_std_V"]
    @info "Reading baseline flucutations rom file: \n$fname"
    close(f)
else
    # read data
    # dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)
    dataset = noise_dataset_str(set.Cinj_fF, set.Temp_K)
    folder = pd_data.datafolder * dataset * "/"
    wvfs_raw, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1);
    @sprintf("Reading %d waveforms ", length(wvfs_raw))

    # shift waveforms
    bl_mean = mean.(wvfs_raw.signal)
    wvfs = shift_waveform.(wvfs_raw, -bl_mean)
    
    # calculate baseline rms/std --> simple proxy for ENC noise without any filtering 
    csa_gain = set.Cinj_fF / set.Cf_fF
    wvfs_bl_all = vcat(wvfs.signal...) / set.noiseAmpFac / csa_gain
    bl_rms_V = sqrt(sum(wvfs_bl_all.^2)/length(wvfs_bl_all))
    bl_std_V = std(wvfs_bl_all)
    bl_rms_e = csa_voltage2charge(set.Cinj_fF, bl_rms_V) # convert to electron equivalent
    bl_std_e = csa_voltage2charge(set.Cinj_fF, bl_std_V)# 
    

    if !ispath(fpath)
        mkpath(fpath)
    end

    jldsave(fname; bl_rms_V, bl_rms_e, bl_std_V, bl_std_e, set)
    println("Writing simple ENC results to file: \n$fname") 

    # plot
    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 14,
    legendtitlefontsize = 16,
    colorbar_titlefontsize = 16,
    legendforegroundcolor = :transparent)

    hist = fit(Histogram, 1e3.*wvfs_bl_all, nbins = 1000)
    gaussian = x -> pdf.(Normal(mean(wvfs_bl_all)*1e3, std(wvfs_bl_all))*1e3, x)
    xs = 1e3.*range(minimum(wvfs_bl_all), stop = maximum(wvfs_bl_all), length = 100)

    plot(hist, linewidth = 0, color = :silver, 
    label = @sprintf("%s shifted baselines \nrms = %.3g mV (%.0f e-)", length(wvfs), bl_rms_V * 1e3, bl_rms_e) )
    plot!(xs, gaussian.(xs)*length(wvfs_bl_all)*step(hist.edges[1]), 
    xlabel = "CSA-in baseline noise (mV)",
    linewidth = 2, color = :dodgerblue, 
    label = false,
    legend = :topleft,
    ylims = (0, maximum(hist.weights)*1.45), #xlims = (-50, 50),
    yformatter = :plain,
    title =  "Baseline noise L1k65n benchtest \n CSA Board-B: Cinj = $(set.Cinj_fF)fF, Cf = $(set.Cf_fF)fF,  T = $(set.Temp_K)K \n",
    titlefontsize = 12,
    size = (600, 450),
    dpi = 300)
    ppath = "$(@__DIR__)$(pd_data.figurefolder)Noise/"
    if !ispath(ppath)
    mkpath(ppath)
    end
    pname = replace(replace(fname, fpath => ""), ".jld2" => ".png")
    savefig(ppath * pname)
    @info "Saving plot to file:  figures/$pname"
    display(plot!())
end 
