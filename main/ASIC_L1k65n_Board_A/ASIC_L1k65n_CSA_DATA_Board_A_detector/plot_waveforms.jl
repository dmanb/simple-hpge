relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
using Unitful
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
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")

# settings 
saveplot = true 

# get data configuration (where to and save stuff)
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

# plot function   
function plot_waveforms(pd_data; saveplot::Bool = false, wvf_sel::Vector{Int} = collect(1:6), qc::Bool = false, qc_par::Symbol = :all)
    # read data 
    folder = pd_data.datafolder
    if !ispath(get(ENV, "ASIC_DATA",".") * folder)
        println("Data folder does not exist: $(get(ENV, "ASIC_DATA",".") * folder)")
    end
    wvfs_ch1, _ = read_folder_csv_oscilloscope(folder; heading = 14, nChannels = 1, nwvfmax = wvf_sel) # make sure you can read all waveforms 

    if qc == true
        # apply quality cuts (qc)
        qc_folder = "$(@__DIR__)$(pd_data.dspfolder)"
        qc_file = qc_folder * "QualityCuts.json"
        if isfile(qc_file)
            cuts = readlprops(qc_file)
        else
            @info "Quality cuts file does not exist: $qc_file"
            return NaN
        end
        cuts_qc = cuts.wvf_keep[qc_par][wvf_sel]
        qc_pars = collect(keys(cuts.wvf_keep))[2:end]
        cuts_qc_pars = [cuts.wvf_keep[par][wvf_sel] for par in qc_pars]
        qc_str = "QC-$(qc_par)"
    else 
        qc_str = ""
    end 
    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 12, titlefontsize = 10,
    legendforegroundcolor = :silver)

    nwvf = length(wvfs_ch1)
    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, nwvf)
    for i = 1:nwvf
        if qc == true 
            wflbl = ""
            if cuts_qc[i] == true
                color = :dodgerblue
                wflbl = "passed"
            else
                color = :red
                cut_lbl = "cut by: "
                for p in eachindex(cuts_qc_pars)
                    if cuts_qc_pars[p][i] == true 
                        if p > 1
                            cut_lbl = cut_lbl * ", "
                        end
                        cut_lbl = cut_lbl * "$(qc_pars[p])"
                    end
                end 
                wflbl = wflbl * cut_lbl#" (cut, $qc_par)"
            end
        else
            wflbl = "Waveform"
            color = :dodgerblue
        end
        p[i] = plot(wvfs_ch1[i].time, wvfs_ch1[i].signal, label = wflbl, color = color)
        plot!(xlabel = "Time ($(unit(collect(wvfs_ch1[1].time)[1])))", ylabel = "Voltage (V)", 
        legend = :topright)
        if qc == true 
            # if !(i == findfirst(cuts_qc) || i == findfirst(.!cuts_qc))
                # only for first cut and uncut waveform 
                # plot!(legend = false)
            # end
        elseif i >1
            plot!(legend = false)
        end 
    end 

    (csaname, board) = CSAname(pd_data)
    ptot = plot(p..., layout = (ceil(Int,nwvf/2), 2), 
            size = (950, 900), 
            left_margin = 5mm, right_margin = 2mm,
            plot_title = "CSA $csaname Board-$board with Germanium detector")

    if saveplot
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)Waveforms/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        wvf_str = replace(replace(replace("$(wvf_sel)", "[" => ""), ", " => "-"), "]" => "")
        fname = fpath * "Waveforms$(qc_str)_$(wvf_str).png"
        savefig(ptot, fname)
        @info "Saving figure to $fname"
    end 
    display(ptot)
end 



#  load quality cuts (qc)
qc_folder = "$(@__DIR__)$(pd_data.dspfolder)"
qc_file = qc_folder * "QualityCuts.json"
cuts_qc = readlprops(qc_file).wvf_keep



# plot cut waveforms. random selection based on each cut parameter 
nwf = 8
qc_par = :all
wvf_cut =  findall(.!cuts_qc[qc_par])  # find  bad (cut) waveforms
wvf_sel = wvf_cut[rand(1:length(wvf_cut), nwf)]
plot_waveforms(pd_data; qc = true, qc_par = qc_par, saveplot = true, wvf_sel = wvf_sel)



# plot 8 good waveforms. random selection 
nwf = 8
saveplot = true
wvf_notcut = findall(cuts_qc.all) # find  good waveforms
wvf_sel = wvf_notcut[rand(1:length(wvf_notcut), nwf)]
plot_waveforms(pd_data; qc = false, saveplot = saveplot, wvf_sel = wvf_sel)


