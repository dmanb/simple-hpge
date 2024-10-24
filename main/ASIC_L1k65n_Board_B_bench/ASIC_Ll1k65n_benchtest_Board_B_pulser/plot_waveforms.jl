relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDSP
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using ArraysOfArrays
using Measures, Unitful 
using Statistics
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")

# plot settings 
saveplot = true 

# get data configuration (where to and save stuff)
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

# plot function  
function plot_benchtest_waveforms(pd_data, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_mV)}; 
            saveplot::Bool = true, showwindows::Bool = false)
    # read data 
    filepath = ENV["ASIC_DATA"] * pd_data.datafolder * "Linearity_Cinj$(set.Cinj_fF)fF/"
    filename = @sprintf("%dk-%dmV-%dfF_000_ALL.csv",set.Temp_K, set.PulserChargeInj_mV, set.Cinj_fF)
    filepath = filepath * filename
    if !isfile(filename)
        println("Data  does not exist: $filename")
    end
    wvfs_ch1, wvfs_ch2, MetaData = read_file_csv_oscilloscope(filepath;  heading = 17, nChannels = 2) 
    
    # get config 
    dsp_config_set = readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config_overwrite.json")
    dsp_config_tmp = readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json")
    dsp_config = DSPConfig(merge(dsp_config_tmp.default, get(dsp_config_set, Symbol("$(set.Cinj_fF)fF"), PropDict())))

    # substract baseline from waveforms
    bl_stats = signalstats(wvfs_ch1[1], leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))
    wvfs_ch1 = shift_waveform(wvfs_ch1[1], -bl_stats.mean)
    bl_stats2 = signalstats(wvfs_ch2[1], leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))
    wvfs_ch2 = shift_waveform(wvfs_ch2[1], -bl_stats2.mean)

    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 10,
    legendforegroundcolor = :silver)

    plt = plot(wvfs_ch1.time, wvfs_ch1.signal, label = "Waveform ($(MetaData["Ch1"]["Channel"]))",
                linewidth = 2.5, color = :dodgerblue)
    plot!(wvfs_ch2.time, wvfs_ch2.signal, label = "Pulser ($(MetaData["Ch2"]["Channel"]))",
            linewidth = 1.5, color = :red2, linestyle = :solid)
    plot!(xlabel = "Time ($(unit(collect(wvfs_ch1.time)[1])))", ylabel = "Voltage (V)", 
        legend = :right) 
   if showwindows == true
     vspan!([ustrip(leftendpoint(dsp_config.bl_window)), ustrip(rightendpoint(dsp_config.bl_window))],
            color = :violet, alpha = 0.3, label = "Baseline window")
    vspan!([ustrip(leftendpoint(dsp_config.tail_window)), ustrip(rightendpoint(dsp_config.tail_window))],
            color = :orange, alpha = 0.3, label = "Baseline window")
   end
    if saveplot
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)Waveforms/$(set.Cinj_fF)fF/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "Waveform_" * replace(replace(filename, ".csv" => ".png"), "_000_ALL" => "")
        savefig(fname)
        @info "Saving figure to $fname"
    end 
    display(plt)
end 


set = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
    Rf_MOhm = 500,  # feedback resistor in Mega Ohm
    Temp_K = 300,  # temperature in Kelvin 
    PulserChargeInj_mV = 10) # pulser charge injected in keV

# plot for 1 setting 
plot_benchtest_waveforms(pd_data, set; saveplot = false, showwindows = true)

# plot for all settings 
if set.Cinj_fF == 500
    if set.Temp_K == 300
        PulserChargesInj_mV =  [10, 30, 50, 70, 100, 130, 150, 170, 200, 230, 250, 280, 300, 320] 
    elseif set.Temp_K == 77
        PulserChargesInj_mV =  [5, 10, 25, 50, 100, 200, 400] 
    end
elseif set.Cinj_fF == 3000
    PulserChargesInj_mV =  [2, 5, 7, collect(10:5:50)...] 
end

Status = fill(false, length(PulserChargesInj_mV))
for i in eachindex(PulserChargesInj_mV)
    @info "Plotting PulserChargeInj_keV = $(PulserChargesInj_mV[i]) ($i/$(length(PulserChargesInj_mV)))"
    set = merge(set, (PulserChargeInj_mV = PulserChargesInj_mV[i],) )
    plot_benchtest_waveforms(pd_data, set; saveplot = saveplot)
end
