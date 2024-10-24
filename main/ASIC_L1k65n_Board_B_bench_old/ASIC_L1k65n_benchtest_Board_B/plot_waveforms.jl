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
include("$(@__DIR__)/../../utils/utils_naming.jl")

# settings 
saveplot = true 
set = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
    Rf_MOhm = 500,  # feedback resistor in Mega Ohm
    Temp_K = 77,  # temperature in Kelvin 
    PulserChargeInj_keV = 6) # pulser charge injected in keV

# get data configuration (where to and save stuff)
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

# plot function  
function plot_benchtest_waveforms(pd_data, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_keV)}; saveplot::Bool = true)
    # read data 
    folder = benchtest_filename_raw(pd_data, set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K) 
    if !ispath(get(ENV, "ASIC_DATA",".") * folder)
        println("Data folder does not exist: $(get(ENV, "ASIC_DATA",".") * folder)")
    end
    nwvf = 6
    wvfs_ch1, wvfs_ch2, MetaData = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 2, nwvfmax = nwvf*2) # make sure you can read all waveforms 

    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 10,
    legendforegroundcolor = :silver)

    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, nwvf)
    for i = 1:nwvf
        p[i] = plot(wvfs_ch1[i].time, wvfs_ch1[i].signal, label = MetaData["Ch1"]["Channel"])
        plot!(wvfs_ch2[i].time, wvfs_ch2[i].signal, label = MetaData["Ch2"]["Channel"])
        plot!(xlabel = "Time ($(unit(collect(wvfs_ch1[1].time)[1])))", ylabel = "Voltage (V)", 
        legend = :right) # legendtitle =  "Waveform $(i)", 
        if i >1
            plot!(legend = false)
        end
    end 
    ptot = plot(p..., layout = (ceil(Int,nwvf/2), 2), 
            size = (1000, 800), 
            left_margin = 5mm, right_margin = 2mm)

    if saveplot
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)Waveforms/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "Waveforms_" * benchtest_filename_plot(set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K)
        savefig(ptot, fname)
        @info "Saving figure to $fname"
    end 
    display(ptot)
end 

# plot for 1 setting 
plot_benchtest_waveforms(pd_data, set; saveplot = saveplot)

# plot for all settings 
PulserChargesInj_keV =  [1000, 1500, 2500, 3000, 5000, 10000,15, 30, 100, 300] 
Status = fill(false, length(PulserChargesInj_keV))
for i in eachindex(PulserChargesInj_keV)
    @info "Plotting PulserChargeInj_keV = $(PulserChargesInj_keV[i]) ($i/$(length(PulserChargesInj_keV)))"
    set = merge(set, (PulserChargeInj_keV = PulserChargesInj_keV[i],) )
    plot_benchtest_waveforms(pd_data, set; saveplot = saveplot)
end