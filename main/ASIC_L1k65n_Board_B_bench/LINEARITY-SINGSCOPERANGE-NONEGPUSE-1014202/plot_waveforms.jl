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
    Temp_K = 300,  # temperature in Kelvin 
    PulserChargeInj_mV = 100) # pulser charge injected in keV

# get data configuration (where to and save stuff)
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

# plot function  
function plot_benchtest_waveforms(pd_data, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_mV)}; saveplot::Bool = true)
    # read data 
    filepath = ENV["ASIC_DATA"] * pd_data.datafolder 
    filename = @sprintf("%dk-%dmV-%dfF_000_ALL.csv",set.Temp_K, set.PulserChargeInj_mV, set.Cinj_fF)
    filepath = filepath * filename
    if !isfile(filename)
        println("Data  does not exist: $filename")
    end
    wvfs_ch1, wvfs_ch2, MetaData = read_file_csv_oscilloscope(filepath;  heading = 17, nChannels = 2) 

    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 10,
    legendforegroundcolor = :silver)

    wvfs_ch1 = wvfs_ch1[1]
    plt = plot(wvfs_ch1.time, wvfs_ch1.signal, label = MetaData["Ch1"]["Channel"])
    plot!(wvfs_ch2.time, wvfs_ch2.signal, label = MetaData["Ch2"]["Channel"])
    plot!(xlabel = "Time ($(unit(collect(wvfs_ch1.time)[1])))", ylabel = "Voltage (V)", 
    legend = :right) 
  

    if saveplot
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)Waveforms/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "Waveform_" * replace(filename, ".csv" => ".png")
        savefig(fname)
        @info "Saving figure to $fname"
    end 
    display(plt)
end 

# plot for 1 setting 
plot_benchtest_waveforms(pd_data, set; saveplot = saveplot)


# plot for all settings 
PulserChargesInj_mV =  [10, 30, 50, 70, 100, 130, 150, 170, 200, 230, 250, 280, 300, 320] 
Status = fill(false, length(PulserChargesInj_mV))
for i in eachindex(PulserChargesInj_mV)
    @info "Plotting PulserChargeInj_keV = $(PulserChargesInj_mV[i]) ($i/$(length(PulserChargesInj_mV)))"
    set = merge(set, (PulserChargeInj_mV = PulserChargesInj_mV[i],) )
    plot_benchtest_waveforms(pd_data, set; saveplot = saveplot)
end
