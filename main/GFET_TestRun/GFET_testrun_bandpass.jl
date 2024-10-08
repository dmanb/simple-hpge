import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP, DSP 
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using Measures 
using Statistics
using ArraysOfArrays
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../src/simple_dsp.jl")

recompute = false
testrun = ""
# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)


# read data 
filename = "GFET_detector_run_test$testrun.hdf5"
wvfs = read_wvfs_h5(filename; folder = pd_data.datafolder, fixtime = true, nwvfs=100)

# plot a waveform 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
xguidefontsize = 16, xtickfontsize = 12,
yguidefontsize = 16, ytickfontsize = 12,
legendfontsize = 14, titlefontsize = 14,
legendforegroundcolor = :silver)

# contruct a digital filter: band stop filter
# rise time is in the order of 100 ns (<--> frequency around 10 MHz)
#   keep frequencies around 1/10 ns --> 1/500 ns: 1e8 Hz to 2e6 Hz
# fall time around 10 µs --> 100 kHz 
sampling_f_Hz = 1/uconvert.(u"s", wvfs[1].time[2])
responsetype = Bandstop(300*1e3, 400*1e3; fs = ustrip(sampling_f_Hz))
designmethod = Butterworth(5)
filter_bandpass = digitalfilter(responsetype, designmethod)

# apply filter and compare 
wvfidx = 5
begin 
    wvfs_filt = filt(filter_bandpass, wvfs[wvfidx].signal)
    # plot result 
    plot(wvfs[wvfidx].time, 1e3.*wvfs[wvfidx].signal, 
        size=(800, 600), 
        color = :dodgerblue,
        label = "Original", 
        xlabel = "Time", 
        ylabel = "Voltage (mV)")
    plot!(wvfs[wvfidx].time, 1e3.*wvfs_filt, label = "Filtered", color = :red2, alpha = 0.2)
end 