import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP
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
# filename = "ASIC_detector_run_20240920_194820"
# filename = "ASIC_detector_run_20240923_133751"
# filename = "ASIC_detector_run_20240923_165944"
filename = "ASIC_detector_run_20240923_165944"
wvfs, metadata = read_wvfs_h5(filename; folder = pd_data.datafolder, fixtime = true)
plot(wvfs[1])


#  # open filesignal
#  fname  = get(ENV, "ASIC_DATA",".") * folder * filename * ".hdf5"
#  fopen = h5open(fname,"r")

#  # read metadata
#  metadata = Dict{String, Any}()
#  for k in keys(attributes(fopen))
#      metadata[string(k)] = read(attributes(fopen), k)
#  end

#  uparse(timestep_unit)
#  timestep = timestep *  uparse(timestep_unit)

#  times = [0u"µs":timestep:(length(voltages[i]) - 1)* timestep for i in eachindex(voltages)] 


 


# shift waveforms 
bl_stats = signalstats.(wvfs, leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))
wvfs = shift_waveform.(wvfs, -bl_stats.mean)

decay_times = dsp_decay_times(wvfs, dsp_config)

default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
xguidefontsize = 16, xtickfontsize = 12,
yguidefontsize = 16, ytickfontsize = 12,
legendfontsize = 14, titlefontsize = 14,
legendforegroundcolor = :silver)
pwvfs = Vector{Plots.Plot{Plots.GRBackend}}(undef, 5)
for i in 1:5#10
    pwvfs[i] =  vspan([55u"µs", 70u"µs"], linecolor = :transparent, fillcolor = :silver, label = "tail window")
    if decay_times[i] != 0u"µs"
        plot!(wvfs[i].time, 1e3.*wvfs[i].signal, label=@sprintf("Waveform %d: τ = %.1f µs", i, ustrip(decay_times[i])), color = :dodgerblue)
    else
        plot!(wvfs[i].time, 1e3.*wvfs[i].signal, label="Waveform $i", color = :red2)
    end
    hline!([0.0], color = :black, linewidth = 2, label = false, linestyle = :dash)
end

plot(pwvfs..., 
    layout=(5,2), 
    size=(1200, 1400), 
    left_margin = 5mm,
    xlabel = "Time (µs)", 
    ylabel = "Voltage (mV)")

# ppath = "$(@__DIR__)" .* pd_data.figurefolder 
# mkpath(ppath)
# savefig(ppath * "GFET_testrun$(testrun)_waveforms.png")

