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
include("$(@__DIR__)/../../../utils/utils_IO.jl")
include("$(@__DIR__)/../../../utils/utils_naming.jl")

shift_wvf = true 
set = (Cinj_fF = 3000, Temp_K = 77)

# read data 
pd_data = readlprops("$(@__DIR__)/data_config.json")
dataset = noise_dataset_str(set.Cinj_fF, set.Temp_K)
folder = pd_data.datafolder * dataset * "/"
if !ispath(get(ENV, "ASIC_DATA",".") * folder)
    println("Data folder does not exist: $(get(ENV, "ASIC_DATA",".") * folder)")
end
wvfs, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1)
@info "load $(length(wvfs)) waveforms"
if shift_wvf == true
    # get baseline mean, std and slope
    dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)
    bl_window                   = dsp_config.bl_window
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
    extraStr = "_shifted"
else
    extraStr = ""
end

# plot 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
 xguidefontsize = 16, xtickfontsize = 12,
 yguidefontsize = 16, ytickfontsize = 12,
 legendfontsize = 14, titlefontsize = 10,
 legendforegroundcolor = :silver)

nwvfs = length(wvfs)
p = Vector{Plots.Plot{Plots.GRBackend}}(undef, nwvfs)
for i in 1:nwvfs
    p[i] =  plot(wvfs[i].time, 1e3 .* wvfs[i].signal, label = false)
    plot!(xlabel = "Time ($(unit(collect(wvfs[1].time)[1])))", ylabel = "Voltage (mV)")
    if shift_wvf == true
        hline!([0], color = :black, linestyle = :dash, label = false)
    end 
end 
plot(p..., layout = (ceil(Int,nwvfs/3), 3), size = (nwvfs*110, nwvfs*70), 
        left_margin = 5mm, bottom_margin = 3mm, right_margin = 2mm, dpi = 300,
        plot_title = "$(extraStr[2:end]) Waveforms, L1k65n CSA Board-B: C = $(set.Cinj_fF)fF,  T = $(set.Temp_K)K",) 
fpath = "$(@__DIR__)" * pd_data.figurefolder * "Waveforms/"
if !ispath(fpath)
    mkpath(fpath)
end
savefig(fpath * "Waveforms$(extraStr)_$dataset.png")