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
include("$(@__DIR__)/../../src/simple_dsp.jl")

PulserVoltage = 0.6500

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         colorbar_titlefontsize = 16,
         legendforegroundcolor = :silver)

# read data 
folder = pd_data.datafolder * @sprintf("ASIC_PulserVoltage_0p%.0fV",PulserVoltage*1000)
wvfs = read_folder_csv(folder; heading = 2)
if PulserVoltage == 0.4
    wvfs = wvfs[1:end-1]
end
p1 = plot(wvfs[1].time, wvfs[1].signal, label = "Waveform 1", xlabel = "Time", ylabel = "Amplitude (a.u.)", linewidth = 2.5, color = :dodgerblue)
for i in 2:5
    plot!(wvfs[i].time, wvfs[i].signal, label = "Waveform $i")
end

# get decay times 
decay_times = dsp_decay_times(wvfs, dsp_config)

## running the dsp_trap_rt_optimization function to get the ENC vs. shaping time plot
trap_rt, trap_ft = get_fltpars(PropDict(),:trap, dsp_config) # default rise-time and flattop-time from config
grid_ft_trap = [0.2, 0.5, 1, 1.5, 2].*u"µs"
ENC = NaN.*zeros(Float64, length(grid_ft_trap), length(dsp_config.e_grid_rt_trap))
for (i, ft) in enumerate(grid_ft_trap)
    try
    local enc_trap_grid = dsp_trap_rt_optimization(wvfs, dsp_config, median(decay_times),; ft = ft)
    ENC[i, :] = [std(enc_trap_grid[j,:]) for j in eachindex(dsp_config.e_grid_rt_trap)]
    catch e
        ENC[i,:] .= NaN
    end
end
ENC[ENC .== 0] .= Inf # if ft > rt, ENC is not calculated and set to 0. remove for better visibility in plot

# plot result: ENC graph 
grid_rt_trap = ustrip.(collect(dsp_config.e_grid_rt_trap))
p2 = heatmap(grid_rt_trap, ustrip.(grid_ft_trap), 1e3.*ENC,  
    ylabel = "Flat-top time (µs)", 
    xlabel = "Shaping time (µs)",
    zformatter = :plain,
    yformatter = x -> @sprintf("%.1f", x),
     yticks = (ustrip.(grid_ft_trap), [@sprintf("%.1f", x) for x in ustrip.(grid_ft_trap)]), 
    colorbar_title = "\n ENC noise " * L"$\times \,10^{3}\,(e^-)$",
    size = (650, 400),
    right_margin = 30mm,
    left_margin = 5mm,
    bottom_margin = 5mm)

# find minimum ENC value: flat top time:
ENC_min,  ENC_idx = findmin(ENC)

p3 = plot(grid_rt_trap, 1e3 .* ENC[ENC_idx[1],:], 
    xlabel="Shaping time (µs)", 
    linewidth = 2.5, color = :dodgerblue,
    ylabel= "\n ENC noise " * L"$\times \,10^{3}\,(e^-)$",
    label = "Flattop time = $(grid_ft_trap[ENC_idx[1]])", 
    legendtitle = "ENC noise based on $(length(wvfs)) waveforms",
    legend=:topright, 
    yformatter = :plain)
vline!([ustrip(trap_rt)], label = "Def. shaping time = $trap_rt", color = :grey, linewidth = 2.5, linestyle = :dashdot)

plot(p1, p2,p3,  layout = (3,1), size = (900, 1500), 
    right_margin = 15mm, 
    left_margin = 5mm,
    dpi = 300,
    thickness_scaling = 1.3,
    plot_title = "ASIC bench test: Pulser voltage = $PulserVoltage V")
# save figure 
fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
mkpath(fpath)
fname = fpath * "benchtest_trapopt_Pulser$(PulserVoltage)V.png"
savefig(fname)