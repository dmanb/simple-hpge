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
include("$(@__DIR__)/../../utils/utils_naming.jl")
include("$(@__DIR__)/../../utils/utils_aux.jl")

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         legendtitlefontsize = 16,
         colorbar_titlefontsize = 16,
         legendforegroundcolor = :silver)

# settings 
PulserChargeInj_keV =  2000 # pulser charge injected in keV
Cinj_fF = 500 # capacitance of the ASIC in femto Farrad 
Rf_MOhm = 500 # feedback resistor in Mega Ohm
Temp_K = 77 # temperature in Kelvin 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# read data 
folder = benchtest_filename_raw(pd_data, Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
wvfs, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1)
p1 = plot(wvfs[1].time, wvfs[1].signal, label = "Raw", 
    xlabel = "Time", ylabel = "Amplitude (a.u.)", linewidth = 2.5, color = :dodgerblue, legend = :right)

# dp pz for example pllot  
bl_stats = signalstats.(wvfs, leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))
wvfs_bl = shift_waveform.(wvfs, -bl_stats.mean)
tail_stats = tailstats.(wvfs_bl, leftendpoint(dsp_config.tail_window), rightendpoint(dsp_config.tail_window))
deconv_flt = InvCRFilter(median(tail_stats.τ))
wvfs_pz = deconv_flt.(wvfs_bl)#[deconv_flt[i](wvfs[i]) for i in eachindex(wvfs)]
plot!(p1, wvfs_pz[1].time, wvfs_pz[1].signal, label = "Pole-zero corrected", legend_title = "Waveform")

 # dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
dsp_par = simple_dsp(wvfs, dsp_config)
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
ENC_eV = csa_voltage2charge.(Cinj_fF, ENC)

# plot result: ENC graph 
grid_rt_trap = ustrip.(collect(dsp_config.e_grid_rt_trap))
p2 = heatmap(grid_rt_trap, ustrip.(grid_ft_trap), ENC_eV,  
    ylabel = "Flat-top time (µs)", 
    xlabel = "Shaping time (µs)",
    zformatter = :plain,
    yformatter = x -> @sprintf("%.1f", x),
     yticks = (ustrip.(grid_ft_trap), [@sprintf("%.1f", x) for x in ustrip.(grid_ft_trap)]), 
    colorbar_title = "\n ENC noise (eV)",
    size = (650, 400),
    right_margin = 30mm,
    left_margin = 5mm,
    bottom_margin = 5mm)

# find minimum ENC value: flat top time:
ENC_min,  ENC_idx = findmin(ENC)

p3 = plot(grid_rt_trap,  ENC_eV[ENC_idx[1],:], 
    xlabel="Shaping time (µs)", 
    linewidth = 2.5, color = :red2,
    ylabel= "\n ENC noise (eV)",
    label = "Flattop time = $(grid_ft_trap[ENC_idx[1]])", 
    legendtitle = "$(length(wvfs)) waveforms",
    legend=:right, 
    yformatter = :plain)
# vline!([ustrip(trap_rt)], label = "Def. shaping time = $trap_rt", color = :grey, linewidth = 2.5, linestyle = :dashdot)

plot(p1, p2,p3,  layout = (3,1), size = (900, 1500), 
    right_margin = 15mm, 
    left_margin = 5mm,
    dpi = 300,
    thickness_scaling = 1.3)

fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
if !ispath(fpath)
    mkpath(fpath)
end
fname = fpath * "TrapFilterOpt_" * benchtest_filename_plot(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K)
savefig(fname)
