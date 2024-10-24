# trap filter:
# flat-top time = "gap time"
# shaping time/ rise time = averaging time
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
include("$(@__DIR__)/../../../src/simple_dsp.jl")
include("$(@__DIR__)/../../../utils/utils_naming.jl")
include("$(@__DIR__)/../../../utils/utils_aux.jl")

Cinj_fF = 500
noiseAmpFac = 100 # noise amplification after 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# read data: waveform only 
folder = folder = pd_data.datafolder
wvfs_raw, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1);
@sprintf("Reading %d waveforms ", length(wvfs_raw))

# get baseline mean, std and slope
bl_window                   = dsp_config.bl_window
bl_stats = signalstats.(wvfs_raw, leftendpoint(bl_window), rightendpoint(bl_window))

# substract baseline from waveforms
wvfs = shift_waveform.(wvfs_raw, -bl_stats.mean)

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         legendtitlefontsize = 16,
         colorbar_titlefontsize = 16,
         legendforegroundcolor = :silver)

# construct trapezoidal filter 
# get config parameters
e_grid_rt_trap              = (0.5:0.5:5) .*u"µs"#dsp_config.e_grid_rt_trap
e_grid_ft_trap              = (0.5:0.5:5) .*u"µs"#dsp_config.e_grid_ft_trap
enc_pickoff_trap            = 18u"µs"#dsp_config.enc_pickoff_trap
trap_rt, trap_ft = get_fltpars(PropDict(), :trap, dsp_config) 
signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))

# loop over different flattop (ft) and shaping (rt) time combinations and caluclate ENC noise
ENC = NaN.*zeros(Float64, length(e_grid_ft_trap), length(e_grid_rt_trap))
for (i, ft) in enumerate(e_grid_ft_trap)
    try
        ENC[i, :] =  [std(signal_estimator.(TrapezoidalChargeFilter(rt, ft).(wvfs), enc_pickoff_trap)) for rt in e_grid_rt_trap]
    catch e
        ENC[i,:] .= NaN
        @info e
    end
end
ENC[ENC .== 0] .= Inf # if ft > rt, ENC is not calculated and set to 0. remove for better visibility in plot
ENC[isnan.(ENC)] .= Inf 
ENC .= ENC./noiseAmpFac
ENC_e = csa_voltage2charge.(Cinj_fF, ENC)  # convert to electronvolt equivalent



# display results
p1 = plot(wvfs[1].time, 1e3 .* wvfs[1].signal, label = "Raw waveform 1", 
    xlabel = "Time", ylabel = "Amplitude (mV)", linewidth = 2.5, color = :dodgerblue, legend = :topright)

    # plot result: ENC graph 

p2 = heatmap(ustrip.(collect(e_grid_rt_trap)), ustrip.(collect(e_grid_ft_trap)), ENC_e,  
    ylabel = "Flat-top time " * L"$\tau_\textrm{ft}$" * " (µs)", 
    xlabel = "Shaping time " * L"$\tau_\textrm{rt}$" * " (µs)", 
    zformatter = :plain,
    yformatter = x -> @sprintf("%.1f", x),
    yticks = (ustrip.(collect(e_grid_ft_trap)), [@sprintf("%.1f", x) for x in ustrip.(collect(e_grid_ft_trap))]), 
    colorbar_title = "\n ENC noise " * L"$(e^-)$",
    size = (650, 400),
    right_margin = 30mm,
    left_margin = 5mm,
    bottom_margin = 5mm)

# find minimum ENC value: flat top time:
ENC_min,  ENC_idx = findmin(ENC_e)

p3 = plot(ustrip.(collect(e_grid_rt_trap)),  ENC_e[ENC_idx[1],:], linestyle = :dash, color = :silver, linewidth = 2, label = false)
scatter!(ustrip.(collect(e_grid_rt_trap)),  ENC_e[ENC_idx[1],:], 
    xlabel="Shaping time (µs)", 
    linewidth = 0, color = :black,
    markerstrokewidth = 0,
    linestyle = :dash,
    markersize = 6,
    ylabel= "\n ENC noise " * L"$(e^-)$" ,
    label = "Best "*  L"$\tau_\textrm{ft}$" * " = $(collect(e_grid_ft_trap)[ENC_idx[1]])", 
    legendtitle = "$(length(wvfs)) waveforms",
    legend=:topright, 
    yformatter = :plain)

plot(p1, p2,p3,  layout = (3,1), size = (900, 1500), 
    right_margin = 15mm, 
    left_margin = 5mm,
    dpi = 300,
    thickness_scaling = 1.3,
    plot_title = "$(pd_data.datafolder) \n ENC pick-off time $(enc_pickoff_trap)", 
    plot_titlefontsize = 12)


fpath = "$(@__DIR__)$(pd_data.figurefolder)"
if !ispath(fpath)
    mkpath(fpath)
end
fname = replace(replace(fpath * "Noise/ENC_filteropt_pickofftime$(enc_pickoff_trap)_ft$(minimum(e_grid_ft_trap))-$(maximum(e_grid_ft_trap))_rt$(minimum(e_grid_rt_trap))-$(maximum(e_grid_rt_trap))"," " => ""), "." => "p") * ".png"
savefig(fname)
display(plot!())