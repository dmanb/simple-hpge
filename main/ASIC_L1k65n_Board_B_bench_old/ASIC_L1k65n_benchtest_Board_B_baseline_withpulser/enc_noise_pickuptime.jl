# do trap filter to understand ENC noise 
# trap filter:
# flat-top time = "gap time"
# shaping time/ rise time = averaging time
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

Cinj_fF = 500
# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# read dat: 1 waveform only 
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
enc_pickoff_trap = [20, collect(50:50:1950)...] .* u"µs"#dsp_config.enc_pickoff_trap
trap_rt, trap_ft = 1.0u"µs", 0.2u"µs" 
signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))

# get ENC for all pickoff times
ENC =  [std(signal_estimator.(TrapezoidalChargeFilter(trap_rt, trap_ft).(wvfs), t)) for t in enc_pickoff_trap]
# ENC[ENC .== 0] .= Inf # if ft > rt, ENC is not calculated and set to 0. remove for better visibility in plot
# ENC[isnan.(ENC)] .= Inf 
ENC_eV = csa_voltage2charge.(Cinj_fF, ENC) # convert to electronvolt equivalent

plot(enc_pickoff_trap, ENC_eV,
    xlabel="Time", 
    ylabel= "\n ENC noise (eV)",
    label =  L"$\tau_\textrm{ft}$" * " = $(trap_ft), " * L"$\tau_\textrm{rt}$" * " = $(trap_rt)",
    legend= :topright, 
    ylims = (115, 165),
    linewidth = 2,
    yformatter = :plain,
    top_margin = 5mm,
    title = "ASIC L1k65n Board B: bench test baseline only \n ENC noise vs. pickoff time in baseline", )
hline!([minimum(ENC_eV)], label = false, linestyle = :dash, linewidth = 2.5, color = :silver)
hline!([maximum(ENC_eV)], label = false, linestyle = :dash, linewidth = 2.5, color = :silver)
hline!([mean(ENC_eV)], label = @sprintf("mean ENC = %.0f eV",mean(ENC_eV)), linestyle = :solid, linewidth = 2.5, color = :silver)


fpath = "$(@__DIR__)$(pd_data.figurefolder)"
if !ispath(fpath)
    mkpath(fpath)
end
fname = replace(replace(fpath * "ENCcurve_baselineonly_pickofftimes_ft$(trap_ft)_rt$trap_rt"," " => ""), "." => "p") * ".png"
savefig(fname)
display(plot!())