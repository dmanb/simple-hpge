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

# settings  
Benchtest_settings = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
                    Rf_MOhm = 500,  # feedback resistor in Mega Ohm
                    Temp_K = 77,  # temperature in Kelvin 
                    PulserChargeInj_keV = 15) # pulser charge injected in keV
saveplot = true 
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

function do_benchtest_encpickofftime(pd_data::PropDict, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_keV)}; saveplot = true)
    dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)
    (csaname, board) = CSAname(pd_data)

    # read dat: 1 waveform only 
    folder = benchtest_filename_raw(pd_data, set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K) 
    wvfs_raw, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1);
    @sprintf("Reading %d waveforms ", length(wvfs_raw))

    # get config parameters
    bl_window                   = dsp_config.bl_window
    enc_pickoff_trap            = dsp_config.enc_pickoff_trap

    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs_raw, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs_raw, -bl_stats.mean)

    # get decay time 
    tail_τ = median(dsp_decay_times(wvfs, dsp_config))

    # deconvolute wavefor: pole-zero correction
    deconv_flt = InvCRFilter(tail_τ)
    wvfs_pz = deconv_flt.(wvfs)

    # construct trapezoidal filter 
    trap_rt, trap_ft = get_fltpars(PropDict(), :trap, dsp_config) 
    signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))
    uflt_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)

    ENC_pickoff_trap = leftendpoint(bl_window):5u"µs":rightendpoint(bl_window)#15u"µs":5u"µs":150u"µs"
    enc_V  = [std(signal_estimator.(uflt_rtft.(wvfs_pz), enc_pickoff_trap)) for enc_pickoff_trap in ENC_pickoff_trap]
    enc_eV = csa_voltage2charge.(set.Cinj_fF, enc_V)
    
    # some plotting defaults 
    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 14,
    legendtitlefontsize = 16,
    colorbar_titlefontsize = 16,
    legendforegroundcolor = :silver)

    # plot waveforms 
    p1 = plot(wvfs[1].time, wvfs[1].signal, label = "Waveform raw",
        xlabel = "Time", ylabel = "Voltage (V)", linewidth = 2.5, color = :dodgerblue, legend = :right)
    plot!(wvfs_pz[1].time, wvfs_pz[1].signal, label = "Waveform pz-corrected", 
        linewidth = 2.5, color = :red2, legend = :right)
    vline!([maximum(ENC_pickoff_trap)], label =  L"Median $\sigma_\textrm{baseline}$" * @sprintf(" = %.2f mV",1e3*median(bl_stats.sigma) ), linestyle = :dash, linewidth = 2.5, color = :silver)

    penc = plot(ENC_pickoff_trap, enc_eV, 
            title = "Trap filter, $(length(wvfs_raw)) wvf: " * L"$\tau_\textrm{rt}$" * " = $(trap_rt), " *  L"$\tau_\textrm{ft}$" * " = $(trap_ft)" ,
            xlabel = "ENC noise pickoff time", 
            ylabel = "ENC (eV)", 
            label = false,
            linewidth = 2.5, color = :dodgerblue, 
            legend = :best,
        )
    enc_eV_flt = filter(x->isfinite(x), enc_eV)
    ENCmin, ENCmax, ENCmean =  [func(enc_eV_flt) for func in [minimum, maximum, mean]]  
    hline!(penc,[ENCmean], label = false, 
            linestyle = :solid, linewidth = 2.5, color = :silver)
    hline!(penc, [ENCmin], label =false, 
            linestyle = :dash, linewidth = 2.5, color = :silver)
    hline!(penc, [ENCmax], label =false, 
            linestyle = :dash, linewidth = 2.5, color = :silver)
    plot!(ylims = (0.9 * ENCmin, 1.1*ENCmax))
    annotate!(0.03 * xlims()[2], ENCmax * 1.05, text(@sprintf("mean = %.0f eV", mean(filter(x->isfinite(x), enc_eV))), 14, :left, :black))
    display(penc)

    if set.PulserChargeInj_keV < 1000
        ustr = "keV"
        PulserChargeInj_lbl = set.PulserChargeInj_keV
    else
        ustr = "MeV"
        PulserChargeInj_lbl = set.PulserChargeInj_keV/1e3
    end

    pall = plot(p1, penc, layout = (2,1), size = (620, 800),
            plot_title = "Bench test: $(csaname)-Board $(board) \n" * "C = $(set.Cinj_fF) fF, R = $(set.Rf_MOhm) MOhm, T = $(set.Temp_K) K, " * L"$Q_\textrm{inj}$" * " = $(PulserChargeInj_lbl) $ustr",
            plot_titlefontsize = 12,
            left_margin = 5mm)  
            
    if saveplot 
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)ENCnoise/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "ENC_pickofftime_" * "$(length(wvfs_raw))wvfs_" * benchtest_filename_plot(set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K) 
        savefig(pall, fname)
        print("Save figure to $fname")
    end
    display(pall)
    return pall, penc
end 

# do for 1 setting 
pall, penc = do_benchtest_encpickofftime(pd_data, Benchtest_settings; saveplot = saveplot)

# do ENC pickoff time study for all different pulser charges
PulserChargesInj_keV =  [1000, 1500, 2500, 3000, 5000, 10000,15, 30, 100, 300] 
for i in eachindex(PulserChargesInj_keV)
    @info "Processing PulserChargeInj_keV = $(PulserChargesInj_keV[i]) ($i/$(length(PulserChargesInj_keV)))"
    Benchtest_settings = merge(Benchtest_settings, (PulserChargeInj_keV = PulserChargesInj_keV[i],) )
    do_benchtest_encpickofftime(pd_data, Benchtest_settings; saveplot = saveplot)
end


# # Other investigations 
# # plot raw baselines 
# bl_window_Idx = findfirst( wvfs[1].time .> rightendpoint(bl_window))
# bl_time = wvfs.time[1][1:bl_window_Idx]
# bl_voltage = [wvfs.signal[i][1:bl_window_Idx] for i = 1:length(wvfs)]
# medians = [median([bl_voltage[w][i] for w = 1:length(wvfs)]) for i = 1:length(bl_time)]
# stds = [std([bl_voltage[w][i] for w = 1:length(wvfs)]) for i = 1:length(bl_time)]

# p1 = plot(bl_time, bl_voltage[1], alpha = 0.3, color = :red2, label = "Raw baseline 1")
# plot!(bl_time, bl_voltage[2], alpha = 0.3, color = :dodgerblue, label = "Raw baseline 2")
# plot!(bl_time, medians, ribbon = stds, label = "median with 1σ band", color = :black, alpha = 0.7, fillcolor = :darkgrey, fillalpha = 0.8)
# display(p1)
# # plot(bl_time, csa_voltage2charge.(Cinj_fF, stds))

# # compare raw with filtered baseline 
# # trap_rt, trap_ft = get_fltpars(PropDict(), :trap, dsp_config) 
# uflt_rtft_def = TrapezoidalChargeFilter(trap_rt, trap_ft)
# uflt_rtft_rt8 = TrapezoidalChargeFilter(8.0u"µs", trap_ft)
# uflt_rtft_rt2 = TrapezoidalChargeFilter(2.5u"µs", trap_ft)
# bl_time_trap = bl_time[findall(x-> x > trap_rt * 2 + trap_ft * 3  , bl_time )][1:10:end]
# pbl = plot(bl_time, bl_voltage[1], alpha = 0.2, color = :silver, label = "Raw", legend_title = "Baseline")
# plot!(bl_time_trap, [signal_estimator(uflt_rtft_def(wvfs_pz[1]), time) for time in bl_time_trap], linewidth = 2, color = :red2, 
#     label = "Trap. filtered: " * L"$t_\textrm{rt}$" *  @sprintf(" = %.1f %s, ", ustrip(uflt_rtft_def.avgtime), string(unit(uflt_rtft_def.avgtime))) * L"$t_\textrm{ft}$" * @sprintf(" = %.1f %s" , ustrip(uflt_rtft_def.gaptime), string(unit(uflt_rtft_def.gaptime))) )
# plot!(bl_time_trap, [signal_estimator(uflt_rtft_rt8(wvfs_pz[1]), time) for time in bl_time_trap], linewidth = 2, color = :dodgerblue, 
#     label = "Trap. filtered: " * L"$t_\textrm{rt}$" *  @sprintf(" = %.1f %s, ", ustrip(uflt_rtft_rt8.avgtime), string(unit(uflt_rtft_rt8.avgtime))) * L"$t_\textrm{ft}$" * @sprintf(" = %.1f %s" , ustrip(uflt_rtft_def.gaptime), string(unit(uflt_rtft_def.gaptime))) )
# plot!(bl_time_trap, [signal_estimator(uflt_rtft_rt2(wvfs_pz[1]), time) for time in bl_time_trap], linewidth = 2, color = :green3, 
#     label = "Trap. filtered: " * L"$t_\textrm{rt}$" *  @sprintf(" = %.1f %s, ", ustrip(uflt_rtft_rt2.avgtime), string(unit(uflt_rtft_rt2.avgtime))) * L"$t_\textrm{ft}$" * @sprintf(" = %.1f %s" , ustrip(uflt_rtft_def.gaptime), string(unit(uflt_rtft_def.gaptime))) )
# hline!([0], color = :darkgrey, linestyle = :dash, linewidth = 1.5, label = false,
#     legend = :outertop, size = (800, 600), xlabel = "Time (µs)", ylabel = "Voltage (V)")
# # # save figure 
# fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
# fname = fpath * "Baseline_raw_trapfilter_"  * benchtest_filename_plot(set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K) 
# savefig(pbl, fname)
# print("Save figure to $fname")
# display(pbl)



# #### other stuff 
# baseline_trap  = [signal_estimator.(uflt_rtft.(wvfs_pz), enc_pickoff_trap) for enc_pickoff_trap in ENC_pickoff_trap]
# medians = [median(baseline_trap[i]) for i = 1:length(ENC_pickoff_trap)]
# stds = [std(baseline_trap[i]) for i = 1:length(ENC_pickoff_trap)]
# pmedian= plot(ENC_pickoff_trap, medians, ribbon = stds,
#     linewidth = 3, color = :dodgerblue, alpha = 1, fillalpha = 0.3, 
#     label = "Trap. filtered baseline. All waveforms median ± 1σ",
#     legendfontsize = 12,
#     ylims = (-0.0005, 0.0007),
#     xlabel = "Time", ylabel = "Voltage (V)", legend = :topleft,
#     size = (700, 420), top_margin = 7mm, left_margin = 3mm, bottom_margin = 3mm,
#     plot_title = "Bench test: $(csaname)-Board $(board) \n" * "C = $Cinj_fF fF, R = $Rf_MOhm MOhm, T = $Temp_K K, " * L"$Q_\textrm{inj}$" * " = $(PulserChargeInj_keV/1e3) MeV",
#     plot_titlefontsize = 12)

# plot!(fill(ustrip(ENC_pickoff_trap[1]), 2), [medians[1] , medians[1] + stds[1]], 
#     arrow=true, color=:orange, linewidth=2, label="ENC noise different pickoff times", legend = :topleft)

# for aidx = 6:5:length(ENC_pickoff_trap)
#     @info aidx
#     plot!(fill(ustrip(ENC_pickoff_trap[aidx]), 2), [medians[aidx] , medians[aidx] + stds[aidx]], 
#         arrow=true, color=:orange, linewidth=2, label=false, legend = :topleft)
# end
# display(pmedian)
# # # save figure 
# fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
# fname = fpath * "Baseline_trapfilter_median_std_"  * benchtest_filename_plot(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
# savefig(pmedian, fname)
# print("Save figure to $fname")
# display(pmedian)




