import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP
using LegendSpecFits
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using ArraysOfArrays
using Measures 
using Measurements
using Measurements: value as mvalue
using Measurements: uncertainty as muncert
using Statistics
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../utils/utils_naming.jl")
include("$(@__DIR__)/../../src/simple_dsp.jl")

# settings 
recompute = false 
Cinj_fF = 500 # capacitance of the ASIC in femto Farrad 
Rf_MOhm = 500 # feedback resistor in Mega Ohm
PulserChargeInj_keV = 300 # pulser charge injected in keV
Temp_K = 77 # temperature in Kelvin 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# read data 
folder = benchtest_filename_raw(pd_data, Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
if !ispath(get(ENV, "ASIC_DATA",".") * folder)
    println("Data folder does not exist: $(get(ENV, "ASIC_DATA",".") * folder)")
end
(wvf, MetaData) = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1, nwvfmax = 1)

# do standard simple dsp
dsp_par = simple_dsp(wvf, dsp_config)[1]
wvf_shift = shift_waveform.(wvf, -dsp_par.blmean)

# reproduce fits of decay time (do fit manually)
function exp_func(t, par)
    return par[1] *  exp( - t * par[2])
end
TailWindow  = findall(x -> leftendpoint(dsp_config.tail_window) <= x <= rightendpoint(dsp_config.tail_window), wvf_shift.time[1])
x = collect(wvf_shift.time[1])[TailWindow] .- dsp_par.t99
y = wvf_shift.signal[1][TailWindow]
v_init = [dsp_par.e_trap, 1/ustrip(dsp_par.tail_τ)] # initual guess of fit parameter 
result, report = chi2fit((t, p1, p2) -> exp_func(t, [p1, p2]), 
                        ustrip.(x), y; v_init = v_init)

# plot 
time = collect(wvf_shift.time[1])
signal = wvf_shift.signal[1]
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
 xguidefontsize = 16, xtickfontsize = 12,
 yguidefontsize = 16, ytickfontsize = 12,
 legendfontsize = 14, titlefontsize = 10,
 legendforegroundcolor = :silver)
ppulse = vspan([time[TailWindow][1], time[TailWindow][end]], color = :silver, 
        alpha = 0.5, label = "Tail fit interval", linewidth = 0)
plot!(time, signal, label = "Waveform",
        xlabel = "Time", 
        ylabel = "Voltage (V)", 
        color = :black,
        legend = :topright)

# plot result from dsp decay time and amplitude
tail_dsp = exp_func.(ustrip.(time .- dsp_par.t99), Ref([dsp_par.e_trap, 1 / ustrip(dsp_par.tail_τ)]))
timeIdx = findfirst(x -> x > dsp_par.t99, time)
plot!(time[timeIdx:end], tail_dsp[timeIdx:end], 
    color = :orange, linewidth = 2, linestyle = :dash, 
    label = @sprintf("Exp. decay fit: τ = %.0f µs", ustrip(dsp_par.tail_τ)))

# plot manual fit result 
tail_fit = exp_func.(ustrip(time .- dsp_par.t99), Ref([mvalue(result.par[1]), mvalue(result.par[2])]))
plot!(time[timeIdx:end], tail_fit[timeIdx:end], 
        color = :dodgerblue, label = "Exp. decay fit: manual", linewidth = 2, linestyle = :dashdot)
xl = xlims()

# plot residuals 
PltStep = 100 # plot only every 10 points 
pres = scatter(time[timeIdx:PltStep:end], signal[timeIdx:PltStep:end] .- tail_dsp[timeIdx:PltStep:end],
                color = :orange, marker = :circle,
                markersize = 2, markerstrokewidth = 0,
                xlabel = "Time", 
                ylabel = "Residuals (V)", 
                legend = false,
                yticks = 1e-4.*collect(-10:10:10),
                xlims = xl)
scatter!(time[timeIdx:PltStep:end], signal[timeIdx:PltStep:end] .- tail_fit[timeIdx:PltStep:end],
                color = :dodgerblue, marker = :diamond, markersize = 2, markerstrokewidth = 0)
hline!([0], color = :black, linewidth = 2, linestyle = :dash, label = false)

# combi plot 
layout = @layout [a{0.65h}; b]
(csaname, board) = CSAname(pd_data)
plot(ppulse, pres, layout = layout, size = (600, 600),
    plot_title = "Bench test: $(csaname)-Board $(board) \n" * "C = $Cinj_fF fF, R = $Rf_MOhm MOhm, T = $Temp_K K", plot_titlefontsize = 10)

# # save figure 
fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
if !ispath(fpath)
    mkpath(fpath)
end
fname = fpath * "DecayTimeResiduals_" * benchtest_filename_plot(Cinj_fF, Rf_MOhm, Temp_K)
savefig(fname)
