import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
# using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops
# using LegendDSP: get_fltpars
# using RadiationDetectorDSP
using TypedTables, StatsBase, PropDicts
using ArraysOfArrays
using Measures 
using Measurements
using Statistics
using Plots
using Printf, LaTeXStrings
using LegendSpecFits

include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../utils/utils_naming.jl")
include("$(@__DIR__)/../../utils/utils_aux.jl")

# settings 
plotFlag = false
Cinj_fF = 500 # capacitance of the ASIC in femto Farrad 
Rf_MOhm = 500 # feedback resistor in Mega Ohm
Temp_K = 300 # temperature in Kelvin 

if Temp_K == 77
    # PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 6, 10000,15, 30, 100, 300]) # pulser charge injected in keV
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 6 ,15, 30, 100, 300]) # pulser charge injected in keV
elseif Temp_K == 300
    # PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 6, 10000,15, 30]) # pulser charge injected in keV
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 6 ,15, 30]) 
end

# get data configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)" # get data folder 
(csaname, board) = CSAname(pd_data)

# read e-trap from dsp file 
e_trap_V = Vector{Vector{Float64}}(undef, length(PulserChargesInj_keV))
e_trap_eV = Vector{Vector{Float64}}(undef, length(PulserChargesInj_keV))
for (i, PulserChargeInj_keV) in enumerate(PulserChargesInj_keV)
    dsp_file = dsp_folder * benchtest_filename_dsp(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
    if isfile(dsp_file)
        dsp_par = readlprops(dsp_file)
        e_trap_V[i] = dsp_par.e_trap
        e_trap_eV[i] = chaStrge_MeV = csa_voltage2charge.(Cinj_fF, e_trap_V[i])
        println("Reading DSP results from file: $dsp_file")
    else 
        println("File does not exist: $dsp_file")
    end
end


# plot defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
xguidefontsize = 16, xtickfontsize = 12,
yguidefontsize = 16, ytickfontsize = 12,
legendfontsize = 14, titlefontsize = 10,
legendforegroundcolor = :silver,)
nbins = 20

# plot histogram of each injected charge 
plt = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(PulserChargesInj_keV)) 
for (i, ChargeInj) in enumerate(PulserChargesInj_keV)
    if PulserChargesInj_keV[i] >= 1200
        factor = 1e-6
        lbl = "MeV"
       xtformat = x -> @sprintf("%.3f",x)
    else
        factor = 1e-3
        lbl = "keV"
        if PulserChargesInj_keV[i] < 10 
            xtformat = x -> @sprintf("%.0f",x)
        else
            xtformat = x -> @sprintf("%.1f",x)
        end 
    end

    plt[i] = stephist(factor .* e_trap_eV[i], bins = nbins, 
    size = (600,400), 
    label = L"$Q_\textrm{inj}$" * @sprintf("= %d %s", 1e3 * factor * ChargeInj, lbl), 
    fill = true, 
    legend = :topleft,
    alpha = 1,
    color = :silver,  
    xformatter = xtformat, 
    xlabel= "Amplitude ($lbl)", ylabel="Occurrence",
    ylims = (0, :auto))
end

nwvf = [length(e_trap_eV[i]) for i in eachindex(PulserChargesInj_keV)]
if all(nwvf .== nwvf[1])
    println("All waveforms have the same length: $(nwvf[1])")
    wvf_str = @sprintf("%d waveforms ", nwvf[1])
else
    println("Waveforms have different lengths: $(nwvf)")
    wvf_str = @sprintf("%d - %d waveforms ", minimum(nwvf), maximum(nwvf))
end

pall = plot(plt..., layout = (6, 2), 
     left_margin = 15mm, size = (900, 1700), 
     thickness_scaling = 0.9,
     dpi = 300,
     plot_title = "Bench test: $(csaname)-Board $(board) \n \n C = $Cinj_fF fF, R = $Rf_MOhm MOhm, T = $Temp_K K,  $wvf_str each",)

# # save figure 
fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
if !ispath(fpath)
    mkpath(fpath)
end
fname = fpath * "LinearityHists_" * benchtest_filename_plot(Cinj_fF, Rf_MOhm, Temp_K)
savefig(fname)


# plot 2: Charge injected vs. amplitude
medians_keV  = 1e-3 .* [median(e_trap_eV[i]) for i in eachindex(PulserChargesInj_keV)]
medians_err_keV  = 1e-3 .* [std(e_trap_eV[i]) for i in eachindex(PulserChargesInj_keV)] ./ sqrt.(nwvf)
result, report = chi2fit(1, 1e-3 .* PulserChargesInj_keV, 1e-3 .* (medians_keV .± medians_err_keV))

# plot(report, dpi = 300, guidefontsize = 9, xguidefontsize = 9, yguidefontsize = 9, legendfontsize = 9, 
            # titlefontsize = 9, left_margin = 2mm, size = (600, 400))
plot(1e-3 .* PulserChargesInj_keV, Measurements.value.(report.f_fit(1e-3 .* PulserChargesInj_keV)), 
    color = :darkgrey, linewidth = 2,
    label = @sprintf("Fit: (%.2f ± %.0e) ", result.par[2].val, result.par[2].err) * L"$\cdot\,  Q_\textrm{inj}\,$" * @sprintf("+ %.0f keV", 1e3*result.par[1].val))

# plot(1e-3 .* PulserChargesInj_keV, 1e-3 .* PulserChargesInj_keV, linestyle = :dash, color = :black, label = "Linear", legend = :topleft)
scatter!(1e-3 .* PulserChargesInj_keV, 1e-3 .* medians_keV, yerr = 1e-3 .* medians_err_keV,
    size = (700, 450),
    markerstrokewidth=0,
    markersize = 5,
    top_margin = 2mm,
    linewidth = 1,
    label = "Bench test: $(csaname)-Board $(board)" ,
    xlabel = "Injected charge " * L"$Q_\textrm{inj}\,$" * "(MeV)", 
    ylabel = "Median pulse amplitude (MeV)",
    title = "CSA Linearity\n C = $Cinj_fF fF, R = $Rf_MOhm MOhm, T = $Temp_K K,  $wvf_str each", titlefontsize = 12)
fname = fpath * "Linearity3MeV_" * benchtest_filename_plot(Cinj_fF, Rf_MOhm, Temp_K)
savefig(fname)
    

# report.f_fit(1e-3 .* PulserChargesInj_keV)