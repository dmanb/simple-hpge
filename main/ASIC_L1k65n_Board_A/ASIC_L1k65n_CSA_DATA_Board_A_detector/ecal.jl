
relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
using Plots
using Printf, LaTeXStrings
using PropDicts
using Unitful, Measurements, Measures
using LegendDataManagement: readlprops
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using Statistics, StatsBase
using RadiationSpectra
using LegendSpecFits

include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

# get configs (data, dsp )= 
pd_data = readlprops("$(@__DIR__)/data_config.json")
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
dsp_file = dsp_folder * "DSPpars.json"
(csa, board) =  CSAname(pd_data)
plt_path = "$(@__DIR__)" *pd_data.figurefolder * "Spectrum/"
ispath(plt_path) || mkpath(plt_path)

# load dsp results 
if isfile(dsp_file)
    dsp_par = readlprops(dsp_file)
    println("Reading DSP results from file: $dsp_file")  
else
    println("DSP results not found") 
end

default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
xguidefontsize = 16, xtickfontsize = 12,
yguidefontsize = 16, ytickfontsize = 12,
legendfontsize = 12, titlefontsize = 10,
legendforegroundcolor = :silver)

# uncalibrated energy
C_fF = 3000
etype = :e_trap
e_uncal  = csa_voltage2charge.(C_fF, dsp_par[etype])  .* EnergyPerEHpair_Ge_eV(77) * 1e-3

# very simple calibration 
co60_lines = [1173.237, 1332.501].*u"keV" # keV
nbins_def =  length(e_uncal)/10 > 100 ?  round(Int, length(e_uncal)/10) : 100
h_uncal = fit(Histogram, e_uncal, nbins=nbins_def)
h_decon, peakpos = RadiationSpectra.peakfinder(h_uncal, σ=1.0, backgroundRemove=true, threshold=50)
@info "$(length(peakpos)) peaks found at $peakpos"
plot(h_uncal, lims = (0, :auto), label = "Uncalibrated energy", xlabel = "Energy (a.u.)", ylabel = "Counts")

# simple calibration 
result, report = chi2fit(1, sort(peakpos) , ustrip.(co60_lines); uncertainty = true)
plot(report, ylabel = "Literature values (keV)", xlabel = "Co-60 peaks data (a.u.)", size = (600, 400))
plot!(Plots.get_thickness_scaling(report), thickness_scaling = 1.2, label = false, margins = 3mm)
ylims!(minimum(ustrip.(co60_lines)) - 40, maximum(ustrip.(co60_lines)) + 40)
xlims!(minimum(peakpos) - 40, maximum(peakpos) + 40)

# simple calibrated spectrum
e_unit = unit(co60_lines[1])
co60_simplecal =  mvalue(result.par[2]) 
e_simple = e_uncal .* co60_simplecal .*e_unit
stephist(e_simple, nbins = nbins_def, fill = true , ylims = (0, :auto), label = "Simple-calibrated energy", xlabel = "Energy (keV)", ylabel = "Counts")

# extract simple-calibrated peaks: 
window_sizes = [(25.0u"keV",25.0u"keV") for _ in 1:length(peakpos)]
binning_peak_window = 5.0u"keV"
peakhists, peakstats, h_calsimple, bin_widths = get_peakhists_th228(e_simple, co60_lines, window_sizes; binning_peak_window = binning_peak_window, e_unit=e_unit)


# actual calibration
#1. fit simple-calibrated 2 gamma peaks 
fit_func = :gamma_def
result_fit, report_fit = LegendSpecFits.fit_peaks(peakhists, peakstats, co60_lines; calib_type = :th228, fit_func = fill(fit_func, length(co60_lines)), uncertainty = true)

# 2. fit calibration curve : convert-back to non-calibrated energies 
μ_fit =  [result_fit[p].centroid./co60_simplecal for p in co60_lines]
pp_fit = co60_lines
result_calib, report_calib = fit_calibration(1, μ_fit, pp_fit; e_expression="e", uncertainty = true)
plot(report_calib, size = (650, 400), xlabel = "Energy (a.u.)", ylabel = "Literature values (keV)")
plot!(Plots.get_thickness_scaling(report), thickness_scaling = 1.2, label = false, margins = 3mm)
ylims!(minimum(ustrip.(co60_lines)) - 40, maximum(ustrip.(co60_lines)) + 40)
xlims!(minimum(mvalue.(μ_fit)) - 40, maximum(mvalue.(μ_fit)) + 40)
plot!(xticks = (round.(Int,mvalue.(μ_fit)), string.(round.(Int,mvalue.(μ_fit)))),
     yticks = (round.(Int,ustrip.(co60_lines)), string.(round.(Int, ustrip.(co60_lines)))))
plot!(title = "Calibration of Ge detector with Co-60 peaks\n" * @sprintf("Fit: %.3f * Energy[a.u.] %.3f keV", ustrip(mvalue(result_calib.par[2])), ustrip(mvalue(result_calib.par[1]))), titlefontsize = 11)
savefig(plt_path * "Co60-EnergyCalibration_$etype.png")


# plot peak fits with calibrated results in title 
co60_cal = e_uncal -> mvalue(result_calib.par[1]) .+ mvalue(result_calib.par[2]) .* ustrip.(e_uncal)
µ_cal = co60_cal.(μ_fit)
fwhm_cal = [result_fit[co60_lines[i]].fwhm /co60_simplecal .* mvalue(result_calib.par[2]) for i in eachindex(co60_lines)]

p = plot(broadcast(k -> plot(report_fit[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), keys(report_fit))..., layout=(length(report_fit), 1), size=(1000,710*length(report_fit)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11, dpi = 300)
plot!(legend = :topleft)
for i in [1, 2]
    if i==1
        titlestr = "Co-60 calibration source \n"
    else 
        titlestr = ""
    end 
    plot!(p[i*2-1], title = titlestr * @sprintf("µ = %.2f keV, ", ustrip(µ_cal[i])) *
                    @sprintf("FWHM = (%.2f", ustrip(fwhm_cal[i])) * L"$\,\pm\,$" * @sprintf("%.2f) keV", ustrip(muncert(fwhm_cal[i]))))
end
display(p)
savefig(plt_path * "Co60-CalibrationPeakFits_$etype.png")


# calibrated spectrum
e_cal = co60_cal.(e_uncal)
stephist(e_cal, nbins = nbins_def, fill = true, color = :silver, ylims = (0, :auto), 
        title = "Calibrated energy spectrum ($(etype))", label = "Co-60", xlabel = "Energy (keV)", ylabel = "Counts", legend = :topleft)
xlims!(470, 1600)
savefig(plt_path * "Spectrum_QC_cal_$etype.png")
