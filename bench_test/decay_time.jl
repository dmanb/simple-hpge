import Pkg

Pkg.activate("$(@__DIR__)/../")

using TypedTables, Dates
using Unitful, Formatting, LaTeXStrings, Measurements, Measures, RadiationDetectorDSP, HDF5
using ArraysOfArrays
using Unitful: uconvert
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDataManagement: readlprops, writelprops
using Plots
using CSV
using RadiationDetectorSignals

@__DIR__
include("$(@__DIR__)/util_funcs.jl")


wvfs, decays = read_csv("../../ASIC_data")



path_config = "$(@__DIR__)/../config/dsp_config.json"
# get DSP configuration data --> Can be modified in .json filepath = " $(@__DIR__)/../config/dsp_config.json"
dsp_meta = readlprops(path_config)
dsp_config = DSPConfig(dsp_meta.default)
dsp_config.bl_window
dsp_config.tail_window

# extract decay times of all waveforms with a simple DSP
decay_times = dsp_decay_times(wvfarray, dsp_config)



# get configuration for decay time extraction
pz_config = dsp_meta.pz.default
min_τ, max_τ = pz_config.min_tau, pz_config.max_tau
nbins = pz_config.nbins
rel_cut_fit = pz_config.rel_cut_fit

# plot decay time distribution
histogram(decay_times[min_τ .< decay_times .< max_τ], bins=:fd, label="Decay Time Distribution", xlabel="Decay Time")
savefig(joinpath(figures_folder, "decay_times.png"))

# at first define cut around peak top to get a better fit
cuts_τ = cut_single_peak(decay_times, min_τ, max_τ,; n_bins=nbins, relative_cut=rel_cut_fit)

# fit the decay time distribution at the peak top with a truncated gaussian
result, report = fit_single_trunc_gauss(decay_times, cuts_τ)

# plot resulting ditribution
plot(report, decay_times, cuts_τ, xlabel="Decay Time [µs]", title="Decay Time Distribution")
savefig(joinpath(figures_folder, "decay_time_fit.png"))

# add pars to pars database while overwriting existing pars
pars_db.pz.fit = result
pars_db.pz.τ = result.μ

@info "Found decay time: $(round(u"µs", result.μ, digits=3))"
# write pars to file
writelprops(pars_file, pars_db)