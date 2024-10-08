import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using Measures 
using Statistics
using ArraysOfArrays
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../src/simple_dsp.jl")

recompute = false
testrun = ""
# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

dsp_path = "$(@__DIR__)" * pd_data.dspfolder
if !ispath(dsp_path)
    mkpath(dsp_path)
end
dsp_file = dsp_path * "dsp_pars_GFET_testrun$testrun.json"

if isfile(dsp_file) && !recompute
    dsp_par = readlprops(dsp_file)
    println("Reading DSP results from file: $dsp_file")
else
    # read data 
    filename = "GFET_detector_run_test$testrun.hdf5"
    wvfs = read_wvfs_h5(filename; folder = pd_data.datafolder, fixtime = true)

    # quality cuts based decay times
    # commentL decay time is set to zero in case tail window of shited waveform has negative values --> remove those with undershoot
    decay_times = dsp_decay_times(wvfs, dsp_config)
    goodIdx  = findall(decay_times .> 0u"µs")
    wvfs = wvfs[goodIdx]

    # dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
    dsp_par = simple_dsp(wvfs, dsp_config)

  # create a PropDict save results as json file
   dsp_pd = PropDict(Dict(varname => getproperty(dsp_par, varname) for varname in columnnames(dsp_par)))
   writelprops(dsp_file, dsp_pd)
   println("Writing DSP results to file: $dsp_file")
end

# # apply some quality cuts
# dsp_par_qc = filter(x -> isfinite(x.e_trap), dsp_par)
# dsp_par_qc = filter(x -> x.tail_τ > 0u"µs", dsp_par_qc)
# filter(x -> x.e_trap > 0, dsp_par_qc)

#########  look at results #########
vars = keys(dsp_par) # overview of variables in dsp_par

# decay time
p1 = stephist(filter(x -> x < 25u"µs",dsp_par.tail_τ), bins=20, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))

# rise time 
p2 = stephist(1e3 .* ustrip.(dsp_par.t90 .- dsp_par.t10), 
    bins=20, xlabel= "Rise time (ns)", ylabel="Occurrence", 
    xformatter = x -> @sprintf("%.0f", x),
    label = L"$t_{90} - t_{10}$", 
    fill = true, color = :red2, ylims = (0, :auto))

# baseline std
p3 = stephist(1e3.*filter(x-> x<= quantile(dsp_par.blsigma, 0.95), dsp_par.blsigma), bins=20, xlabel= "Baseline σ (mV)", ylabel="Occurrence", 
label = false, fill = true, color = :orange, ylims = (0, :auto))

# baseline 
p4 = stephist(1e3.*filter(x-> x<= quantile(dsp_par.blmean, 0.99), dsp_par.blmean), bins=20, 
        xlabel= "Baseline (mV)", ylabel="Occurrence", 
        label = false, fill = true, color = :forestgreen, ylims = (0, :auto))

# energy with trapezoidal filter 
p5 = stephist(1e3.*filter(x-> isa(x, Float64), dsp_par.e_trap), bins=500, 
    xlabel= "Amplitude (mV)", ylabel="Occurrence", 
    label = @sprintf("Trap filter"),# \nσ = %.1f mV",std(1e3.*filter(x-> isa(x, Float64), dsp_par.e_trap))), 
    fill = true, 
    color = :dodgerblue, 
    ylims = (0, :auto))
# savefig(ppath * "GFET_testrun$(testrun)_etrap_uncal.png")

# energy with cusp filter
p6 = stephist(1e3.*filter(x-> isa(x, Float64), dsp_par.e_cusp), bins=20, 
    xlabel= "Amplitude (mV)", ylabel="Occurrence", 
    label = @sprintf("Cusp filter \nσ = %.1f mV",std(1e3.*filter(x-> isa(x, Float64), dsp_par.e_cusp))), fill = true, 
    color = :red2, ylims = (0, :auto))

# overview 
plot(p1, p2, p3, p4, p5,  layout = (5,1), size = (700, 1850), 
    left_margins = 7mm, thickness_scaling = 1.3,
    plot_title =  "GFET test run $(testrun)")
# save figure 
ppath = "$(@__DIR__)" .* pd_data.figurefolder 
mkpath(ppath)
# savefig(ppath * "GFET_testrun$(testrun)_dsp.png")