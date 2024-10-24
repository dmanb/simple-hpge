import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDSP
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendDSP: get_fltpars
using RadiationDetectorDSP
using IntervalSets, TypedTables, StatsBase, PropDicts
using Measures 
using Statistics
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../src/simple_dsp.jl")

PulserVoltage = 0.5
recompute = false 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

dsp_path = "$(@__DIR__)" * pd_data.dspfolder
if !ispath(dsp_path)
    mkpath(dsp_path)
end
dsp_file = dsp_path * "dsp_pars_Pulser$(PulserVoltage)V.json"

if isfile(dsp_file) && !recompute
    dsp_par = readlprops(dsp_file)
    println("Reading DSP results from file: $dsp_file")
else
    # some plotting defaults 
    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
            xguidefontsize = 16, xtickfontsize = 12,
            yguidefontsize = 16, ytickfontsize = 12,
            legendfontsize = 14, titlefontsize = 14,
            legendforegroundcolor = :silver)

    # read data 
    folder = pd_data.datafolder * @sprintf("ASIC_PulserVoltage_0p%.0fV",PulserVoltage*1000)
    wvfs = read_folder_csv(folder; heading = 2)

    if PulserVoltage == 0.4
        wvfs = wvfs[1:end-1]
    end

    # dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
    dsp_par = simple_dsp(wvfs, dsp_config)

    # create a PropDict save results as json file
    dsp_pd = PropDict(Dict(varname => getproperty(dsp_par, varname) for varname in columnnames(dsp_par)))
    writelprops(dsp_file, dsp_pd)
    println("Writing DSP results to file: $dsp_file")
end 

#########  look at results #########
vars = keys(dsp_par) # overview of variables in dsp_par

# decay time
p1 = stephist(dsp_par.tail_τ, bins=20, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))

# rise time 
p2 = stephist(1e3 .* ustrip.(dsp_par.t90 .- dsp_par.t10), 
    bins=20, xlabel= "Rise time (ns)", ylabel="Occurrence", 
    xformatter = x -> @sprintf("%.0f", x),
    label = L"$t_{90} - t_{10}$", 
    fill = true, color = :red2, ylims = (0, :auto))

# baseline std
p3 = stephist(1e3.*filter(x-> x<= quantile(dsp_par.blsigma, 0.95), dsp_par.blsigma), bins=20, xlabel= "Baseline σ (mV)", ylabel="Occurrence", 
label = false, fill = true, color = :orange, ylims = (0, :auto))

# energy with trapezoidal filter 
p4 = stephist(1e3.*dsp_par.e_trap, bins=20, 
    xlabel= "Amplitude (mV)", ylabel="Occurrence", 
    label = @sprintf("Trap filter \nσ = %.1f mV",std(1e3.*dsp_par.e_trap)), 
    fill = true, 
    color = :dodgerblue, 
    ylims = (0, :auto))

# energy with cusp filter
p5 = stephist(1e3.*dsp_par.e_cusp, bins=20, 
    xlabel= "Amplitude (mV)", ylabel="Occurrence", 
    label = @sprintf("Cusp filter \nσ = %.1f mV",std(1e3.*dsp_par.e_cusp)), fill = true, 
    color = :red2, ylims = (0, :auto))

# overview 
plot(p1, p2, p3, p4, layout = (4,1), size = (700, 1850), 
    left_margins = 7mm, thickness_scaling = 1.3,
    plot_title =  "ASIC bench test: Pulser voltage = $PulserVoltage V")

# save figure 
fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
mkpath(fpath)
fname = fpath * "benchtest_dsp_Pulser$(PulserVoltage)V.png"
savefig(fname)

nSamples = 5:5:length(dsp_par.e_trap)
stds = [std(dsp_par.e_trap[1:n]) for n in nSamples]
p5 = plot(nSamples, 1e3*stds, xlabel = "Number of samples", ylabel = "Standard deviation(mV)", linewidth = 2, label = "Trap filter", color = :dodgerblue, dpi = 250)
vline!([50], linestyle = :dash, label = "50 samples", title =  "ASIC bench test: Pulser voltage = $PulserVoltage V")
fname = fpath * "benchtest_dsp_Pulser$(PulserVoltage)V_samples.png"
savefig(fname)


