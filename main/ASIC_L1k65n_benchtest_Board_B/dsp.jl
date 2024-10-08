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
plotFlag = true
PulserChargesInj_keV =  [1000, 1500, 2500, 3000, 5000, 6, 10000,15, 30, 100, 300] # pulser charge injected in keV
# PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 6 ,15, 30, 100, 300]) 
recompute = false 
Cinj_fF = 500 # capacitance of the ASIC in femto Farrad 
Rf_MOhm = 500 # feedback resistor in Mega Ohm
Temp_K = 300#77 # temperature in Kelvin 

Status = fill(false, length(PulserChargesInj_keV))

for (i, PulserChargeInj_keV) in enumerate(PulserChargesInj_keV)
    try 
    # get dsp configuration
    dpath = "$(@__DIR__)/data_config.json"
    pd_data = readlprops(dpath)
    dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

    # get data folder 
    dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
    dsp_file = dsp_folder * benchtest_filename_dsp(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 

    if isfile(dsp_file) && !recompute
        dsp_par = readlprops(dsp_file)
        println("Reading DSP results from file: $dsp_file")
        Status[i] = true
    else
        # read data 
        folder = benchtest_filename_raw(pd_data, Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
        wvfs, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1)

        # dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
        dsp_par = simple_dsp(wvfs, dsp_config)

        # create a PropDict save results as json file
        dsp_pd = PropDict(Dict(varname => getproperty(dsp_par, varname) for varname in columnnames(dsp_par)))
        if !ispath(dsp_folder)
            mkpath(dsp_folder)
        end
        writelprops(dsp_file, dsp_pd)
        println("Writing DSP results to file: $dsp_file")
        Status[i] = true       
    end 
    catch e
        # println("Error: $e")
    end
end 
print(Status)

if plotFlag
    for (i, PulserChargeInj_keV) in enumerate(PulserChargesInj_keV)  
        # get data folder 
        dsp_file = dsp_folder * benchtest_filename_dsp(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
        if !isfile(dsp_file) 
            continue
        end
        dsp_par = readlprops(dsp_file)
        println("Reading DSP results from file: $dsp_file")

        # #########  look at results #########
        # some plotting defaults 
        default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
        xguidefontsize = 16, xtickfontsize = 12,
        yguidefontsize = 16, ytickfontsize = 12,
        legendfontsize = 14, titlefontsize = 10,
        legendforegroundcolor = :silver)

        vars = keys(dsp_par) # overview of variables in dsp_par
        nbins = 0.25*length(dsp_par.t90) < 20 ? 20 :  round(Int, 0.25*length(dsp_par.t90)) 

        # decay time
        p1 = stephist(dsp_par.tail_τ, bins= nbins, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))

        # rise time 
        p2 = stephist(1e3 .* ustrip.(dsp_par.t90 .- dsp_par.t10), 
            bins= nbins, xlabel= "Rise time (ns)", ylabel="Occurrence", 
            xformatter = x -> @sprintf("%.0f", x),
            label = L"$t_{90} - t_{10}$", 
            fill = true, color = :red2, ylims = (0, :auto))

        # baseline std
        p3 = stephist(1e3.*filter(x-> x<= quantile(dsp_par.blsigma, 0.95), dsp_par.blsigma), bins=20, xlabel= "Baseline σ (mV)", ylabel="Occurrence", 
        label = false, fill = true, color = :orange, ylims = (0, :auto))

        # energy with trapezoidal filter  in Volts 
        charge_MeV = 1e-6 .* csa_voltage2charge.(Cinj_fF, dsp_par.e_trap)
        p4 = stephist(charge_MeV, bins = nbins, 
            xlabel= "Amplitude (MeV)", ylabel="Occurrence", 
            label = @sprintf("Trap. filter: σ = %.1f keV", 1e3*std(charge_MeV)), 
            fill = true, 
            color = :dodgerblue, 
            ylims = (0, :auto))

        # overview 
        plot(p1, p2, p3, p4, layout = (2,2), size = (1200, 800), 
            left_margins = 2mm , thickness_scaling = 1.4,
            dpi = 300)

        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "DSPpars_" * benchtest_filename_plot(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K)
        savefig(fname)
    end
end 