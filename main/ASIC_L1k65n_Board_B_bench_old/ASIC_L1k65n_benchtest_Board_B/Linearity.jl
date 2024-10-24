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
using Measurements: value as mvalue
using Measurements: uncertainty as muncert
using Statistics
using Plots
using Printf, LaTeXStrings
using LegendSpecFits

include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../utils/utils_naming.jl")
include("$(@__DIR__)/../../utils/utils_aux.jl")
include("$(@__DIR__)/../../bench_test/benchtest_linearity.jl")

# plotting options!
VinMode = "V" # "eV" or "V"  or "e" 
VoutMode = "V" # "eV" or "V"  or "e" 
saveplot = true
maxCharge_keV = NaN
# settings 

set = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
            Rf_MOhm = 500,  # feedback resistor in Mega Ohm
            Temp_K = 300)  # temperature in Kelvin)

if set.Temp_K == 77
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 6, 10000,15, 30, 100, 300]) # pulser charge injected in keV
elseif set.Temp_K == 300
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 6, 10000,15, 30]) # pulser charge injected in keV
end

### LOAD DATA 
# get data configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)" # get data folder 
(csaname, board) = CSAname(pd_data)

# read e-trap from dsp file 
e_trap_V = Vector{Float64}(undef, length(PulserChargesInj_keV))
e_trap_V_err = Vector{Float64}(undef, length(PulserChargesInj_keV))
pulser_trap_V = Vector{Float64}(undef, length(PulserChargesInj_keV))
nwvf = Vector{Int}(undef, length(PulserChargesInj_keV))
# e_trap_eV = Vector{Vector{Float64}}(undef, length(PulserChargesInj_keV))
for (i, PulserChargeInj_keV) in enumerate(PulserChargesInj_keV)
    dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, PulserChargeInj_keV, set.Temp_K) 
   
    if isfile(dsp_file)
        dsp_par = readlprops(dsp_file)
        e_trap_V[i] = mean(dsp_par.e_trap)
        e_trap_V_err[i] = std(dsp_par.e_trap) ./ sqrt(length(dsp_par.e_trap)) # error of the mean
        # e_trap_eV[i] = chaStrge_MeV = csa_voltage2charge.(set.Cinj_fF, e_trap_V[i])
        nwvf[i] = length(dsp_par.e_trap)
        println("Reading DSP results from file: $dsp_file")
    else 
        println("File does not exist: $dsp_file")
    end

    pulser_file = replace(dsp_file, ".json" => "_pulser.json")
    if isfile(pulser_file)
        pulser_par = readlprops(pulser_file)
        pulser_trap_V[i] = mean(pulser_par.e_trap)
        println("Reading pulser DSP results from file: $pulser_file")
    else 
        println("File does not exist: $pulser_file")
    end
end

maxCharge_keV = 5000#NaN
benchtest_linearity(e_trap_V, pulser_trap_V, set, pd_data; VinMode = "eV", VoutMode = "eV", 
                    maxCharge_keV = maxCharge_keV)

if saveplot 
    # save plot 
    fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)/Linearity/"
    if !ispath(fpath)
        mkpath(fpath)
    end
    fname = fpath * "Linearity_" * "Axes$(VinMode)$(VoutMode)_" * benchtest_filename_plot(set.Cinj_fF, set.Rf_MOhm, set.Temp_K)
    if !isnan(maxCharge_keV)
        fname = replace(fname, "Linearity_" => "Linearity_max$(maxCharge_keV/1e3)MeV_")
    end
    savefig(fname)
    @info "Saved plot to $fname"  
end
