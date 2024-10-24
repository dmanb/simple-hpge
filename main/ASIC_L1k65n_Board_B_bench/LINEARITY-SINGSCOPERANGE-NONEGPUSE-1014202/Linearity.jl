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
plotFlag = true

# settings 
set = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
            Rf_MOhm = 500,  # feedback resistor in Mega Ohm
            Temp_K = 300)  # temperature in Kelvin)
            PulserChargesInj_mV =  [10, 30, 50, 70, 100, 130, 150, 170, 200, 230, 250, 280, 300, 320]

# get data configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)" # get data folder 


# read e-trap from dsp file 
e_trap_V = Vector{Float64}(undef, length(PulserChargesInj_mV))
pulser_trap_V = Vector{Float64}(undef, length(PulserChargesInj_mV))
for (i, PulserChargeInj_mV) in enumerate(PulserChargesInj_mV)
    set = merge(set, (PulserChargeInj_mV = PulserChargesInj_mV[i],) )
    dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, PulserChargeInj_mV, set.Temp_K; PulserChargeInj_Mode = "mV")
    if isfile(dsp_file)
        dsp_par = readlprops(dsp_file)
        e_trap_V[i] = dsp_par.e_trap[1]
        # e_trap_eV[i] = csa_voltage2charge.(set.Cinj_fF, e_trap_V[i])
        println("Reading DSP results from file: $dsp_file")
    else 
        println("File does not exist: $dsp_file")
    end

    pulser_file = replace(dsp_file, ".json" => "_pulser.json")
    if isfile(pulser_file)
        pulser_par = readlprops(pulser_file)
        pulser_trap_V[i] = pulser_par.e_trap[1]
        println("Reading pulser DSP results from file: $pulser_file")
    else 
        println("File does not exist: $pulser_file")
    end
end

plot(e_trap_V, pulser_trap_V)

maxCharge_keV = NaN#4000
VinMode = "eV" # "V" or "e" or "eV"
VoutMode = "eV" # "V" or "e" or "eV"
begin 
    p77 = benchtest_linearity(e_trap_V, pulser_trap_V, set, pd_data; VinMode = VinMode, VoutMode = VoutMode, maxCharge_keV = maxCharge_keV)
    fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)Linearity/"
    if !ispath(fpath)
        mkpath(fpath)
    end
    fname = fpath * "Linearity_Axes-Vin$(VinMode)-Vout$(VoutMode)_" * benchtest_filename_plot(set.Cinj_fF, set.Rf_MOhm, set.Temp_K)
    if !isnan(maxCharge_keV)
        fname = replace(fname, "Linearity_" => "Linearity_max$(maxCharge_keV/1e3)MeV_")
    end
    savefig(fname)
end 
