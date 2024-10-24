import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDataManagement
using LegendDataManagement: readlprops
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

# plotting options!
saveplot = true
maxCharge_keV = NaN

# settings 
set = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
            Rf_MOhm = 500,  # feedback resistor in Mega Ohm
            Temp_K = 300)  # temperature in Kelvin)

if set.Temp_K == 77
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 10000, 15, 6,  30, 100, 300]) # pulser charge injected in keV
elseif set.Temp_K == 300
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 10000, 15, 6, 30]) # pulser charge injected in keV
end

### LOAD DATA 
# get data configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)" # get data folder 
(csaname, board) = CSAname(pd_data)


# read dsp files and extract relevant parameters
e_trap_V = Vector{Float64}(undef, length(PulserChargesInj_keV)) # Vout [V]
risetime = Vector{Quantity}(undef, length(PulserChargesInj_keV))
risetime_err = Vector{Quantity}(undef, length(PulserChargesInj_keV))
pulser_trap_V = Vector{Float64}(undef, length(PulserChargesInj_keV)) # Vin [V]
nwvf = Vector{Int}(undef, length(PulserChargesInj_keV))
for (i, PulserChargeInj_keV) in enumerate(PulserChargesInj_keV)
    dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, PulserChargeInj_keV, set.Temp_K)  
    if isfile(dsp_file)
        dsp_par = readlprops(dsp_file)
        e_trap_V[i] = median(dsp_par.e_trap)
        risetime_err[i] = std(dsp_par.t90 .- dsp_par.t10) ./ sqrt(length(dsp_par.t90))
        risetime[i] = median(dsp_par.t90 .- dsp_par.t10)
        risetime_err[i] = std(dsp_par.t90 .- dsp_par.t10) ./ sqrt(length(dsp_par.t90))
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

# convert some units for plotting 
risetime = uconvert.(u"ns",risetime) # convert to ns
pulser_trap_MeV =  1e-6 .*csa_voltage2charge.(set.Cinj_fF, pulser_trap_V) .* EnergyPerEHpair_Ge_eV(set.Temp_K)

scatter(pulser_trap_MeV, risetime, 
    ylabel  = "Risetime", xlabel = "Injected pulse (MeV)",
    label =  L"$t_{90} -  t_{10}$")
