relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
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

include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

# plotting options!
saveplot = true
PulserChargesInj_mV =  [5, 10, 25, 50, 100, collect(200:200:2000)...] 

# get data configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)" # get data folder 

### define DATA  loader funciton 
function load_risetime(set, PulserChargesInj_mV)
    # read dsp files and extract relevant parameters
    e_trap_V = Vector{Float64}(undef, length(PulserChargesInj_mV)) # Vout [V]
    risetime = Vector{Quantity}(undef, length(PulserChargesInj_mV))
    pulser_trap_V = Vector{Float64}(undef, length(PulserChargesInj_mV)) # Vin [V]

    for (i, PulserChargeInj_mV) in enumerate(PulserChargesInj_mV)
        dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, PulserChargeInj_mV, set.Temp_K; PulserChargeInj_Mode = "mV")
        if isfile(dsp_file)
            dsp_par = readlprops(dsp_file)
            e_trap_V[i] = dsp_par.e_trap[1]
            risetime[i] = dsp_par.t90[1] .- dsp_par.t10[1]
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

    # convert some units for plotting 
    risetime = uconvert.(u"ns",risetime) # convert to ns
    pulser_trap_MeV =  1e-6 .*csa_voltage2charge.(set.Cinj_fF, pulser_trap_V) .* EnergyPerEHpair_Ge_eV(77)

    return risetime, pulser_trap_MeV
end 

# settings 
set_cold = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
            Rf_MOhm = 500,  # feedback resistor in Mega Ohm
            Temp_K = 77)  # temperature in Kelvin)
            set_warm = merge(set_cold, (Temp_K = 300,))

if set_cold.Cinj_fF == 500
    PulserChargesInj_mV_warm =  [10, 30, 50, 70, 100, 130, 150, 170, 200, 230, 250, 280, 300, 320] 
    PulserChargesInj_mV_cold =  [5, 10, 25, 50, 100, 200, 400] 
elseif set_cold.Cinj_fF == 3000
    PulserChargesInj_mV_warm  =  [2, 5, 7, collect(10:5:50)...] 
    PulserChargesInj_mV_cold =  PulserChargesInj_mV_warm
end
rt_cold, pulser_trap_MeV_cold = load_risetime(set_cold, PulserChargesInj_mV_cold)
rt_warm, pulser_trap_MeV_warm = load_risetime(set_warm, PulserChargesInj_mV_warm)

# plot 
maxCharge_MeV = NaN

(csaname, board) = CSAname(pd_data)
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 10,
    legendforegroundcolor = :silver)

begin 
    p1 = scatter(pulser_trap_MeV_cold, rt_cold, 
        ylabel  = "Risetime: " * L"$t_{90}$" * "- " * L"$ t_{10}$", xlabel = "Injected pulse (MeV)",
        title = "CSA $csaname Board-$board\n C = $(set_cold.Cinj_fF) fF, R = $(set_cold.Rf_MOhm) MOhm",
        label = "T = $(set_cold.Temp_K) K",
        markerstrokewidth = 0,
        legend = :topright,
        markersize = 4,
        color = :dodgerblue)
    scatter!(pulser_trap_MeV_warm, rt_warm, 
        label = "T = $(set_warm.Temp_K) K",
        markerstrokewidth = 0,
        markersize = 4,
        color = :red2,
        alpha = 1,
        dpi = 300)

    if !isnan(maxCharge_MeV)
        xlims!(0, maxCharge_MeV)
    end

    if saveplot
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)Risetime/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "RiseTime_vs_Vin_Cinj$(set_cold.Cinj_fF)fF_Rf$(set_cold.Rf_MOhm)M.png"
        if !isnan(maxCharge_MeV)
            fname = replace(fname, ".png" => "_max$(maxCharge_MeV)MeV.png")
        end
        savefig(fname)
        @info "Saving figure to $fname"
    end 
    display(p1)
end 