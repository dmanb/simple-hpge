relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
using Unitful
using LegendDataManagement
using LegendDataManagement: readlprops
using TypedTables, StatsBase, PropDicts
using LegendDSP, RadiationDetectorDSP
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
saveplot = false
showResiduals = true 

# get data configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)
dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)" # get data folder 

function exp_func(t, par)
    return par[1] *  exp( - t * par[2])
end

### define DATA  loader funciton 
function load_decaytime(set, PulserChargesInj_mV)
    # read dsp files and extract relevant parameters
    e_trap_V = Vector{Float64}(undef, length(PulserChargesInj_mV)) # Vout [V]
    baseline_V = Vector{Float64}(undef, length(PulserChargesInj_mV)) # Vout [V]
    decaytime = Vector{Quantity}(undef, length(PulserChargesInj_mV))
    pulser_trap_V = Vector{Float64}(undef, length(PulserChargesInj_mV)) # Vin [V]
    t99 = Vector{Quantity}(undef, length(PulserChargesInj_mV)) # Vin [V]
    tail = Vector{Function}(undef, length(PulserChargesInj_mV)) # 
    for (i, PulserChargeInj_mV) in enumerate(PulserChargesInj_mV)
        dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, PulserChargeInj_mV, set.Temp_K; PulserChargeInj_Mode = "mV")
        if isfile(dsp_file)
            dsp_par = readlprops(dsp_file)
            e_trap_V[i] = dsp_par.e_trap[1]
            decaytime[i] = dsp_par.tail_τ[1]
            t99[i] = dsp_par.t99[1]
            baseline_V[i] = dsp_par.blmean[1]
            tail[i] = Base.Fix2(exp_func, [dsp_par.tailoffset[1], 1 / ustrip(dsp_par.tail_τ[1])])
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
    pulser_trap_MeV =  1e-6 .*csa_voltage2charge.(set.Cinj_fF, pulser_trap_V) .* EnergyPerEHpair_Ge_eV(77)
   
    return decaytime, e_trap_V, t99, baseline_V, pulser_trap_MeV, tail
end 


# settings 
set_cold = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
            Rf_MOhm = 500,  # feedback resistor in Mega Ohm
            Temp_K = 77)  # temperature in Kelvin)
set_warm = merge(set_cold, (Temp_K = 300,))
if set_cold.Cinj_fF == 500
    PulserChargesInj_mV =  [10, 50, 100, 200] 
elseif set_cold.Cinj_fF == 3000
    PulserChargesInj_mV =  [2, 5, 7, collect(10:5:50)...] 
end

# load dsp 
τ_cold, e_trap_V_cold, t99_cold, baseline_cold, pulser_trap_MeV_cold, tail_cold = load_decaytime(set_cold, PulserChargesInj_mV)
τ_warm, e_trap_V_warm, t99_warm, baseline_warm, pulser_trap_MeV_warm, tail_warm = load_decaytime(set_warm, PulserChargesInj_mV)

# load waveforms
path = ENV["ASIC_DATA"] * pd_data.datafolder * "Linearity_Cinj$(set_cold.Cinj_fF)fF/"
filenames_cold = [@sprintf("%dk-%dmV-%dfF_000_ALL.csv", set_cold.Temp_K, p, set_cold.Cinj_fF) for p in PulserChargesInj_mV]
filenames_warm = [@sprintf("%dk-%dmV-%dfF_000_ALL.csv", set_warm.Temp_K, p, set_warm.Cinj_fF) for p in PulserChargesInj_mV]
filepath_cold = path .* filenames_cold
filepath_warm = path .* filenames_warm 
wvfs_cold = Vector{}(undef, length(filenames_cold))
for (i, file) in enumerate(filepath_cold)
    local wvfs_ch1, _,  = read_file_csv_oscilloscope(file;  heading = 17, nChannels = 1)#, - baseline_cold[i]
    wvfs_cold[i] = wvfs_ch1[1]
end
wvfs_warm = Vector{}(undef, length(filenames_warm))
for (i, file) in enumerate(filepath_warm)
    local wvfs_ch1, _,  = read_file_csv_oscilloscope(file;  heading = 17, nChannels = 1)
    wvfs_warm[i] = wvfs_ch1[1]
end

# get tail from dsp 


# PLOT 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 10,
    legendforegroundcolor = :silver)



i = 4
begin
    t99diff =   t99_cold[i] - t99_warm[i]
    time = collect(wvfs_cold[i].time)
    tail_dsp_cold = tail_cold[i].(ustrip.(time .- t99_cold[i]))#exp_func.(ustrip.( time .- t99_cold[i]), Ref([e_trap_V_cold[i], 1 / ustrip(τ_cold[i])]))
    timeIdx = findfirst(x -> x > ustrip.(t99_cold .-t99diff)[1] + 10 , ustrip.(time)) 
    time_warm = collect(wvfs_warm[i].time)
    # tail_warm = exp_func.(ustrip.( time_warm .- t99_warm[i]), Ref([e_trap_V_warm[i], 1 / ustrip(τ_warm[i])]))
    tail_dsp_warm = tail_warm[i].(ustrip.(time .- t99_warm[i]))
    timeIdx_warm = findfirst(x -> x > ustrip.(t99_warm)[1]  + 10 , ustrip.(time_warm)) 

    pulserMeV  = mean([pulser_trap_MeV_cold[i], pulser_trap_MeV_warm[i]])
    if pulserMeV < 1 
        pulserStr = @sprintf("%.0f keV", 1e3 * pulserMeV)
        pulsersavestr = @sprintf("%.0fkeV", 1e3 * pulserMeV)
    else
        pulserStr = @sprintf("%.0f MeV", pulserMeV)
        pulsersavestr = @sprintf("%.0fMeV", pulserMeV)
    end
    pulserlbl = "Injected pulse: $pulserStr"

   
    
    plts = plot(collect(wvfs_cold[i].time) .-t99diff, wvfs_cold[i].signal .- baseline_cold[i], linewidth = 2.5, color = :dodgerblue, 
            label =  @sprintf("τ = %.0f µs ; ", ustrip(τ_cold[i])) * "$(set_cold.Temp_K) K" ,
            legendtitle = pulserlbl, legendtitlefontsize = 13)
   # plot!(time[timeIdx:end]  , tail_dsp_cold[timeIdx:end], linewidth = 2, linestyle = :dash, color = :black, label = false )
    vline!([τ_cold[i] .+ t99_cold[i] - t99diff], linewidth = 1.5, linestyle = :dash, color = :dodgerblue, label = false)
    plot!(collect(wvfs_warm[i].time), wvfs_warm[i].signal .- baseline_warm[i], linewidth = 2.5, color = :red2, label =  @sprintf("τ = %.0f µs ; ", ustrip(τ_warm[i])) * "$(set_warm.Temp_K) K" )
    vline!([τ_warm[i] .+ t99_warm[i]], linewidth = 1.5, linestyle = :dash, color = :red2, label = false)
    plot!(xlabel = "Time (µs)", ylabel = "Voltage (V)", legend = :topright, dpi = 300)
    xl = xlims() 

    if showResiduals 
        PltStep = 100
        pres = plot(time[timeIdx:PltStep:end], 1e3.*(tail_dsp_cold[timeIdx:PltStep:end] .- (wvfs_cold[i].signal[timeIdx:PltStep:end].- baseline_cold[i])),
                    marker = :circle,
                    markersize = 0, markerstrokewidth = 0,
                    linewidth = 2.5,
                    xlabel = "Time", 
                    ylabel = "Residuals (mV)", 
                    legend = false,
                    color = :dodgerblue,
                    xlims = xl)
        plot!(time[timeIdx:PltStep:end], 1e3.*(tail_dsp_warm[timeIdx:PltStep:end] .- (wvfs_warm[i].signal[timeIdx:PltStep:end].- baseline_warm[i])),
            marker = :circle,
            linewidth = 2.5, markersize = 0, markerstrokewidth = 0,
            alpha = 0.5,
            legend = false,
            yticks = collect(-4:2:4),
            ylims = (-5, 5),
            color = :red2,
            xlims = xl)            
        hline!([0], color = :black, linewidth = 2, linestyle = :dash, label = false)

        layout = @layout [a{0.65h}; b]
        ptot = plot(plts, pres, layout = layout, size = (600, 600))
        resStr = "Residuals_"
    else
        resStr = ""
    end 

    if saveplot
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)Decaytime/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "DecayTime_Temperature_$(resStr)$(set_cold.Cinj_fF)fF_Rf$(set_cold.Rf_MOhm)M_pulser$pulsersavestr.png"
        savefig(fname)
        @info "Saving figure to $fname"
    end 

end 
display(ptot)