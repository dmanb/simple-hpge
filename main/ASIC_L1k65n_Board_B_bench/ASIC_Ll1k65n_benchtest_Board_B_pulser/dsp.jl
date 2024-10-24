relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
import Pkg
Pkg.activate("$(@__DIR__)/$relPath/")
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
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

# settings 
plotFlag = true
recompute = false 
pd_data = readlprops("$(@__DIR__)/data_config.json")

# function definitions 
function do_benchtest_dsp(pd_data, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_mV)}; recompute = false, pulserFlag = true) 
    # get dsp config 
    dsp_config_set = readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config_overwrite.json")
    dsp_config_tmp = readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json")
    dsp_config = DSPConfig(merge(dsp_config_tmp.default, get(dsp_config_set, Symbol("$(set.Cinj_fF)fF"), PropDict())))
    try  
        dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
        dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_mV, set.Temp_K; PulserChargeInj_Mode = "mV")
        pulser_file = replace(dsp_file, ".json" => "_pulser.json")

        if isfile(dsp_file) && !recompute
                dsp_pd = readlprops(dsp_file)
                pulser_pd = readlprops(pulser_file)
                println("Reading DSP results from file: $dsp_file")
            else
                # read data 
                filepath = ENV["ASIC_DATA"] * pd_data.datafolder * "Linearity_Cinj$(set.Cinj_fF)fF/"
                filename = @sprintf("%dk-%dmV-%dfF_000_ALL.csv",set.Temp_K, set.PulserChargeInj_mV, set.Cinj_fF)
                filepath = filepath * filename
                wvfs, pulser, _ = read_file_csv_oscilloscope(filepath;  heading = 17, nChannels = 2) 

                # dsp: this does: baseline shift, pole-zero correction, t0 determination, energy filters,....
                dsp_par = simple_dsp(wvfs, dsp_config)
                # create a PropDict save results as json file
                dsp_pd = PropDict(Dict(varname => getproperty(dsp_par, varname) for varname in columnnames(dsp_par)))
                if !ispath(dsp_folder)
                    mkpath(dsp_folder)
                end
                writelprops(dsp_file, dsp_pd)
                println("Writing DSP results to file: $dsp_file")   

                if pulserFlag == true 
                    dsp_pulser = simple_dsp(pulser, dsp_config)
                    # create a PropDict save results as json file
                    pulser_pd = PropDict(Dict(varname => getproperty(dsp_pulser, varname) for varname in columnnames(dsp_pulser)))
                    if !ispath(dsp_folder)
                        mkpath(dsp_folder)
                    end
                   
                    writelprops(pulser_file, pulser_pd)
                    println("Writing pulser DSP results to file: $pulser_file") 
                end
            end 
            return dsp_pd, pulser_pd   
    catch e
        println("Error: $e")
        return NaN
    end
end


Benchtest_settings = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
            Rf_MOhm = 500,  # feedback resistor in Mega Ohm
            Temp_K = 77,  # temperature in Kelvin 
            PulserChargeInj_mV = 10)

# do dsp for 1 benchmark settings
recompute = true
dsp_par, pulser_par = do_benchtest_dsp(pd_data, Benchtest_settings; recompute= recompute, pulserFlag = true );
typeof(dsp_par)

# loop over different pulser charges
if Benchtest_settings.Cinj_fF == 500
    if Benchtest_settings.Temp_K == 300
        PulserChargesInj_mV =  [10, 30, 50, 70, 100, 130, 150, 170, 200, 230, 250, 280, 300, 320] 
    elseif Benchtest_settings.Temp_K == 77
        PulserChargesInj_mV =  [5, 10, 25, 50, 100, 200, 400] 
    end
elseif Benchtest_settings.Cinj_fF == 3000
    PulserChargesInj_mV =  [2, 5, 7, collect(10:5:50)...] 
end
Status = fill(false, length(PulserChargesInj_mV))
for i in eachindex(PulserChargesInj_mV)
    @info "Processing PulserChargeInj_mV = $(PulserChargesInj_mV[i]) ($i/$(length(PulserChargesInj_mV)))"
    Benchtest_settings = merge(Benchtest_settings, (PulserChargeInj_mV = PulserChargesInj_mV[i],) )
    dsp_par = do_benchtest_dsp(pd_data, Benchtest_settings; recompute = recompute, pulserFlag = true)[1];
    if isa(dsp_par, PropDict)
        Status[i] = true
    end
end
print(Status)