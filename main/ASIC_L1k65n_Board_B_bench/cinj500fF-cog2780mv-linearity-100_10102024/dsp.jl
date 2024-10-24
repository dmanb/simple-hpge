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
recompute = false 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

# function definitions 
function do_benchtest_dsp(pd_data, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_mV)}; recompute = false, pulserFlag = true) 
    dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)
    try 
        # get data folder 
        dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
        dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_mV, set.Temp_K; PulserChargeInj_Mode = "mV")
        pulser_file = replace(dsp_file, ".json" => "_pulser.json")

        if isfile(dsp_file) && !recompute
                dsp_par = readlprops(dsp_file)
                dsp_pulser = readlprops(pulser_file)
                println("Reading DSP results from file: $dsp_file")
            else
                # read data 
                filepath = ENV["ASIC_DATA"] * pd_data.datafolder 
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
            return dsp_par, dsp_pulser   
    catch e
        println("Error: $e")
        return NaN
    end
end


Benchtest_settings = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
            Rf_MOhm = 500,  # feedback resistor in Mega Ohm
            Temp_K = 300,  # temperature in Kelvin 
            PulserChargeInj_mV = 100)

# do dsp for 1 benchmark settings
dsp_par, pulser_par = do_benchtest_dsp(pd_data, Benchtest_settings; recompute= recompute, pulserFlag = true );
typeof(dsp_par)




# loop over different pulser charges
PulserChargesInj_mV =  [5, 10, 25, 50, 100, collect(200:200:2000)...] 
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
