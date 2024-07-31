############################################################
# PACKAGES #
############################################################

Pkg.activate("$(@__DIR__)/../")
import Pkg
using TypedTables, Dates
using RadiationDetectorSignals
using Unitful, Formatting, LaTeXStrings, Measures, Measurements
using Measurements: value as mvalue, uncertainty as muncertainty
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP: get_fltpars
using ArraysOfArrays
using Distributed
using RadiationDetectorDSP
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendHDF5IO, LegendDSP, LegendSpecFits, LegendDataTypes
using IntervalSets, TypedTables, StatsBase, PropDicts
using Unitful, Formatting, Printf, Measures, Dates
using Distributed, ProgressMeter, TimerOutputs
using RadiationDetectorDSP
using ArraysOfArrays
using HDF5
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP: get_fltpars
using TypedTables, Dates
using Measurements: value as mvalue, uncertainty as muncertainty
using LegendDataManagement: readlprops, writelprops
using CSV

############################################################
############################################################


## path = path from your current directory to the folder with the data
## heading = how many rows to skip before reading data, optional input, default to 3 since that is our format
function read_folder(relative_path::String; heading = 1)

    ## creates path to data folder #
    data_folder = joinpath(@__DIR__, relative_path)
    # print(isdir(data_folder))

    # arrays needed for data collection ## 
    times = Vector{typeof(0.0u"µs":1.0u"µs":1.0u"µs")}()
    voltages = Vector{Vector{Float64}}()
    wvfarray = ArrayOfRDWaveforms

    ## for each file in the directory with ASIC in the title, it will read the .csv file
    for file in readdir(data_folder)
        if occursin("ASIC", file)
             
            path = joinpath(data_folder, file)
            df = CSV.read(path, DataFrame; delim=',', header = 2)
            # Table(file)
            # println("Column names: ", names(df))

            # creating time, voltage from .csv file
            asic_voltages = df[:, "ASIC Voltage (V)"]
            csv_time = uconvert.(u"µs", (df[:, "Time (s)"] .- df[1, "Time (s)"])*u"s")

            ## formatting time properly for the RDWaveform object, then appending 
            timestep = csv_time[2] - csv_time[1]
            time = 0u"µs":timestep:(length(asic_voltages) - 1)*timestep
            push!(times, time)
        
            # putting vector of voltages into voltages object
            push!(voltages, asic_voltages)
        end
    end
    wvfarray = ArrayOfRDWaveforms((times, voltages))

    # gives us an array of RDWaveforms, one for each .csv file. 
    return wvfarray
end


## modify the path as required for location of config file
function get_decay_times(wvfs)
    path_config = "$(@__DIR__)/../config/dsp_config.json"
    # get DSP configuration data --> Can be modified in .json filepath = " $(@__DIR__)/../config/dsp_config.json"
    dsp_meta = readlprops(path_config)
    dsp_config = DSPConfig(dsp_meta.default)
    dsp_config.bl_window
    dsp_config.tail_window  

    ## extract decay times of all waveforms with a simple DSP
    decay_times = Vector{Float64}()
    decay_times = dsp_decay_times(wvfs, dsp_config)
    return decay_times
end




function simple_dsp(wvfs, decays)
    ### config stuff
    path_config = "$(@__DIR__)/../config/dsp_config.json"
    # get DSP configuration data --> Can be modified in .json filepath = " $(@__DIR__)/../config/dsp_config.json"
    dsp_meta = readlprops(path_config)
    config = DSPConfig(dsp_meta.default)
    config.bl_window
    config.tail_window


    bl_window                = config.bl_window
    t0_threshold             = config.t0_threshold
    tail_window              = config.tail_window
    inTraceCut_std_threshold = config.inTraceCut_std_threshold
    sg_flt_degree            = config.sg_flt_degree
    current_window           = config.current_window
    qdrift_int_length        = config.qdrift_int_length
    lq_int_length            = config.lq_int_length

    #### all filter parameters set to default since we do not have access to full library 
    trap_rt, trap_ft = get_fltpars(PropDict(),:trap, config)
    cusp_rt, cusp_ft = get_fltpars(PropDict(), :cusp, config)
    zac_rt, zac_ft   = get_fltpars(PropDict(), :zac, config)
    sg_wl            = get_fltpars(PropDict(), :sg, config)

    # get CUSP and ZAC filter length and flt scale
    flt_length_zac              = config.flt_length_zac
    zac_scale                   = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))
    flt_length_cusp             = config.flt_length_cusp
    cusp_scale                  = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # get number of samples the waveform is saturated at low and high of FADC range
    # bit_depth = config.kwargs_pars.fc_bit_depth # of FlashCam FADC
    # sat_low, sat_high = 0, 2^bit_depth - bit_depth
    # sat_stats = saturation.(wvfs, sat_low, sat_high)

    # set tau for CUSP filter to very high number to switch of CR filter
    τ_cusp = 10000000.0u"µs"
    τ_zac = 10000000.0u"µs"


    ################## ACTUAL WAVEFORM FILTERING AND RECONSTRUCTION, ANALYSIS BEGINS HERE ##################
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
    
    # get wvf maximum
    wvf_max = maximum.(wvfs.signal)
    wvf_min = minimum.(wvfs.signal)
    
    # extract decay times
    tail_stats = tailstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))
    
    # deconvolute waveform 
    # --> wvfs = wvfs_pz
    deconv_flt = InvCRFilter(decays[1])
    wvfs = deconv_flt.(wvfs)
    
    # get tail mean, std and slope
    pz_stats = signalstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))
    
    # t0 determination
    t0 = get_t0(wvfs, t0_threshold; flt_pars=config.kwargs_pars.t0_flt_pars, mintot=config.kwargs_pars.t0_mintot)
    
    # if all waveforms are saturated set threshold to 1.0 to avoid numerical problems
    # replace!(wvf_max, zero(wvf_max[1]) => one(wvf_max[1]))
    
    # get threshold points in rise
    t10 = get_threshold(wvfs, wvf_max .* 0.1; mintot=config.kwargs_pars.tx_mintot)
    t50 = get_threshold(wvfs, wvf_max .* 0.5; mintot=config.kwargs_pars.tx_mintot)
    t80 = get_threshold(wvfs, wvf_max .* 0.8; mintot=config.kwargs_pars.tx_mintot)
    t90 = get_threshold(wvfs, wvf_max .* 0.9; mintot=config.kwargs_pars.tx_mintot)
    t99 = get_threshold(wvfs, wvf_max .* 0.99; mintot=config.kwargs_pars.tx_mintot)
        
    drift_time = uconvert.(u"ns", t90 - t0)
    
    # get Q-drift parameter
    qdrift = get_qdrift(wvfs, t0, qdrift_int_length; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length)
    
    # get LQ parameter
    lq  = get_qdrift(wvfs, t80, lq_int_length; pol_power=config.kwargs_pars.int_interpolation_order, sign_est_length=config.kwargs_pars.int_interpolation_length)
    
    # robust energy reconstruction with long, middle and short rise and flat-top times
    uflt_10410 = TrapezoidalChargeFilter(10u"µs", 4u"µs")
    e_10410  = maximum.((uflt_10410.(wvfs)).signal)

    uflt_535 = TrapezoidalChargeFilter(5u"µs", 3u"µs")
    e_535  = maximum.((uflt_535.(wvfs)).signal)
    
    uflt_313 = TrapezoidalChargeFilter(3u"µs", 1u"µs")
    e_313  = maximum.((uflt_313.(wvfs)).signal)
    
    # signal estimator for precise energy reconstruction
    signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))
    
    # get trap energy of optimized rise and flat-top time
    uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
    
    e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))
    
    # get cusp energy of optimized rise and flat-top time
    uflt_cusp_rtft = CUSPChargeFilter(cusp_rt, cusp_ft, τ_cusp, flt_length_cusp, cusp_scale)
    
    e_cusp = signal_estimator.(uflt_cusp_rtft.(wvfs), t50 .+ (flt_length_cusp /2))
    
    # get zac energy of optimized rise and flat-top time
    uflt_zac_rtft = ZACChargeFilter(zac_rt, zac_ft, τ_zac, flt_length_zac, zac_scale)
    
    e_zac = signal_estimator.(uflt_zac_rtft.(wvfs), t50 .+ (flt_length_zac /2))
    
    # extract current with optimal SG filter length with second order polynominal and first derivative
    wvfs_sgflt_deriv = SavitzkyGolayFilter(sg_wl, sg_flt_degree, 1).(wvfs)
    a_sg = get_wvf_maximum.(wvfs_sgflt_deriv, leftendpoint(current_window), rightendpoint(current_window))
    
    a_60 = get_wvf_maximum.(SavitzkyGolayFilter(60u"ns", sg_flt_degree, 1).(wvfs), leftendpoint(current_window), rightendpoint(current_window))
    a_100 = get_wvf_maximum.(SavitzkyGolayFilter(100u"ns", sg_flt_degree, 1).(wvfs), leftendpoint(current_window), rightendpoint(current_window))
    a_raw = get_wvf_maximum.(DifferentiatorFilter(1).(wvfs), leftendpoint(current_window), rightendpoint(current_window))
    
    # get in-trace pile-up
    inTrace_pileUp = get_intracePileUp(wvfs_sgflt_deriv, inTraceCut_std_threshold, bl_window; mintot=config.kwargs_pars.intrace_mintot)
        
    # get position of current rise
    thres = maximum.(wvfs_sgflt_deriv.signal) .* 0.5
    # replace!(thres, zero(thres[1]) => one(thres[1]))
    
    t50_current = get_threshold(wvfs_sgflt_deriv, thres; mintot=config.kwargs_pars.tx_mintot)
    
    # invert waveform for DC tagging
    # wvfs --> wvfs_pz_inv
    wvfs = multiply_waveform.(wvfs, -1.0)
    
    # get inverted waveform maximum for long and short filter
    e_10410_max_inv  = maximum.(uflt_10410.(wvfs).signal)
    
    e_313_max_inv  = maximum.(uflt_313.(wvfs).signal)
    
    # t0 determination
    t0_inv = get_t0(wvfs, t0_threshold; mintot=config.kwargs_pars.t0_mintot)
    
    # output Table 
    return TypedTables.Table(blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        tailmean = pz_stats.mean, tailsigma = pz_stats.sigma, tailslope = pz_stats.slope, tailoffset = pz_stats.offset,
        t0 = t0, t10 = t10, t50 = t50, t80 = t80, t90 = t90, t99 = t99,
        t50_current = t50_current, 
        drift_time = drift_time,
        tail_τ = tail_stats.τ, tail_mean = tail_stats.mean, tail_sigma = tail_stats.sigma,
        e_max = wvf_max, e_min = wvf_min,
        e_10410 = e_10410, e_535 = e_535, e_313 = e_313,
        e_10410_inv = e_10410_max_inv, e_313_inv = e_313_max_inv,
        t0_inv = t0_inv,
        e_trap = e_trap, e_cusp = e_cusp, e_zac = e_zac, 
        qdrift = qdrift, lq = lq,
        a_sg = a_sg, a_60 = a_60, a_100 = a_100, a_raw = a_raw,
        )
end 




