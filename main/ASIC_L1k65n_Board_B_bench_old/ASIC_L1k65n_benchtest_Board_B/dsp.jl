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

Benchtest_settings = (Cinj_fF = 500, # capacitance of the ASIC in femto Farrad 
                    Rf_MOhm = 500,  # feedback resistor in Mega Ohm
                    Temp_K = 77,  # temperature in Kelvin 
                    PulserChargeInj_keV = 1500) # pulser charge injected in keV

# function definitions 
function do_benchtest_dsp(pd_data, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_keV)}; recompute = false, pulserFlag = true) 
    dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)
    try 
        # get data folder 
        dsp_folder = "$(@__DIR__)$(pd_data.dspfolder)"
        dsp_file = dsp_folder * benchtest_filename_dsp(set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K)
        pulser_file = replace(dsp_file, ".json" => "_pulser.json")

        if isfile(dsp_file) && !recompute
                dsp_par = readlprops(dsp_file)
                dsp_pulser = readlprops(pulser_file)
                println("Reading DSP results from file: $dsp_file")
            else
                # read data 
                folder = benchtest_filename_raw(pd_data, set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K) 
                wvfs, pulser, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 2)

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
                else 
                    dsp_pulser = PropDict()
                end
            end 
            return (dsp_par , dsp_pulser) 
    catch e
        println("Error: $e")
        return NaN
    end
end

function plot_benchtest_dsp(dsp_par::Union{PropDict, Table}, pulser_par::Union{PropDict, Table}, pd_data::PropDict, set::NamedTuple{(:Cinj_fF, :Rf_MOhm, :Temp_K, :PulserChargeInj_keV)}; saveplot = true)
    # some plotting defaults 
    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 10,
    legendforegroundcolor = :silver)

    nbins = 0.25*length(dsp_par.t90) < 20 ? 20 :  round(Int, 0.25*length(dsp_par.t90)) 

    # decay time
    p1 = stephist(dsp_par.tail_τ, bins= nbins, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))
    annotate!(xlims()[2] - 0.95*(xlims()[2] - xlims()[1]), 0.9*ylims()[2], 
                text(@sprintf("mean = %.1f µs", ustrip(mean(dsp_par.tail_τ))), :black, :left))
                # text(@sprintf("mean = %.1f µs \n std = %.1f µs", ustrip(mean(dsp_par.tail_τ)), ustrip(std(dsp_par.tail_τ))), :darkgret, :left))
    # rise time 
   
    timerange = 1e3 * ustrip(maximum(dsp_par.t90 .- dsp_par.t10) - minimum(dsp_par.t90 .- dsp_par.t10));
    if timerange < 1
        xt = round.(range(1e3 * 0.95 * minimum(ustrip.(dsp_par.t90 .- dsp_par.t10)), 1e3 * 1.05 * maximum(ustrip.(dsp_par.t90 .- dsp_par.t10)), 4), digits = 2)
    elseif timerange < 3
        xt = (floor(Int, 2* 1e3 .* minimum(ustrip.(dsp_par.t90 .- dsp_par.t10)))/2):0.5:ceil(Int, 2 * 1e3 .* maximum(ustrip.(dsp_par.t90 .- dsp_par.t10)))/2
    else
        if timerange < 5 
            tstep = 1
        elseif timerange < 15
            tstep = 2
        else
            tstep = 5
        end
        xt = (floor(Int, 1e3 .* minimum(ustrip.(dsp_par.t90 .- dsp_par.t10)))):tstep:ceil(Int, 1e3 .* maximum(ustrip.(dsp_par.t90 .- dsp_par.t10)))
    end

    p2 = stephist(1e3 .* ustrip.(dsp_par.t90 .- dsp_par.t10), 
        bins= nbins, xlabel= "Rise time (ns)", ylabel="Occurrence", 
        xticks = xt, xlims = (minimum(xt)-diff(xt)[1]/2, maximum(xt)+diff(xt)[1]/2),
        label = @sprintf(""),
        fill = true, color = :red2, ylims = (0, :auto))
    annotate!(xlims()[2] - 0.95*(xlims()[2] - xlims()[1]), 0.9*ylims()[2], 
        text(@sprintf("mean = %.1f ns", 1e3 * ustrip(mean(dsp_par.t90 .- dsp_par.t10))), :black, :left))

    # baseline std in keV
    #1e3.*filter(x-> x<= quantile(dsp_par.blsigma, 0.99), dsp_par.blsigma
    p3 = stephist(1e-3 .* csa_voltage2charge.(set.Cinj_fF, dsp_par.blsigma), bins=20, xlabel= "Baseline σ (keV)", ylabel="Occurrence", 
        xformatter = x -> @sprintf("%.2f", x),
        label = false, fill = true, color = :orange, ylims = (0, :auto))
    annotate!(xlims()[2] - 0.95*(xlims()[2] - xlims()[1]), 0.9*ylims()[2], 
        text(@sprintf("mean = %.1f keV", ustrip(mean(1e-3 .* csa_voltage2charge.(set.Cinj_fF, dsp_par.blsigma)))), :black, :left))

    # energy with trapezoidal filter  in Volts 
    charge_MeV = 1e-6 .* csa_voltage2charge.(set.Cinj_fF, dsp_par.e_trap)
    pulser_MeV = 1e-6 .* csa_voltage2charge.(set.Cinj_fF, pulser_par.e_trap)
    if mean(charge_MeV)<0.9
        charge_plt = charge_MeV * 1e3
        pulser_plt = pulser_MeV * 1e3
        ustr = "keV"
        PulserChargeInj_lbl = set.PulserChargeInj_keV
        xt = round.(range(minimum(charge_plt), maximum(charge_plt), length = 3), digits = 0)
        xtp = round.(range(minimum(pulser_plt), maximum(pulser_plt), length = 3), digits = 1)
    else
        charge_plt = charge_MeV
        pulser_plt = pulser_MeV
        ustr = "MeV"
        PulserChargeInj_lbl = set.PulserChargeInj_keV/1e3
        xt = round.(range(minimum(charge_plt), maximum(charge_plt), length = 3), digits = 3)
        xtp = round.(range(minimum(pulser_plt), maximum(pulser_plt), length = 3), digits = 3)
    end
    p4 = stephist(charge_plt, bins = nbins, 
        xlabel= "Amplitude ($ustr)", ylabel="Occurrence", 
        label = false, #, @sprintf("Trap. filter: σ = %.1f keV", 1e3*std(charge_MeV)), 
        fill = true, 
        color = :dodgerblue, 
        xticks = xt, xlims = (minimum(xt)-diff(xt)[1]/2, maximum(xt)+diff(xt)[1]/2),
        ylims = (0, :auto))
    annotate!(xlims()[2] - 0.95*(xlims()[2] - xlims()[1]), 0.9*ylims()[2], 
        text(@sprintf("mean = %.2f %s", ustrip(mean(charge_plt)), ustr), :black, :left))
    
    p5 = stephist(pulser_plt, bins = nbins, 
        xlabel= "Pulser amplitude ($ustr)", ylabel="Occurrence", 
        label = false, #, @sprintf("Trap. filter: σ = %.1f keV", 1e3*std(charge_MeV)), 
        fill = true, 
        xticks = xtp, xlims = (minimum(xtp)-diff(xtp)[1]/2, maximum(xtp)+diff(xtp)[1]/2),
        color = :green3,
        ylims = (0, :auto))
    annotate!(xlims()[2] - 0.95*(xlims()[2] - xlims()[1]), 0.9*ylims()[2], 
        text(@sprintf("mean = %.2f %s", ustrip(mean(pulser_plt)), ustr), :black, :left))

    # labels 
    (csaname, board) = CSAname(pd_data)
    
    # overview 
    ptot = plot(p1, p2, p3, p4, p5,
        layout = (3,2), size = (1200, 1400), 
        top_margin = 5mm, left_margins = 1mm , right_margin = 5mm,
        thickness_scaling = 1.3,
        plot_title = "Bench test: $(csaname)-Board $(board) \n" * "C = $(set.Cinj_fF) fF, R = $(set.Rf_MOhm) MOhm, T = $(set.Temp_K) K, " * L"$Q_\textrm{inj}$" * " = $(PulserChargeInj_lbl) $ustr",
        plot_titlefontsize = 13,
        plot_titlefontcolor = :darkgrey,
        dpi = 300)

    if saveplot
        # # save figure 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)DSP/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "DSPpars_" * benchtest_filename_plot(set.Cinj_fF, set.Rf_MOhm, set.PulserChargeInj_keV, set.Temp_K)
        savefig(ptot, fname)
        @info "Saving figure to $fname"
    end 
    display(ptot)
    return ptot 
end

# do dsp for 1 benchmark settings
dsp_par, pulser_par = do_benchtest_dsp(pd_data, Benchtest_settings; recompute = recompute, pulserFlag = true);
ptot = plot_benchtest_dsp(dsp_par, pulser_par, pd_data, Benchtest_settings; saveplot = true)

# loop over different pulser charges
if Benchtest_settings.Temp_K == 77
    # PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 6, 10000,15, 30, 100, 300]) # pulser charge injected in keV
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 6 ,15, 30, 100, 300]) # pulser charge injected in keV
elseif Benchtest_settings.Temp_K == 300
    # PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 5000, 6, 10000,15, 30]) # pulser charge injected in keV
    PulserChargesInj_keV =  sort([1000, 1500, 2500, 3000, 6 ,15, 30]) 
end
Status = fill(false, length(PulserChargesInj_keV))
for i in eachindex(PulserChargesInj_keV)
    @info "Processing PulserChargeInj_keV = $(PulserChargesInj_keV[i]) ($i/$(length(PulserChargesInj_keV)))"
    Benchtest_settings = merge(Benchtest_settings, (PulserChargeInj_keV = PulserChargesInj_keV[i],) )
    local dsp_par, dsp_pulser = do_benchtest_dsp(pd_data, Benchtest_settings; recompute = recompute, pulserFlag = true);
    if isa(dsp_par, PropDict) | isa(dsp_par, Table)
        Status[i] = true
        if plotFlag == true 
            plot_benchtest_dsp(dsp_par, dsp_pulser, pd_data, Benchtest_settings; saveplot = true)
        end
    end
end
print(Status)




# if plotFlag
#     for (i, PulserChargeInj_keV) in enumerate(PulserChargesInj_keV)  
#         # get data folder 
#         dsp_file = dsp_folder * benchtest_filename_dsp(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
#         if !isfile(dsp_file) 
#             continue
#         end
#         dsp_par = readlprops(dsp_file)
#         println("Reading DSP results from file: $dsp_file")

#         # #########  look at results #########
#         # some plotting defaults 
#         default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
#         xguidefontsize = 16, xtickfontsize = 12,
#         yguidefontsize = 16, ytickfontsize = 12,
#         legendfontsize = 14, titlefontsize = 10,
#         legendforegroundcolor = :silver)

#         vars = keys(dsp_par) # overview of variables in dsp_par
#         nbins = 0.25*length(dsp_par.t90) < 20 ? 20 :  round(Int, 0.25*length(dsp_par.t90)) 

#         # decay time
#         p1 = stephist(dsp_par.tail_τ, bins= nbins, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))

#         # rise time 
#         p2 = stephist(1e3 .* ustrip.(dsp_par.t90 .- dsp_par.t10), 
#             bins= nbins, xlabel= "Rise time (ns)", ylabel="Occurrence", 
#             xformatter = x -> @sprintf("%.0f", x),
#             label = L"$t_{90} - t_{10}$", 
#             fill = true, color = :red2, ylims = (0, :auto))

#         # baseline std
#         p3 = stephist(1e3.*filter(x-> x<= quantile(dsp_par.blsigma, 0.95), dsp_par.blsigma), bins=20, xlabel= "Baseline σ (keV)", ylabel="Occurrence", 
#         label = false, fill = true, color = :orange, ylims = (0, :auto))

#         # energy with trapezoidal filter  in Volts 
#         charge_MeV = 1e-6 .* csa_voltage2charge.(Cinj_fF, dsp_par.e_trap)
#         p4 = stephist(charge_MeV, bins = nbins, 
#             xlabel= "Amplitude (MeV)", ylabel="Occurrence", 
#             label = @sprintf("Trap. filter: σ = %.1f keV", 1e3*std(charge_MeV)), 
#             fill = true, 
#             color = :dodgerblue, 
#             ylims = (0, :auto))

#         # overview 
#         plot(p1, p2, p3, p4, layout = (2,2), size = (1200, 800), 
#             left_margins = 2mm , thickness_scaling = 1.4,
#             dpi = 300)

#         # # save figure 
#         fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)/DSP/"
#         if !ispath(fpath)
#             mkpath(fpath)
#         end
#         fname = fpath * "DSPpars_" * benchtest_filename_plot(Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K)
#         savefig(fname)
#     end
# end 