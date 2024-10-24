# trap filter:
# flat-top time = "gap time"
# shaping time/ rise time = averaging time
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
using Measures, Measurements
using Measurements: value as mvalue
using Statistics
# using Plots
using CairoMakie
using Printf, LaTeXStrings
using JLD2

include("$(@__DIR__)/$relPath/utils/utils_IO.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")

# settings 
recompute = false
# set = (Cinj_fF = 3000, Temp_K = 300, noiseAmpFac = 48, Cf_fF = 500, extraStr = "_LDO")
# set = (Cinj_fF = 3000, Temp_K = 300, noiseAmpFac = 48, Cf_fF = 500, extraStr = "_400ns")
# set = (Cinj_fF = 3000, Temp_K = 300, noiseAmpFac = 48, Cf_fF = 500, extraStr = "_LDO_SUPERCAP")
# set = (Cinj_fF = 3000, Temp_K = 300, noiseAmpFac = 1, Cf_fF = 500, extraStr = "_LDO_NO_AMP")
# set = (Cinj_fF = 3000, Temp_K = 300, noiseAmpFac = 100, Cf_fF = 500, extraStr = "")
# set = (Cinj_fF = 3000, Temp_K = 77, noiseAmpFac = 100, Cf_fF = 500, extraStr = "")
set = (Cinj_fF = 500, Temp_K = 300, noiseAmpFac = 100, Cf_fF = 500, extraStr = "")


dataset = noise_dataset_str(set.Cinj_fF, set.Temp_K, set.noiseAmpFac) * set.extraStr
if set.extraStr == "_400ns"
    trap_rts = 10 .^ (range(log10(0.2), stop = log10(100), length = 100)) .*u"µs" 
else
    trap_rts = 10 .^ (range(log10(0.1), stop = log10(100), length = 100)) .*u"µs" 
end
trap_ft = 0.6u"µs"

# load or compute ENC 
pd_data = readlprops("$(@__DIR__)/data_config.json")
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

fpath = "$(@__DIR__)$(pd_data.noisefolder)"
fname = fpath * replace("ENC_scan_$(replace(dataset, "NOISE_" =>""))_ft$(round(ustrip(trap_ft), digits = 1))_rt-$(length(trap_rts))steps-min$(round(minimum(ustrip.(trap_rts)), digits = 1))-max$(round(ustrip(maximum(trap_rts)), digits =1))", "." => "p") * ".jld2"
if isfile(fname) & !recompute 
    f = jldopen(fname)
    ENC = f["ENC"]
    ENC_e = f["ENC_e"]
    @info "Reading ENC from file: \n$fname"
    close(f)
else
    # read data
    folder = pd_data.datafolder * dataset * "/"
    wvfs_raw, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1);
    @sprintf("Reading %d waveforms ", length(wvfs_raw))

    # shift waveforms
    bl_window                   = dsp_config.bl_window
    bl_stats = signalstats.(wvfs_raw, leftendpoint(bl_window), rightendpoint(bl_window))
    wvfs = shift_waveform.(wvfs_raw, -bl_stats.mean)
    
    # cut waveforms 1 
    wvfs_timestep = diff(wvfs_raw.time[1])[1]
    sample_timestep =  wvfs_timestep#2*maximum(trap_rts) + trap_ft + 5 * wvfs_timestep
    sample_timestep = round(Int, sample_timestep/wvfs_timestep) * wvfs_timestep # make sure it is a multiple of wvfs_timestep
    sample_binstep = round(Int,sample_timestep/ wvfs_timestep)
    # enc_pickoff  = range(sample_timestep, step = sample_timestep, stop = wvfs.time[1][end] - sample_timestep)
    # nSamples = length(wvfs) * length(enc_pickoff)

    # construct trapezoidal filter and calculate ENC noise and convert to charge 
    csa_gain = set.Cinj_fF / set.Cf_fF
    ENC = [std(filter(isfinite, vcat(map(x-> x.signal[1:sample_binstep:end], TrapezoidalChargeFilter(rt, trap_ft).(wvfs))...))) for rt in trap_rts]
    # signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))
    # ENC = [std(filter(isfinite, vcat([signal_estimator.(TrapezoidalChargeFilter(rt, trap_ft).(wvfs), time) for time in enc_pickoff]...)) ) for rt in trap_rts]
    ENC .= ENC./set.noiseAmpFac./csa_gain
    ENC_e = csa_voltage2charge.(Ref(set.Cinj_fF), ENC)  # convert to electronvolt equivalent

    if !ispath(fpath)
        mkpath(fpath)
    end

    jldsave(fname; ENC =  ENC, ENC_e = ENC_e, trap_rt = trap_rts, trap_ft = trap_ft)
    println("Writing ENC results to file: \n$fname") 
end 

# save 
ppath = "$(@__DIR__)$(pd_data.figurefolder)Noise/"
if !ispath(ppath)
    mkpath(ppath)
end

# find exponent of power law function
func_power = (x, exponent, scale)-> x^exponent * scale # power law function. in log-log scale, slope == exponent, offset == scale 
slope_exp = 0.5 # expected
timeanchor_p = trap_rts[findfirst(ENC_e .== minimum(ENC_e))] + 20u"µs"
idx = findfirst(trap_rts .>= timeanchor_p)
slope_data = (log10(ENC_e[end]) - log10(ENC_e[idx])) / (log10(ustrip(trap_rts[end]))  - log10(ustrip(trap_rts[idx]))) 
offset_data = ustrip.(ENC_e[idx:end] ./ trap_rts[idx:end].^slope_data)[1]
slope_p = 0.5 # parallel slope  --> exponent in power law function (lin-lin scale)
offset_p = median(ustrip.(ENC_e[idx:end] ./ trap_rts[idx:end].^slope_p))

# serial
timeanchor_s = trap_rts[findfirst(ENC_e .== minimum(ENC_e))] - 2u"µs"
idx_s = findfirst(trap_rts .>= timeanchor_s) 
offset_s = median(ustrip.(ENC_e[1:idx_s] ./ trap_rts[1:idx_s].^(-0.5)))

x = [collect(0.05:0.05:0.9)...,collect(1:1000)...]
yrange = maximum(ENC_e) - minimum(ENC_e)
xrange_dec = log10(ustrip(maximum(trap_rts)))  - log10(ustrip(minimum(trap_rts)))
# plot
Makie_theme(; fs = 20)


# setup 
begin 
    fig = Figure()
    ax = Axis(fig[1, 1],  
            xlabel= latexstring("\\textrm{Shaping time} \\,\\, \\tau_\\textrm{s} \\, \\textrm{(µs)}"), 
            ylabel= L"$\textrm{ENC noise} \,\, (e^-$)",
            xscale = log10,
            yscale = log10,
            title =  "Baseline noise L1k65n benchtest, CSA Board-B: \nCf = $(set.Cf_fF)fF,  $(replace(replace(replace(replace(dataset, "NOISE_" =>""),"_" => ", " ),"inj"=>"inj = "),"fF," => "fF, T =") )",
            xticks = 10.0 .^ collect(-1:1:3),
            yticks =  10.0 .^collect(-1:3), #sort([10, 50, 200:100:1000...,round(Int, minimum(ENC_e))]),#[ceil(Int, ymin/50)*50, round(Int, minimum(ENC_e)), ceil(Int, maximum(ENC_e)/100)*100],#[20:20:200..., 300:100:1000...],
            xtickformat = x -> [latexstring("10^{$(round(Int,log10(i)))}") for i in x],
            ytickformat = x -> [latexstring("10^{$(round(Int,log10(i)))}") for i in x], #log_tick_formatter,  
            yminorticks = [0.2:0.1:0.9..., 2:1:9..., 20:10:90..., 200:100:900..., 2000:1000:9000...],#IntervalsBetween(10), 
            xminorticks = IntervalsBetween(9),
            ygridvisible = true,
            xgridvisible = true,
            xminorgridvisible = true,
            yminorgridvisible = true,
            xminorgridstyle = :solid,
            yminorgridstyle = :solid,
            xgridwidth = 1.5,
            ygridwidth = 1.5,
            xminorgridwidth = 1,
            yminorgridwidth = 1,)

    # plot data
    ymin = minimum(ENC_e) > 70 ? 50 : 40 
    xlims!(ax, 0.1, 100)
    ylims!(ax, ymin , ceil(Int, maximum(ENC_e)/100)*100 + 50) 
    # lmin = hlines!(ax, [minimum(ENC_e)], xmin = 0, xmax = (log10(ustrip(trap_rts[findfirst(ENC_e .== minimum(ENC_e))]))  - log10(ustrip(minimum(trap_rts)))) / xrange_dec, color = :silver, linestyle = :solid, linewidth = 2)
    penc = scatter!(ax, ustrip.(trap_rts), ENC_e, markersize = 6, color = :dodgerblue)
  
    pname = replace(replace(fname, fpath => ""), ".jld2" => "_loglog.png")
    save(ppath * pname, fig)
    @info "Saving plot to file:  figures/$pname"

    line_p = lines!(x , func_power.(x, slope_p, offset_p) ,
                color = :Navy, linestyle = :solid, linewidth = 2)
    line_s = lines!(x , func_power.(x, -0.5, offset_s) ,
                color = :Navy, linestyle = :solid, linewidth = 2)

    # legend and other text 
    rt_min = ustrip(trap_rts[findfirst(ENC_e .== minimum(ENC_e))])
    rt_s = 10^quantile([log10(0.1), log10(rt_min)], 0.4)#(log10(0.1) - 0.3 * (log10(rt_min) + log10(0.1)) ) 
    rt_p = 10^quantile([log10(rt_min), log10(100)], 0.45)    #(log10(rt_min) + 0.3 * (log10(100) - log10(rt_min)) )     #mean([log10(100), log10(rt_min)])#^( 0.5 * (log10(100) - log10(rt_min)) )
    text!(rt_p, func_power(rt_p, slope_p, offset_p), 
        text =  latexstring("\$ \\textrm{ENC}^2_\\textrm{p} \\, \\propto \\tau_\\textrm{s}\$"), 
        rotation = angle =  2* atan(0.5) , color = :Navy, 
         offset = (0, -45), align = (:left, :bottom) )
    
    text!(rt_s, func_power(rt_s, -0.5, offset_s), 
        text =  latexstring("\$ \\textrm{ENC}^2_\\textrm{s} \\, \\propto 1/\\tau_\\textrm{s}\$"), 
        rotation = angle =  2*pi - 2* atan(0.5) , 
        color = :Navy,  offset = (0, -40),
        align = (:left, :bottom) )
    
    axislegend(ax, [penc, line_p], [latexstring("\\textrm{Data} \\, (\\tau_\\textrm{ft} = $(round(ustrip(trap_ft), digits = 1))\\,\\textrm{µs})"),
                                    L"\textrm{Expectation from noise model}"],
        position = :ct, markersize = 10)

    pname = replace(replace(fname, fpath => ""), ".jld2" => "_loglog_components.png")
    save(ppath * pname, fig)
    @info "Saving plot to file:  figures/$pname"
    # show plot
    fig
end


fig