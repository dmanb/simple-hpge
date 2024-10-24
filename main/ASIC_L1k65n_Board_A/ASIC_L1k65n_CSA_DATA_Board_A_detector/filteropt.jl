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
include("$(@__DIR__)/$relPath/utils/utils_naming.jl")
include("$relPath/utils/utils_aux.jl")
include("$relPath/utils/utils_plot.jl")

include("$relPath/src/filteropt_trap.jl")
Cinj_fF = 3000
Cf_fF = 500
noiseAmpFac = 1
nwvf_max = 1000

# get data configuration (where to and save stuff)
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

# load quality cuts (qc) and select only good waveforms
qc_folder = "$(@__DIR__)$(pd_data.dspfolder)"
qc_file = qc_folder * "QualityCuts.json"
cuts_qc = readlprops(qc_file).wvf_keep.all
wvfs_sel = findall(cuts_qc)[1:nwvf_max]

# read data 
folder = pd_data.datafolder
if !ispath(get(ENV, "ASIC_DATA",".") * folder)
    println("Data folder does not exist: $(get(ENV, "ASIC_DATA",".") * folder)")
end
wvfs_raw, _ = read_folder_csv_oscilloscope(folder; heading = 14, nChannels = 1, nwvfmax = wvfs_sel) # make sure you can read all waveforms 

# get dsp configuration
dsp_config = DSPConfig(readlprops("$(@__DIR__)" * pd_data.configfolder * "dsp_config.json").default)

# shift waveforms
bl_window                   = dsp_config.bl_window
bl_stats = signalstats.(wvfs_raw, leftendpoint(bl_window), rightendpoint(bl_window))
wvfs = shift_waveform.(wvfs_raw, -bl_stats.mean)

# calculate ENC noise  
wvfs_timestep = diff(wvfs_raw.time[1])[1]
sample_timestep =  wvfs_timestep
sample_timestep = round(Int, sample_timestep/wvfs_timestep) * wvfs_timestep # make sure it is a multiple of wvfs_timestep
sample_binstep = round(Int,sample_timestep/ wvfs_timestep)
bl_start = findfirst(wvfs[1].time .>= leftendpoint(bl_window))
bl_stop = findfirst(wvfs[1].time .>= rightendpoint(bl_window))-1

# construct trapezoidal filter and calculate ENC noise and convert to charge 
trap_rts = 10 .^ (range(log10(1), stop = log10(20), length = 20)) .*u"µs" 
trap_ft = 0.6u"µs"

ENC = [std(filter(isfinite, vcat(map(x-> x.signal[bl_start:sample_binstep:bl_stop], TrapezoidalChargeFilter(rt, trap_ft).(wvfs))...))) for rt in trap_rts]
csa_gain = Cinj_fF / Cf_fF
ENC .= ENC ./noiseAmpFac ./ csa_gain
ENC_e = csa_voltage2charge.(Ref(Cinj_fF), ENC)  # convert to electronvolt equivalent

axFlag = :log 
xunit = "e"#"a.u."


begin 
    if xunit == "e"
        plotx = ENC_e
        xustr = L"$(e^-)$" 
    elseif xunit == "a.u."
        plotx = ENC
        xustr = "(a.u.)"
    end

    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
            xguidefontsize = 16, xtickfontsize = 12,
            yguidefontsize = 16, ytickfontsize = 12,
            legendfontsize = 12, titlefontsize = 14,
            legendtitlefontsize = 16,
            colorbar_titlefontsize = 16,
            legendforegroundcolor = :silver)

    yrange = diff([minimum(plotx), maximum(plotx)])[1]
    rt = ustrip.(collect(trap_rts))
    scatter(rt, plotx, 
        xlabel="Shaping time (µs)", 
        linewidth = 2, 
        markerstrokewidth = 0,
        color = :dodgerblue,
        linestyle = :solid,
        linecolor = :dodgerblue,
        markercolor = :dodgerblue,
        markersize = 2,
        ylabel= "\n ENC noise " * xustr,
        label = "Flat-top time = $(trap_ft)", 
        legend=:top, # 
        title =  "",#"Baseline noise L1k65n benchtest, CSA Board-B: \nCf = $(Cf_fF)fF,  $(replace(replace(replace(replace(dataset, "NOISE_" =>""),"_" => ", " ),"inj"=>"inj = "),"fF," => "fF, T =") )",
        titlefontsize = 11,
       # yformatter = log_tick_formatter_plain,
      #  yticks = [20:20:200..., 300:100:1000...],# ylims =  (floor(Int, (minimum(plotx) - 0.2 * yrange)/10)*10, ceil(Int, (maximum(plotx) + 0.05 * yrange)/10)*10),
       # xlims = (0.5, 100),#(minimum(rt) - 0.25, maximum(rt) + 0.25),
        dpi = 300,
        margins = 3mm)

    hline!([minimum(plotx)], color = :red2, linestyle = :dash, linewidth = 2, label = @sprintf("min = %.0f ", minimum(plotx)) * L"e^-" *
        "using " * L"$\tau_\textrm{rt} = $" * @sprintf("%.1f µs", ustrip(trap_rts[findfirst(plotx .== minimum(plotx))])) )
    if axFlag == :log
        plot!(yscale = :log10)
        plot!(xscale = :log10)
        axstr = "_loglog"
    else 
        axstr = ""
    end 
    # ppath = "$(@__DIR__)$(pd_data.figurefolder)Noise/"
    # if !ispath(ppath)
    #     mkpath(ppath)
    # end
    # pname = replace(replace(fname, fpath => ""), ".jld2" => "$(axstr).png")
    # savefig(ppath * pname)
    # @info "Saving plot to file:  figures/$pname"
    display(plot!())
end
