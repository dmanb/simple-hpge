

function filteropt_trap(pd_data::PropDict, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig, Cinj_fF::Int; saveplot = false)
    decay_times = dsp_decay_times(wvfs, dsp_config)

    ## running the dsp_trap_rt_optimization function to get the ENC vs. shaping time plot
    grid_ft_trap = [0.2, 0.5, 1, 1.5, 2].*u"µs"
    ENC = NaN.*zeros(Float64, length(grid_ft_trap), length(dsp_config.e_grid_rt_trap))
    for (i, ft) in enumerate(grid_ft_trap)
        try
        local enc_trap_grid = dsp_trap_rt_optimization(wvfs, dsp_config, median(decay_times),; ft = ft)
        ENC[i, :] = [std(enc_trap_grid[j,:]) for j in eachindex(dsp_config.e_grid_rt_trap)]
        catch e
            ENC[i,:] .= NaN
        end
    end
    ENC[ENC .== 0] .= Inf # if ft > rt, ENC is not calculated and set to 0. remove for better visibility in plot
    ENC[isnan.(ENC)] .= Inf 
    ENC_eV = csa_voltage2charge.(Cinj_fF, ENC)

    # Plot
    # some plotting defaults 
    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
            xguidefontsize = 16, xtickfontsize = 12,
            yguidefontsize = 16, ytickfontsize = 12,
            legendfontsize = 14, titlefontsize = 14,
            legendtitlefontsize = 16,
            colorbar_titlefontsize = 16,
            legendforegroundcolor = :silver)

    # # plot result: ENC graph 
    grid_rt_trap = ustrip.(collect(dsp_config.e_grid_rt_trap))
    p1 = heatmap(grid_rt_trap, ustrip.(grid_ft_trap), ENC_eV,  
        ylabel = "Flat-top time (µs)", 
        xlabel = "Shaping time (µs)",
        zformatter = :plain,
        yformatter = x -> @sprintf("%.1f", x),
        yticks = (ustrip.(grid_ft_trap), [@sprintf("%.1f", x) for x in ustrip.(grid_ft_trap)]), 
        colorbar_title = "\n ENC noise (eV)",
        size = (650, 400),
        right_margin = 30mm,
        left_margin = 5mm,
        bottom_margin = 5mm)

    # find minimum ENC value: flat top time:
    _,  ENC_idx = findmin(ENC_eV)

    p2 = plot(grid_rt_trap,  ENC_eV[ENC_idx[1],:], 
        xlabel="Shaping time (µs)", 
        linewidth = 2.5, color = :red2,
        ylabel= "\n ENC noise (eV)",
        label = "Flattop time = $(grid_ft_trap[ENC_idx[1]])", 
        legendtitle = "$(length(wvfs)) waveforms",
        legend=:best, 
        yformatter = :plain)
    # vline!([ustrip(trap_rt)], label = "Def. shaping time = $trap_rt", color = :grey, linewidth = 2.5, linestyle = :dashdot)

    pENC = plot(p1, p2,  layout = (3,1), size = (900, 1500), 
        right_margin = 15mm, 
        left_margin = 5mm,
        dpi = 300,
        thickness_scaling = 1.3)

    if saveplot 
        fpath = "$(@__DIR__)$(pd_data.figurefolder)$(pd_data.datafolder)ENCnoise/"
        if !ispath(fpath)
            mkpath(fpath)
        end
        fname = fpath * "TrapFilterOpt"
        savefig(fname)
    end 
    display(pENC)
    return pENC
end