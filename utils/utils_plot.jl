function plot_dsp_par(par::Vector{Float64}; nbins::Union{Float64, Int} = NaN,
    ylbl::String = "Occurrence", xlbl::String = "Parameter", xunit = "", color = :silver, lbl::Union{Bool, String} = false,
    txt = true)

    if isnan(nbins)
        nbins = 0.25*length(par) < 20 ? 20 :  round(Int, 0.25*length(par)) 
    elseif isa(nbins, Float64)
        nbins = round(Int, nbins)
    end
    if isempty(xunit)
        xunit = "a.u."
    end

    # some plotting defaults 
    default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
    xguidefontsize = 16, xtickfontsize = 12,
    yguidefontsize = 16, ytickfontsize = 12,
    legendfontsize = 14, titlefontsize = 10,
    legendforegroundcolor = :silver)

    plt = stephist(par, nbins = nbins,
            ylabel = ylbl,
            xlabel = xlbl * " ($xunit)",
            label = lbl,
            fill =  true, 
            color = color,
            ylims = (0, :auto),
            margin = 3mm,
            dpi = 300)
    if txt == true 
        annotate!(xlims()[2] - 0.95*(xlims()[2] - xlims()[1]), 0.9*ylims()[2], 
        text(@sprintf("mean = %.1f %s", ustrip(mean(par)), xunit), :black, :left))
    end

    return plt
end

function plot_dsp_par(dsp::PropDict, par::Symbol; kwargs...)
    x = getproperty(dsp, par)
    plot_dsp_par(x; kwargs...)
end

function log_tick_formatter(val)
    exponent = floor.(Int, log10.(val))
    coefficient = round.(val ./ 10 .^ exponent, digits=1)
    return  [latexstring("\$ $(coefficient[i])\\cdot\\textrm{10}^\\textrm{$(exponent[i])}\$") for i in eachindex(val)]
end

function log_tick_formatter_plain(val)
    return [@sprintf("%.0f ",val[i]) for i in eachindex(val)]
end

function Makie_theme(; fs = 20)
        PlotTheme = Theme(
        size = (600, 420),
        fontsize = fs,
        dpi = 300,
        margin = 3mm,
        Fonts = (
            regular="Helvetica", math="Helvetica"),
        Axis = (
            titlefont = :regular,
            xgridvisible = false,
            ygridvisible = false,
            xlabelsize= fs + 6,
            ylabelsize= fs + 6,
            xtickalign=1,
            ytickalign=1,
        ),
        Legend = (
            framecolor = :silver,
            labelsize = fs + 2,
            position = :lt,
        )
    )
    set_theme!(PlotTheme)
end