using Measurements
using Measurements: value as mvalue
using Measurements: uncertainty as muncert

function benchtest_linearity(e_trap_V::Union{Vector{Float64}, Vector{Measurement{Float64}}}, 
        pulser_trap_V::Vector{Float64}, set, pd_data::PropDict; VinMode::String = "MeV",
        VoutMode::String = "V", maxCharge_keV::Union{Int64, Float64} = NaN)

    if !isnan(maxCharge_keV)
        pulser_trap_keV = 1e-3 * csa_voltage2charge.(set.Cinj_fF, pulser_trap_V) * EnergyPerEHpair_Ge_eV(77)
        KeepIdx = findall(x-> x <= maxCharge_keV, pulser_trap_keV)
        pulser_trap_V = pulser_trap_V[KeepIdx]
        e_trap_V = e_trap_V[KeepIdx]
    end

    if isa(e_trap_V, Vector{Measurement{Float64}})
        y = measurement.(e_trap_V, yerr)
        e_trap_V_err = muncert.(e_trap_V)
        e_trap_V = mvalue.(e_trap_V)
    end

    # prepare for fit
    if VinMode == "V"
        x = 1e3 .* pulser_trap_V
        xlbl = L"$V_\textrm{in}\,$" 
        xustr = "mV"
    elseif VinMode == "e"
        x = 1e-6 .* csa_voltage2charge.(set.Cinj_fF, pulser_trap_V) 
        xlbl = L"$Q_\textrm{in}\,$" 
        xustr = L"$\textrm{M}e^-$"
    elseif VinMode == "eV"
        x = 1e-6 .* csa_voltage2charge.(set.Cinj_fF, pulser_trap_V) * EnergyPerEHpair_Ge_eV(77)
        xlbl = L"$Q_\textrm{in}\,$" 
        xustr = "MeV"
    end
    if VoutMode == "V"
        y =  1e3 .* e_trap_V
        if @isdefined(e_trap_V_err)
            yerr =  1e3 .* e_trap_V_err 
        end
        ylbl =  L"$V_\textrm{out}\,$" 
        yustr = "mV"
        yresustr = "mV"
    elseif VoutMode == "e"
        y = 1e-6 .* csa_voltage2charge.(set.Cf_fF,  e_trap_V)  
        if @isdefined(e_trap_V_err)
            yerr =  1e-6 .* csa_voltage2charge.(set.Cf_fF,  e_trap_V_err)
        end
        ylbl = L"$V_\textrm{out}\,$"
        yustr = L"$\textrm{M}e^-$"#L"$\, 10^6 \cdot e^-\,$"
        yresustr = L"$\textrm{k}e^-$"
    elseif VoutMode == "eV"
        y = 1e-6 .* csa_voltage2charge.(set.Cf_fF,  e_trap_V)   * EnergyPerEHpair_Ge_eV(77)
        if @isdefined(e_trap_V_err)
            yerr =  1e-6 .* csa_voltage2charge.(set.Cf_fF,  e_trap_V_err)  * EnergyPerEHpair_Ge_eV(77)
        end
        ylbl = L"$V_\textrm{out}\,$"
        yustr = "MeV"
        yresustr = "keV"
    end

    if @isdefined(y_err)
        result, report = chi2fit(1, x, measurement.(y, yerr ) ; v_init = [1, 1])
        result.gof.converged ? @info("Fit converged") : nothing
    else
        result, report = chi2fit(1, x, y ; v_init = [1, 1])
    end

    @sprintf("Result: slope:  (%.2f ± %.0e), offset (%.2f ± %.0e) %s ", result.par[2].val, result.par[2].err, 1e3 * result.par[1].val, 1e3 * result.par[1].err, yustr)
    if result.par[1].val > 0 
        fitlabel = @sprintf("Fit: (%.2f ± %.0e) ", result.par[2].val, result.par[2].err) * L"$\cdot\,$" * xlbl* @sprintf("+ %.1g %s", result.par[1].val, yustr)
    else
        fitlabel = @sprintf("Fit: (%.2f ± %.0e) ", result.par[2].val, result.par[2].err) * L"$\cdot\,$" * xlbl * @sprintf("- %.1g %s", abs(result.par[1].val), yustr)
    end

    (csaname, board) = CSAname(pd_data)

     # plot defaults 
     default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
     xguidefontsize = 16, xtickfontsize = 12,
     yguidefontsize = 16, ytickfontsize = 12,
     legendfontsize = 14, titlefontsize = 10,
     legendforegroundcolor = :silver,) 
     yfit = Measurements.value.(report.f_fit(x))
     plin = plot(x, yfit, 
         color = :red2, linewidth = 3,
         label = fitlabel)
     ms = 4 
     scatter!(x, y,
         size = (700, 450),
         markerstrokewidth=0,
         markersize = ms,
         top_margin = 2mm,
         linewidth = 1,
         legendfontsize = 12,
         color = :black,
         marker = :circle,
         label = "Bench test: $(csaname)-Board $(board)" ,
         ylabel = ylbl * "($yustr)")
     xl = xlims()
 
     pres = hline([0], color = :darkgrey, linewidth = 2, linestyle = :dash, label = false)
     if VoutMode == "eV" 
        hspan!([-1, 1], color = :silver, alpha = 0.8, label = false)
        hline!([0], color = :darkgrey, linewidth = 2, linestyle = :dash, label = false)
        # hline!([-1], color = :darkgrey, linewidth = 2, linestyle = :dash, label = false)
        # hline!([1], color = :darkgrey, linewidth = 2, linestyle = :dash, label = false)
     end

     if VoutMode == "eV" || VoutMode == "e"
        yres = (y .- yfit) .* 1e3
    else
        yres = y .- yfit
    end
     scatter!(x,  yres ,
             xlabel = "Injected pulse " * xlbl * "($xustr)", 
             ylabel = "Residuals ($(yresustr))",
             color = :black,
             markerstrokewidth = 0,
             markersize = ms,
             label = false, 
             xlims = xl)
 
     resrange = abs(diff([minimum(yres), maximum(yres)])[1])
     yl = maximum(abs.(yres)) + 0.2 * resrange
     ylims!(-yl, yl)
 
     layout = @layout [a{0.65h}; b]
     p99 = plot(plin, pres, layout = layout, size = (600, 600),
         plot_title = title = "CSA Linearity $(csaname)-Board $(board)\n" * L"$C_\textrm{inj}$" * " = $(set.Cinj_fF) fF, R = $(set.Rf_MOhm) MOhm, T = $(set.Temp_K) K", 
         plot_titlefontsize = 10, 
         top_margin = 1mm,
         right_margin = 3mm,
         bottom_margin = -2mm)
    display(p99)
    return p99
end