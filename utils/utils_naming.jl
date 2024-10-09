

function benchtest_filename_dsp(Cinj_fF::Int, Rf_MOhm::Int, PulserChargeInj_keV::Union{Int, Float64}, Temp_K::Int)
    dsp_file =  @sprintf("DSPpar_Cinj%04dfF_Rf%3dM", Cinj_fF, Rf_MOhm)
    if PulserChargeInj_keV >= 1000
        PulserChargeInj_MeV = PulserChargeInj_keV / 1000
        if isinteger(PulserChargeInj_MeV)
            dsp_file = dsp_file * @sprintf("_ChargeInj%dMeV", PulserChargeInj_MeV)
        else 
            dsp_file = dsp_file *  @sprintf("_ChargeInj%.1fMeV", PulserChargeInj_MeV) |> x -> replace(x, "." => "p")
        end
    else 
        dsp_file = dsp_file * @sprintf("_ChargeInj%dkeV", PulserChargeInj_keV)
    end 
    dsp_file = dsp_file * @sprintf("_Temp%03dK_", Temp_K) * ".json"
    return dsp_file
end

function benchtest_filename_raw(pd_data, Cinj_fF::Int, Rf_MOhm::Int, PulserChargeInj_keV::Union{Int, Float64}, Temp_K::Int)
    datafolder = pd_data.datafolder * @sprintf("Cinj%04dfF_Rf%3dM_Pulser/", Cinj_fF, Rf_MOhm)
    if PulserChargeInj_keV >= 1000
        PulserChargeInj_MeV = PulserChargeInj_keV / 1000
        if isinteger(PulserChargeInj_MeV)
            datafolder = datafolder * @sprintf("%dMeV", PulserChargeInj_MeV)
        else 
            datafolder = datafolder *  @sprintf("%.1fMeV", PulserChargeInj_MeV) |> x -> replace(x, "." => "p")
        end
    else
        datafolder = datafolder * @sprintf("%dkeV", PulserChargeInj_keV)
    end
    datafolder = datafolder * @sprintf("_%2dK/",Temp_K)
    return datafolder
end

function benchtest_filename_plot(Cinj_fF::Int, Rf_MOhm::Int, PulserChargeInj_keV::Union{Int, Float64}, Temp_K::Int)
    fname =  @sprintf("Cinj%04dfF_Rf%3dM", Cinj_fF, Rf_MOhm)
    if PulserChargeInj_keV >= 1000
        PulserChargeInj_MeV = PulserChargeInj_keV / 1000
        if isinteger(PulserChargeInj_MeV)
            fname = fname * @sprintf("_ChargeInj%dMeV", PulserChargeInj_MeV)
        else 
            fname = fname *  @sprintf("_ChargeInj%.1fMeV", PulserChargeInj_MeV) |> x -> replace(x, "." => "p")
        end
    else
        fname = fname * @sprintf("_ChargeInj%dkeV", PulserChargeInj_keV)
    end
    fname = fname * @sprintf("_Temp%03dK", Temp_K) * ".png"
    return fname
end

function benchtest_filename_plot(Cinj_fF::Int, Rf_MOhm::Int, Temp_K::Int)
    fname =  @sprintf("Cinj%04dfF_Rf%3dM", Cinj_fF, Rf_MOhm)
    fname = fname * @sprintf("_Temp%03dK", Temp_K) * ".png"
    return fname
end

function CSAname(pd_data)
    csaname = pd_data.datafolder[1:findfirst("_CSA",pd_data.datafolder)[1]-1]
    board = string(pd_data.datafolder[findfirst("Board_",pd_data.datafolder)[end]+1])
    return csaname, board
end