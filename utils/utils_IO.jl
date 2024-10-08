using HDF5
using CSV, DataFrames
using Unitful
using RadiationDetectorSignals
# """
# OLD
#     read_wvf_csv(folder::String, file::String; heading = 2)
# read 1 csv file and return a RDWaveform
# input:
# folder::String: folder name. relative to ENV["ASIC_DATA"]
# file::String: file name
# """
# function read_wvf_csv(folder::String, filename::String;  heading = 2)
#     fpath = get(ENV, "ASIC_DATA",".") * joinpath(folder, filename)
#     df = CSV.read(fpath, DataFrame; delim=',', header = heading)
#     asic_voltages = df[:, "ASIC Voltage (V)"]
#     csv_time = uconvert.(u"µs", (df[:, "Time (s)"] .- df[1, "Time (s)"])*u"s")
#     timestep = csv_time[2] - csv_time[1]
#     time = 0u"µs":timestep:(length(asic_voltages) - 1)*timestep
#     RDWaveform(time, asic_voltages)
# end

"""
    read_csv_metadata(filepath::String; heading::Int = 17, nChannels::Int = 2)
read metadata from a csv file (as taken from oscilloscope)
"""
function read_csv_metadata(filepath::String; heading::Int = 17, nChannels::Int = 2)
    header_fileds = CSV.read(filepath, DataFrame; delim=',', limit = heading-1, header = false)
    MetaData1_lbls = collect(header_fileds.Column1)
    MetaData1_vals = collect(header_fileds.Column2)
    IdxKeep = map(x -> !ismissing(x), MetaData1_vals)
    MetaData = Dict("Ch1" => Dict(zip(String.(MetaData1_lbls[IdxKeep]), String.(MetaData1_vals[IdxKeep]))))
    timestep = uconvert(u"µs", parse(Float64, MetaData["Ch1"]["Sample Interval"]) .* uparse(MetaData["Ch1"]["Horizontal Units"]))
    if nChannels == 2
        MetaData2_lbls = collect(header_fileds.Column4)
        MetaData2_vals = collect(header_fileds.Column5)
        IdxKeep = map(x -> !ismissing(x), MetaData2_vals)
        MetaData["Ch2"] = Dict(zip(String.(MetaData2_lbls[IdxKeep]), String.(MetaData2_vals[IdxKeep])))
    end 
    return MetaData, timestep 
end


function read_folder_csv_oscilloscope(folder::String; heading::Int = 17, nwvfmax::Union{Int, Float64} = NaN, nChannels::Int = 2)
    files = readdir(get(ENV, "ASIC_DATA",".") * folder)
    filter!(x -> occursin(".csv", x), files)
    if !isnan(nwvfmax)
        files = files[1:nwvfmax]
    end

    fnames = [get(ENV, "ASIC_DATA",".") * joinpath(folder, file) for file in files]
    data = [CSV.read(fname, DataFrame; delim=',', header = heading) for fname in fnames]

    # Get MetaData and timestep 
    MetaData, timestep = read_csv_metadata(fnames[1], heading = heading, nChannels = nChannels)
    data_ch1 = [df[!,  MetaData["Ch1"]["Channel"]] for df in data]

    # make sure than vertical unit is always in Volts!
    if uparse(MetaData["Ch1"]["Vertical Units"]) != u"V"
        data_ch1 = [ustrip.(uconvert.(u"V", wvf .* uparse(MetaData["Ch1"]["Vertical Units"]))) for wvf in data_ch1]
    end 

    if nChannels == 2
        data_ch2 = [df[!, MetaData["Ch2"]["Channel"]] for df in data]
        if uparse(MetaData["Ch2"]["Vertical Units"]) != u"V"
            data_ch2 = [ustrip.(uconvert.(u"V", wvf .* uparse(MetaData["Ch2"]["Vertical Units"]))) for wvf in data_ch2]
        end 
    end 

    times = fill(0u"µs":timestep:(length(data_ch1[1]) - 1)*timestep, length(data_ch1[1]))
    wvfs_ch1 = ArrayOfRDWaveforms([RDWaveform(time, wvfs) for (time, wvfs) in zip(times, data_ch1)])

    if nChannels == 1
        return wvfs_ch1,  MetaData
    elseif nChannels ==2 
        wvfs_ch2 = ArrayOfRDWaveforms([RDWaveform(time, wvfs) for (time, wvfs) in zip(times, data_ch2)])
        return wvfs_ch1, wvfs_ch2, MetaData
    end 
end

function read_file_csv_oscilloscope(filepath::String;  heading::Int = 17, nChannels::Int = 2)
    data = CSV.read(filepath, DataFrame; delim=',', header = heading)
    MetaData, timestep = read_csv_metadata(fnames[1], heading = heading, nChannels = nChannels)
    data_ch1 = data[!, MetaData["Ch1"]["Channel"]]
    # Get MetaData and timestep 
    if nChannels == 2
        data_ch2 = data[!, MetaData["Ch2"]["Channel"]]
    end 
    time = fill(0u"µs":timestep:(length(data_ch1) - 1)*timestep, length(data_ch1))
    wvfs_ch1 = RDWaveform(time, data_ch1) 
    if nChannels == 1
        return wvfs_ch1, MetaData
    elseif nChannels ==2 
        wvfs_ch2 = RDWaveform(time, data_ch2) 
        return wvfs_ch1, wvfs_ch2, MetaData
    end 
end


### HDF5
function read_wvfs_h5(filename::String; folder::String = "", nwvfs::Union{Int, Float64} = NaN, fixtime::Bool = false)
    if endswith(filename, ".hdf5")
        filename = replace(filename, ".hdf5" => "")
    end

    # open file
    fname  = get(ENV, "ASIC_DATA",".") * folder * filename * ".hdf5"
    fopen = h5open(fname,"r")

    # read metadata
    metadata = Dict{String, Any}()
    for k in keys(attributes(fopen))
        metadata[string(k)] = read(attributes(fopen), k)
    end

    # # read pulses
    # data_group = fopen[replace(filename, "run" => "test")]
    # nwvfs = ifelse(isnan(nwvfs), length(keys(data_group)), nwvfs)
    # voltages = [read(data_group[w])[2,:] for w in keys(data_group)[1:nwvfs]]
    # times_abs = [uconvert.(u"µs", read(data_group[w])[1,:]*u"ms") for w in keys(data_group)[1:nwvfs]]
    # timesteps = [t[2] - t[1] for t in times_abs]
    # if fixtime == true
    #     times = fill(0u"µs":timesteps[1]:(length(voltages[1]) - 1)*timesteps[1], length(voltages))
    # else
    #     times = [0u"µs":timesteps[i]:(length(voltages[i]) - 1)*timesteps[i] for i in eachindex(voltages)] 
    # end
    signal = fopen["signal"]
    nwvfs = ifelse(isnan(nwvfs), length(keys(signal)), nwvfs)
    voltages_unit = read(attributes(signal[keys(signal)[1]]), "Unit")
    voltages = [read(signal[w]) for w in keys(signal)[1:nwvfs]]
    voltages_unit = read(attributes(signal[keys(signal)[1]]), "Unit") 

    time = fopen["timestep"][keys(signal)[1]]
    timestep_unit = uparse(read(attributes(time), "Unit"))
    timestep = uconvert.(u"µs", read(time) * timestep_unit)
    times = [0u"µs":timestep:(length(voltages[i]) - 1)* timestep for i in eachindex(voltages)] 

    return ArrayOfRDWaveforms((times, voltages)), metadata 
end 
