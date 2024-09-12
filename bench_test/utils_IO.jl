using HDF5
using CSV, DataFrames
using Unitful
using RadiationDetectorSignals
"""
    read_wvf_csv(folder::String, file::String; heading = 2)
read 1 csv file and return a RDWaveform
input:
folder::String: folder name. relative to ENV["ASIC_DATA"]
file::String: file name
"""
function read_wvf_csv(folder::String, filename::String;  heading = 2)
    fpath = get(ENV, "ASIC_DATA",".") * joinpath(folder, filename)
    df = CSV.read(fpath, DataFrame; delim=',', header = heading)
    asic_voltages = df[:, "ASIC Voltage (V)"]
    csv_time = uconvert.(u"µs", (df[:, "Time (s)"] .- df[1, "Time (s)"])*u"s")
    timestep = csv_time[2] - csv_time[1]
    time = 0u"µs":timestep:(length(asic_voltages) - 1)*timestep
    RDWaveform(time, asic_voltages)
end

"""
    read_folder_csv(folder::String; kwargs...)
read all csv files in a folder and return an ArrayOfRDWaveforms
folder is relative to ENV["ASIC_DATA"]
"""
function read_folder_csv(folder::String; kwargs...)
    # read all files
    files = readdir(get(ENV, "ASIC_DATA",".") * folder)
    filter!(x -> occursin(".csv", x), files)
    ArrayOfRDWaveforms(read_wvf_csv.(folder, files; kwargs...))
end

function read_wvfs_h5(filename::String; folder::String = "", nwvfs::Int = NaN)
    fname  = get(ENV, "ASIC_DATA",".") * folder *filename
    fopen = h5open(fname,"r")
    nwvfs = ifelse(isnan(nwvfs), length(keys(fopen)), nwvfs)
    voltages = [read(fopen[w])[2,:] for w in keys(fopen)[1:nwvfs]]
    times_abs = [uconvert.(u"µs", read(fopen[w])[1,:]*u"s") for w in keys(fopen)[1:nwvfs]]
    timesteps = [t[2] - t[1] for t in times_abs]
    times = [0u"µs":timesteps[i]:(length(voltages[i]) - 1)*timesteps[i] for i in eachindex(voltages)] 
    ArrayOfRDWaveforms((times, voltages))
end 
