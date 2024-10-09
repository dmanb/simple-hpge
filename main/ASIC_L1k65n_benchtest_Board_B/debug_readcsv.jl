import Pkg
Pkg.activate("$(@__DIR__)/../../")
using Unitful
using LegendDataManagement
using LegendDataManagement: readlprops
using TypedTables, StatsBase, PropDicts
using ArraysOfArrays
using Measures 
using Plots
using Printf, LaTeXStrings
include("$(@__DIR__)/../../utils/utils_IO.jl")
include("$(@__DIR__)/../../utils/utils_naming.jl")

# settings 
PulserChargeInj_keV =  100 # pulser charge injected in keV
Cinj_fF = 500 # capacitance of the ASIC in femto Farrad 
Rf_MOhm = 500 # feedback resistor in Mega Ohm
Temp_K = 300 # temperature in Kelvin 

# get dsp configuration
dpath = "$(@__DIR__)/data_config.json"
pd_data = readlprops(dpath)

# read data: doesnt work for this configuration
folder = benchtest_filename_raw(pd_data, Cinj_fF, Rf_MOhm, PulserChargeInj_keV, Temp_K) 
wvfs_raw, _ = read_folder_csv_oscilloscope(folder; heading = 17, nChannels = 1);

###### debug step by step:
nwvfmax = NaN
heading = 17
files = readdir(get(ENV, "ASIC_DATA",".") * folder)
filter!(x -> occursin(".csv", x), files)
if !isnan(nwvfmax)
    files = files[1:nwvfmax]
end

fnames = [get(ENV, "ASIC_DATA",".") * joinpath(folder, file) for file in files]
data = [CSV.read(fname, DataFrame; delim=',', header = heading) for fname in fnames]

# Get MetaData and timestep 
nChannels = 2
MetaData, timestep = read_csv_metadata(fnames[1], heading = heading, nChannels = nChannels)
data_ch1 = [df[!,  MetaData["Ch1"]["Channel"]] for df in data if any(map(x-> occursin(MetaData["Ch1"]["Channel"],x), names(df))) ]


for (i, df) in enumerate(data)
    try 
        df[!,  MetaData["Ch1"]["Channel"]]
    catch e
        println("error for file $i ($(fnames[i])): $e")
    end
end
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