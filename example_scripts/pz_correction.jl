import Pkg
Pkg.activate("$(@__DIR__)/../")
using Unitful
using LegendDSP
using RadiationDetectorDSP
using Statistics
using LegendDataManagement
using LegendDataManagement: readlprops
using Plots
include("$(@__DIR__)/../utils/utils_IO.jl")

# get dsp configuration
dsp_config_path = "$(@__DIR__)/../config/dsp_config.json"
dsp_config = DSPConfig(readlprops(dsp_config_path).default)

# some plotting defaults 
default(size=(600, 400), legend = :best, grid=:off, frame=:semi,
         xguidefontsize = 16, xtickfontsize = 12,
         yguidefontsize = 16, ytickfontsize = 12,
         legendfontsize = 14, titlefontsize = 14,
         legendforegroundcolor = :silver)

# read data 
folder = "ASIC_data_07312024/ASIC_PulserVoltage_0p400V/"
wvfs = read_folder_csv(folder; heading = 2)
plot(wvfs[1].time, wvfs[1].signal, label = "Raw waveform", ylabel = "Volts (V)", xlabel = "Time")

# extract decay times 
if PulserVoltage == 0.4
    wvfs = wvfs[1:end-1]
end
decay_times = dsp_decay_times(wvfs,dsp_config)
stephist(decay_times, bins=20, xlabel= "Decay time", ylabel="Occurrence", label = false, fill = true, color = :silver, ylims = (0, :auto))

# do pole-zero correction using median decay time for all waveforms
deconv_flt = InvCRFilter(median(decay_times))
# # wvfs_pz = ArrayOfRDWaveforms([deconv_flt[i](wvfs[i]) for i in eachindex(wvfs)])
wvfs_pz = deconv_flt.(wvfs)
plot(wvfs[1].time, wvfs[1].signal, label = "raw waveform", ylabel = "Volts (V)", xlabel = "Time")
plot!(wvfs_pz[1].time, wvfs_pz[1].signal, label = "pole-zero corrected",  legend = :right)

