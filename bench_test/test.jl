Pkg.activate("$(@__DIR__)/../")
import Pkg
using TypedTables, Dates
using RadiationDetectorSignals
using Unitful, Formatting, LaTeXStrings, Measures, Measurements
using Measurements: value as mvalue, uncertainty as muncertainty
using Plots, StatsBase, PropDicts
using LegendHDF5IO, LegendDSP, LegendSpecFits
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP: get_fltpars
using ArraysOfArrays
using Distributed
using RadiationDetectorDSP
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendHDF5IO, LegendDSP, LegendSpecFits, LegendDataTypes
using IntervalSets, TypedTables, StatsBase, PropDicts
using Unitful, Formatting, Printf, Measures, Dates
using Distributed, ProgressMeter, TimerOutputs
using RadiationDetectorDSP
using ArraysOfArrays
using HDF5
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using LegendDSP: get_fltpars
using TypedTables, Dates
using Measurements: value as mvalue, uncertainty as muncertainty
using LegendDataManagement: readlprops, writelprops
using CSV


times = Vector{Vector{typeof(0.0u"µs":1.0u"µs":1.0u"µs")}}()
voltages = Vector{Vector{Float64}}()
decay_times = Vector{Float64}()
wvfarray = ArrayOfRDWaveforms

@__DIR__
path = "$(@__DIR__)/../../ASIC_data/ASIC_test_4p500V.csv"
file = CSV.File(path; header = 3)
Table(file)

csv_voltage = file["ASIC Voltage (V)"]
csv_time = uconvert.(u"µs", (file["Time (s)"] .- file["Time (s)"][1])*u"s")
csv_pulser_voltage = file["Pulser Voltage (V)"]
## formatting time properly for the RDWaveform object, then appending 
timestep = csv_time[2] - csv_time[1]
time = 0u"µs":timestep:(length(csv_voltage) - 1)*timestep

waveform = RDWaveform(time, csv_voltage)
waveform

push!(voltages, csv_voltage)
push!(times, time)
wvfsarray = ArrayOfRDWaveforms

push!(wvfsarray, waveform)

plot(time, csv_voltage, label = "ASIC & Pulser Voltage vs. Time", xlabel = "Time", ylabel = "Voltage (V)", title = "ASIC & Pulser Voltage vs. Time", legend=:right)
plot!(time, csv_pulser_voltage, label = "Pulser Voltage vs. Time")