
using Measurements 
"""
    csa_voltage2charge(Cinj_fF::Int, Voltage_V::Float64)
converts voltage in V to charge in e based on capacitance Cinj_fF in fF
"""
function csa_voltage2charge(Cinj_fF::Int, Voltage_V::Union{Float64, Measurements.Measurement}) 
    charge_C = Voltage_V * (Cinj_fF * 1e-15)
    Coulomb2eV = 1 / (1.602 * 1e-19)
    charge_e = charge_C * Coulomb2eV
    return charge_e
end

"""
    bandgap energy of Germanium and energy per electron-hole pair in Germanium
    formulas from Frank Edzards PhD thesis 
"""
Bandgap_Ge_eV = (T) -> 0.774 - 4.77e-4 * T^2 / (T + 235)  # temperature T in Kelvin , returns bandgap energy of Germanium in eV
EnergyPerEHpair_Ge_eV = (T) -> 2.2 * Bandgap_Ge_eV(T)  + 1.99 * Bandgap_Ge_eV(T)^(3/2) * exp(4.75 * Bandgap_Ge_eV(T)/T)

