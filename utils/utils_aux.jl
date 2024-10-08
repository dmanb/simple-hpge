

"""
    csa_voltage2charge(Cinj_fF::Int, Voltage_V::Float64)
converts voltage in V to charge in eV based on capacitance Cinj_fF in fF
"""
function csa_voltage2charge(Cinj_fF::Int, Voltage_V::Float64) 
    charge_C = Voltage_V * (Cinj_fF * 1e-15)
    Coulomb2eV = 1 / (1.602 * 1e-19)
    charge_eV = charge_C * Coulomb2eV
    return charge_eV
end

