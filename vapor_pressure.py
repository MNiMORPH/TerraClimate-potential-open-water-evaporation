import numpy as np

# Latent heat of vaporization for water
Lv = 2.5E6
# Specific gas constant of water vapor
Rv = 461. # J / (kg K)


# VAPOR PRESSURE

def compute_saturation_vapor_pressure(Tair_degC):
    Tair = Tair_degC + 273.15
    # e_sat
    return 611 * np.exp( Lv/Rv * (1/273.15 - 1/Tair))

def compute_vapor_pressure(e_sat, RH):
    # e_a
    return RH * e_sat

def compute_vpd(Tair_degC, e_a):
    e_sat = compute_saturation_vapor_pressure(Tair_degC)
    return e_sat - e_a


# SATURATION VAPOR PRESSURE (CLAUSIUS-CLAYPERON) SLOPE

def compute_Delta_e_sat(Tair_degC):
    # d(e_sat)/dT
    return 25083 / (Tair_degC + 237.3)**2 \
            * np.exp(17.3*Tair_degC / (Tair_degC +273.3 ))


# SPECIFIC HEAT CAPACITY OF AIR AT CONSTANT PRESSURE
# from Zotarelli et al. -- I'd like to double-check this (and LW rad approx)
def compute_atmospheric_pressure(elevation):
    return 101300 * np.exp( (293 - 0.0065*elevation)/293. )
