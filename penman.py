import numpy as np
from scipy.optimize import root_scalar, root
from scipy.special import lambertw
from scipy.interpolate import interp1d

# Using TerraClimate

solar_constant = 1360.8 # W/m^2
stefan_boltzmann_constant = 5.67E-8
von_karman_constant = 0.407
wind_measurement_elevation = 2. # m


def _computeClearSkySolarRadiation(elevation, julian_day, latitude):

    # Step 1: top-of-atmosphere ("extraterrestrial") radiation
    # Change 365 to number of days in year?
    inverse_relative_earth_sun_distance = 1 + 0.033 * np.cos(2*np.pi/365. * julian_day )
    solar_declination = 0.409 * np.sin( 2*np.pi/365. * julian_day - 1.39 )
    
    sunset_hour_angle = np.arccos( -np.tan(np.pi/180. * latitude) 
                                   * np.tan(solar_declination) )
    
    extraterrestrial_radiation = solar_constant * inverse_relative_earth_sun_distance \
                                 * ( sunset_hour_angle
                                     * np.sin(np.pi/180. * latitude)
                                     * np.sin(solar_declination)
                                     + np.cos(np.pi/180. * latitude)
                                     * np.cos(solar_declination)
                                     * np.sin(sunset_hour_angle) )

    clear_sky_solar_radiation = (0.75 + 2E-5 * elevation) * extraterrestrial_radiation
    return clear_sky_solar_radiation

def computeNetLongwaveRadiation(elevation, julian_day, latitude,
                                Tmax_degC, Tmin_degC, vapor_pressure,
                                incoming_solar_radiation,
                                clear_sky_solar_radiation):
    """
    Following Zotarelli et al., 2010
    If this is for the land surface, it may need some modification to be
    appropriate for lakes and their heating/cooling
    """
    Rnl = stefan_boltzmann_constant \
              * ( (Tmax_degC + 273.16)**4 + (Tmin_degC + 273.16)**4 ) / 2. \
              * (.34 - .14 * vapor_pressure**.5) \
              * (1.35*incoming_solar_radiation/clear_sky_solar_radiation - 0.35)
    return Rnl
    
def computeNetShortwaveRadiation(incoming_solar_radiation, albedo):
    return incoming_solar_radiation * (1-albedo)



# Hersbach 2011

def ustar_fcn(ustar, v):
    if ustar != 0:
        kappa = .407
        g = 9.805
        z = 2.
        ac = 1.8E-2
        am = 0.11
        nu_air = 1.5E-5
        z0_viscous_term = am * nu_air / ustar
        z0_charnock_turbulent_term = ac * ustar**2 / g
        z0 = z0_viscous_term + z0_charnock_turbulent_term
        # return v * kappa / np.log( g * z / (ac * ustar**2) ) - ustar
        return v * kappa / np.log( z / z0 ) - ustar
    else:
        return np.inf
    

def create_shear_velocity_lookup_table(v_array = np.arange(0, 55, 0.1)):
    ustar_array = []
    for v in v_array:
        rf = root_scalar(ustar_fcn, args=(v), x0=.5, x1=10.)
        ustar_array.append(rf.root)
    ustar_array = np.array(ustar_array)
    v_array = np.array(v_array)
    ustar_array[v_array == 0] = 1E-12 # safe tiny number
    return v_array, ustar_array

wind_velocities, shear_velocities = create_shear_velocity_lookup_table()
shear_velocities[0] = 0 # fix from root finder

def compute_z0_from_ustar(wind_velocities, shear_velocities):
        kappa = .407
        z = 2.
        return z / np.exp( kappa * wind_velocities / shear_velocities )

z0 = compute_z0_from_ustar(wind_velocities, shear_velocities)
plt.plot(wind_velocities, shear_velocities)


# SEEMS TO ONLY WORK FOR A SMALL RANGE OF VALUES, AND NOT WELL AT THAT.
# QUITE CONFUSED.
def solve_shear_velocity_analytical(v_array = np.arange(0, 50, 0.1)):
    ustar_array = []
    for v in v_array:
        ustar_array.append( -0.2035*v / lambertw(-0.01372 * np.abs(v)) )
    ustar_array = np.array(ustar_array)
    return v_array, ustar_array

wind_velocities_analytical, shear_velocities_analytical = solve_shear_velocity_analytical()


ustar_interp = interp1d(wind_velocities, shear_velocities)


#plt.plot(wind_velocities_analytical, shear_velocities_analytical)


def compute_shear_velocity(wind_speed):
    return ustar_interp(wind_speed)
    



