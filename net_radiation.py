import numpy as np

solar_constant = 1360.8 # W/m^2
stefan_boltzmann_constant = 5.67E-8
von_karman_constant = 0.407
wind_measurement_elevation = 2. # m

def _computeClearSkySolarRadiation(elevation, julian_day, latitude):

    # Step 1: top-of-atmosphere ("extraterrestrial") radiation
    # Change 365 to number of days in year?
    inverse_relative_earth_sun_distance = 1 + 0.033 * np.cos(2*np.pi/365. * julian_day )
    #solar_declination = 0.409 * np.sin( 2*np.pi/365. * julian_day - 1.39 )
    # Rosenberg 1983 -- basically the same answer.
    solar_declination = 0.4101 * np.cos( 2*np.pi/365. * (julian_day - 172) )

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

def computeNetRadiation(elevation, julian_day, latitude,
                                Tmax_degC, Tmin_degC, vapor_pressure,
                                incoming_solar_radiation, albedo):
    cs = _computeClearSkySolarRadiation(elevation, julian_day, latitude)
    lw = computeNetLongwaveRadiation(elevation, julian_day, latitude,
                                Tmax_degC, Tmin_degC, vapor_pressure,
                                incoming_solar_radiation, cs)
    sw = computeNetShortwaveRadiation(incoming_solar_radiation, albedo)
    return lw + sw
