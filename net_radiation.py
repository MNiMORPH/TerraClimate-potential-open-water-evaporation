import numpy as np
import sunpos
import datetime

solar_constant = 1360.8 # W/m^2
stefan_boltzmann_constant = 5.67E-8
von_karman_constant = 0.407
wind_measurement_elevation = 2. # m

def _computeClearSkySolarRadiation(elevation, date, latitude, nlons):

    # Step 1: top-of-atmosphere ("extraterrestrial") radiation
    # Numerically integrate with dt = 1 minute
    time_start = datetime.time(0,0,0)
    _datetime = datetime.datetime.combine(date, time_start)
    _timedelta = datetime.timedelta(minutes=1) # dt = 1 minute
    # Assuming that it is a numpy array
    # Daily mean extraterrestrial radiation
    extraterrestrial_radiation = np.zeros(latitude.shape)
    while _datetime.date() == date:
        extraterrestrial_radiation += \
                        sunpos.extraterrestrial_solar_radiation(
                            _datetime, latitude, np.zeros(latitude.shape)
                            )
        _datetime += _timedelta
    extraterrestrial_radiation /= (24.*60.)

    # Same values across latitudes
    extraterrestrial_radiation = np.tile( np.expand_dims(extraterrestrial_radiation, 2), nlons)

    # Clear sky input a f(elevation)
    clear_sky_solar_radiation = (0.75 + 2E-5 * elevation) * extraterrestrial_radiation
    return clear_sky_solar_radiation

def computeNetLongwaveRadiation(elevation,
                                Tmax_degC, Tmin_degC, vapor_pressure,
                                incoming_solar_radiation,
                                clear_sky_solar_radiation):
    """
    Following Zotarelli et al., 2010
    If this is for the land surface, it may need some modification to be
    appropriate for lakes and their heating/cooling
    """
    isr = incoming_solar_radiation
    cssr = clear_sky_solar_radiation
    # 0.7 is from that 2017 paper I think. Anyway, seems reasonable.
    # Or, let's get a fill value from similar places
    cloud_function = np.nan * np.zeros(cssr.shape)
    # Cutoff before measured vs. modelled issues start to cause problems
    _cf_cutoff = 15.
    cloud_function[cssr>_cf_cutoff] = 1.35*isr[cssr>_cf_cutoff]/cssr[cssr>_cf_cutoff] - 0.35
    _cfm_near = np.nanmean(cloud_function[(cssr> _cf_cutoff) * (cssr< _cf_cutoff+5)])
    cloud_function[cssr<=_cf_cutoff] = _cfm_near

    # Looks like this is positive = loss -- near-surface air temperature
    # related to postive radiation.
    Rnl = stefan_boltzmann_constant \
              * ( (Tmax_degC + 273.16)**4 + (Tmin_degC + 273.16)**4 ) / 2. \
              * (.34 - .14E-3 * vapor_pressure**.5) \
              * cloud_function
    return -Rnl # Negative to note taht it is positive incoming

def computeNetShortwaveRadiation(incoming_solar_radiation, albedo):
    return incoming_solar_radiation * (1-albedo)

def computeNetRadiation(elevation, date, latitude, nlons,
                                Tmax_degC, Tmin_degC, vapor_pressure,
                                incoming_solar_radiation, albedo):
    cs = _computeClearSkySolarRadiation(elevation, date, latitude, nlons)
    lw = computeNetLongwaveRadiation(elevation,
                                Tmax_degC, Tmin_degC, vapor_pressure,
                                incoming_solar_radiation, cs)
    sw = computeNetShortwaveRadiation(incoming_solar_radiation, albedo)
    return lw + sw
