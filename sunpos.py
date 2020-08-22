# Modified & expanded from https://github.com/Neon22/sun-angle

import numpy as np

twopi = 2 * np.pi
deg2rad = np.pi / 180.
solar_constant = 1367. # W/m^2

def calc_datetime( _datetime ):
    # JD - 2.4M - 51545 (latter for astronomer's almanac)
    return _datetime.toordinal() + _datetime.hour/24. \
           + _datetime.minute/60./24. + _datetime.second/60./60./24. \
           + 1721424.5 - 2400000 - 51545

def sun_position(_datetime, lats, lons):

    datetime_aa = calc_datetime( _datetime )

    # Ecliptic coordinates
    # Mean longitude
    mnlong = 280.46 + 0.9856474 * datetime_aa
    mnlong = mnlong % 360
    if mnlong < 0: mnlong += 360

    # Mean anomaly
    mnanom = 357.528 + 0.9856003 * datetime_aa
    mnanom = mnanom % 360
    if mnanom < 0: mnanom += 360
    mnanom = mnanom * deg2rad

    # Ecliptic longitude and obliquity of ecliptic
    eclong = mnlong + 1.915 * np.sin(mnanom) + 0.02 * np.sin(2 * mnanom)
    eclong = eclong % 360
    if eclong < 0: eclong += 360
    oblqec = 23.439 - 0.0000004 * datetime_aa
    eclong = eclong * deg2rad
    oblqec = oblqec * deg2rad

    # Celestial coordinates
    # Right ascension and declination
    num = np.cos(oblqec) * np.sin(eclong)
    den = np.cos(eclong)
    ra = np.arctan(num / den)
    if den < 0: ra += np.pi
    if den >= 0 and num < 0: ra += twopi
    dec = np.arcsin(np.sin(oblqec) * np.sin(eclong))

    # Local coordinates
    # Greenwich mean sidereal time
    gmst = 6.697375 + 0.0657098242 * datetime_aa + _datetime.hour
    gmst = gmst % 24
    if gmst < 0: gmst += 24

    # Local mean sidereal time
    lmst = gmst + lons / 15.0
    lmst = lmst % 24
    #if lmst < 0: lmst += 24
    lmst[lmst < 0] += 24
    lmst = lmst * 15 * deg2rad

    # Hour angle
    ha = lmst - ra
    #if ha < -np.pi: ha += twopi
    #if ha > np.pi: ha -= twopi
    ha[ha < -np.pi] += twopi
    ha[ha > np.pi] -= twopi

    # Latitude to radians
    lats = lats * deg2rad

    # Solar zenith angle
    zenithAngle = np.arccos(np.sin(lats) * np.sin(dec) + np.cos(lats) * np.cos(dec) * np.cos(ha))
    # Solar azimuth
    az = np.arccos(((np.sin(lats) * np.cos(zenithAngle)) - np.sin(dec)) / (np.cos(lats) * np.sin(zenithAngle)))
    # Solar elevation
    el = np.arcsin(np.sin(dec) * np.sin(lats) + np.cos(dec) * np.cos(lats) * np.cos(ha))

    el = el / deg2rad
    az = az / deg2rad
    lats = lats / deg2rad

    # Azimuth correction for Hour Angle
    #if ha > 0:
    #    az += 180
    #else:
    #    az = 540 - az
    #az = az % 360
    az[ha > 0] += 180
    az[ha <= 0] = 540 - az[ha <= 0]
    az %= 360

    return(az, el)

def extraterrestrial_solar_radiation(_datetime, lats, lons):
    az, el_pos = sun_position(_datetime, lats, lons)
    el_pos[el_pos < 0] = 0 # no negative solar radiation
    return solar_constant * np.sin(el_pos * deg2rad)
