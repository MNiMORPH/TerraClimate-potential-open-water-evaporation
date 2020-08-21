import numpy as np
import netCDF4 as nc4
import net_radiation
import atmospheric_parameters
import wind_shear_velocity

# Using TerraClimate
# 2.5 arcminute (1/24 degree) resolution: ~5 km N-S

# Import step
# ... load files here or with a CLI
years = range(1958, 2019)
months_zero_indexed = range(12)
TerraClimateDir = '/media/andy/data1/TerraClimate/'

def extract_data(varnc, varname, varmonth_zero_indexed=None):
    if varmonth_zero_indexed is None:
        var = varnc.variables[varname][:]
    else:
        var = varnc.variables[varname][varmonth_zero_indexed]
    fv = var.fill_value
    var = var.data
    var[var == fv] = np.nan
    return var

# Get lats and lons from one file
srad_nc = nc4.Dataset(TerraClimateDir+'TerraClimate_srad_1958.nc')
lats = extract_data(srad_nc, 'lat')
lons = extract_data(srad_nc, 'lon')
LONS, LATS = np.meshgrid (lons, lats)

# Shear velocity of winds: tool to compute from velocity
ustar_interp = wind_shear_velocity.create_lookup_table_one_step()

# Elevation
elevation_nc = nc4.Dataset(TerraClimateDir+'Gebco_2020_2_point_5_arcminute.nc')
elevation = extract_data(elevation_nc, 'value')
elevation = elevation[::-1]

# Heat capacity of air
specific_heat_capacity_of_air = 1.005 # approx. constant at 1 atm
                                      # Humidity minor impact below 40C or so
                                      # But this is an approximation!
cp = specific_heat_capacity_of_air # Easier

# Water density
rho_w = 1000.

# Latent heat of vaporization for water
Lv = 2.5E6
DeltaH_vap = Lv # to make me happier

# Ratio of molecular weight of water vapor to dry air
epsilon = 0.622

# Days in month, for weighting
days_in_month = [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# Evaporation array
E = np.zeros(elevation.shape)

for year in years[:1]:
    print(year)
    # Incoming solar radiation (monthly average)
    srad_nc = nc4.Dataset(TerraClimateDir+'TerraClimate_srad_'+str(year)+'.nc')
    # Maximum daily temperature (monthly average)
    tmax_nc = nc4.Dataset(TerraClimateDir+'TerraClimate_tmax_'+str(year)+'.nc')
    # Minimum daily temperature (monthly average)
    tmin_nc = nc4.Dataset(TerraClimateDir+'TerraClimate_tmin_'+str(year)+'.nc')
    # Wind speed (monthly average)
    ws_nc = nc4.Dataset(TerraClimateDir+'TerraClimate_ws_'+str(year)+'.nc')
    # Vapor pressure (monthly average)
    vap_nc = nc4.Dataset(TerraClimateDir+'TerraClimate_vap_'+str(year)+'.nc')

    # Now compute for each month
    for month in months_zero_indexed:
        print(month+1)
        # Data
        srad = extract_data(srad_nc, 'srad', month)
        tmax = extract_data(tmax_nc, 'tmax', month)
        tmin = extract_data(tmin_nc, 'tmin', month)
        ws = extract_data(ws_nc, 'ws', month)
        vap = extract_data(vap_nc, 'vap', month) * 1000. # kPa to Pa
        julian_day = np.sum(days_in_month[:month]) + days_in_month[month]/2.
        #elevation = 2000. # placeholder
        #julian_day = 205 # placeholder
        #vap = .03*101325 # placeholder
        albedo = 0.06

        # Calculations:
        # Net Radiation
        Rn = net_radiation.computeNetRadiation(elevation, julian_day, LATS,
                                                tmax, tmin, vap, srad, albedo)

        # Shear velocity of winds
        ustar = ustar_interp(ws)

        # Vapor-pressure deficit
        # We don't have max and min humidity
        VPD = atmospheric_parameters.compute_vpd( (tmax+tmin)/2., vap )

        # Atmospheric pressure
        P = atmospheric_parameters.compute_atmospheric_pressure(elevation)

        # Atmospheric density (ignoring temperature + humidity effects)
        rho_a = atmospheric_parameters.compute_atmospheric_density(elevation,
                                        (tmax + tmin)/2.)

        # Clausius-Clayperon phase-change slope
        Delta = ( atmospheric_parameters.compute_Delta_e_sat( tmax )
                  + atmospheric_parameters.compute_Delta_e_sat( tmin ) ) / 2.

        _E = (Rn + cp*rho_a*ustar**2/(Delta*ws) * VPD) \
             / ( rho_w*Lv  + P*cp*rho_w/epsilon )
        _E[_E<0] = 0 # ignore condensation; I think it's spurious (Antarctica?)
        E += _E*days_in_month[month]

# REMEMBER TO UPDATE FOR FULL TIME SERIES
E /= (365.25*len(years[:1]))
