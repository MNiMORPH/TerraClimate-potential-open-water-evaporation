import numpy as np
from scipy.optimize import root_scalar, root
from scipy.special import lambertw
from scipy.interpolate import interp1d

# Shear velocity from the wind, including wave generation on the lake surface

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

#wind_velocities, shear_velocities = create_shear_velocity_lookup_table()
#shear_velocities[0] = 0 # fix from root finder

def compute_z0_from_ustar(wind_velocities, shear_velocities):
        kappa = .407
        z = 2.
        return z / np.exp( kappa * wind_velocities / shear_velocities )

#z0 = compute_z0_from_ustar(wind_velocities, shear_velocities)
#plt.plot(wind_velocities, shear_velocities)


# SEEMS TO ONLY WORK FOR A SMALL RANGE OF VALUES, AND NOT WELL AT THAT.
# QUITE CONFUSED.
def solve_shear_velocity_analytical(v_array = np.arange(0, 50, 0.1)):
    ustar_array = []
    for v in v_array:
        ustar_array.append( -0.2035*v / lambertw(-0.01372 * np.abs(v)) )
    ustar_array = np.array(ustar_array)
    return v_array, ustar_array

wind_velocities_analytical, shear_velocities_analytical = solve_shear_velocity_analytical()


#ustar_interp = interp1d(wind_velocities, shear_velocities)


#plt.plot(wind_velocities_analytical, shear_velocities_analytical)

def compute_shear_velocity__lookup_interp(wind_speed):
    return ustar_interp(wind_speed)




# An iterative non-zero-finding approach to getting ustar

def compute_ustar_using_z0(wind_velocity, z0):
        kappa = .407
        z = 2.
        ustar = kappa * wind_velocity / np.log(z/z0)
        return ustar

def compute_water_surface_z0_using_ustar(ustar):
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
        return z0
    else:
        return np.inf

def iteratively_solve_ustar_z0(wind_velocity, z0_0=1E-3, tolerance=1E-6):
    z0 = z0_0 # initial guess
    z0_old = 1E12
    ustar_old = 1E12
    z0_err = 1E12
    ustar_err = 1E12
    while (z0_err > tolerance) and (ustar_err > tolerance):
        ustar = compute_ustar_using_z0(wind_velocity, z0)
        z0 = compute_water_surface_z0_using_ustar(ustar)
        ustar_err = np.abs( (ustar_old - ustar) / np.mean((ustar_old, ustar)) )
        z0_err = np.abs( (z0_old - z0) / np.mean((z0_old, z0)) )
        z0_old = z0
        ustar_old = ustar
        #print( ustar, z0, ustar_err, z0_err)
    return ustar, z0

def create_shear_velocity_lookup_table_iterative(v_array = np.arange(0, 59.6, 0.1)):
    ustar_array = []
    z0_last = 1E-5 # initial guess, updates through range
    for v in v_array:
        #print (v)
        ustar, z0 = iteratively_solve_ustar_z0(v, z0_0 = z0_last)
        ustar_array.append(ustar)
        if np.isfinite(z0):
            z0_last = z0
    ustar_array = np.array(ustar_array)
    v_array = np.array(v_array)
    ustar_array[v_array == 0] = 0.# 1E-12 # safe tiny number
    return v_array, ustar_array

def linearly_interpolate_ustar(v_array, ustar_array):
    wind_velocities_iter, shear_velocities_iter = create_shear_velocity_lookup_table_iterative()
    # ustar_interp
    return interp1d(wind_velocities_iter, shear_velocities_iter, bounds_error=False, fill_value="extrapolate")

#plt.plot(wind_velocities, shear_velocities, linewidth=3)
#plt.plot(wind_velocities_iter, shear_velocities_iter)
#plt.xlabel('Wind velocity [m/s]')
#plt.ylabel('Shear velocity including roughness\nfeedback with surface waves [m/s]')

def create_lookup_table_one_step():
    # ustar_interp
    wind_velocities_iter, shear_velocities_iter = create_shear_velocity_lookup_table_iterative()
    return linearly_interpolate_ustar( wind_velocities_iter, shear_velocities_iter )

def compute_shear_velocity__lookup(wind_speed, ustar_interp):
    return ustar_interp(wind_speed)
