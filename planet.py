import astropy.constants as co
import astropy.units as u


def make_light_curve(cadence, tspan, noise, ground_based=False,
                     night_length=8.):
    """
    Generate a light curve with a specified cadence, time-span and white noise
    properties.

    Args:
        cadence (float): The time (in days) between observations.
        tspan (float): The total time span of observations in days.
        noise (float): The white noise. Assumed to be one sigma.
        ground_based (bool): True if ground based, False (default) if space-
            based. If true, day sections will be removed from light curve.
        night_length (float): If ground_based, the length of each night in
            hours. Default is 8.


    Returns:
        time (array): The time array.
        flux (array): The flux array.

    """
    time = np.arange(0, tspan, cadence)
    flux = np.ones_like(t) + np.random.randn(len(t)) * noise

    if ground_based:
        for night in range(tspan):
            mask = (time < night_length) * (night_length)

    return time, flux


def inject_transit(time, flux, period_days, radius tint=0)
    """
    Inject a planet transit into a light curve.

    Args:
        time (array): The time array in days.
        flux (array): The flux array.
        period_days (float): The planet's orbital period in days.
        radius (float): The planet's radius in units of stellar radii.
        tint (Optional[float]): The integration time (in minutes).

    Returns:
        new_flux (array): The flux array with injected planet transit.

    """

    period_s = period_days * 24 * 3600 * u.s
    G = 6.67e-11
    M = 2e30
    a = (period_s**2*co.G*co.M_sun/(4*np.pi**2))**(1./3)/co.R_sun

    params = batman.TransitParams()  # object to store transit parameters
    params.t0 = np.random.uniform(0, 20)  # time of inferior conjunction
    params.per = period_days  # orbital period
    params.rp = radius  # planet radius (in units of stellar radii)
    params.a = a  # semi-major axis (in units of stellar radii)
    params.inc = 90.  # orbital inclination (in degrees)
    params.ecc = 0.  # eccentricity
    params.w = 90.  # longitude of periastron (in degrees)
    params.limb_dark = "nonlinear"  # limb darkening model
    params.u = [0.5, 0.1, 0.1, -0.1]  # limb darkening coefficients

    m = batman.TransitModel(params, time) #initializes model
    new_flux = flux + m.light_curve(params)  #calculates light curve

    return new_flux
