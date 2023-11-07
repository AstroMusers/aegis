import numpy as np


class Exoplanet:
    def __init__(self, name, semi_major_axis, radius, mass, density, magnetic_field, star_mass, star_mass_loss,
                 distance):
        self.name = name
        self.semi_major_axis = semi_major_axis
        self.radius = radius
        self.mass = mass
        self.density = density
        self.magnetic_field = magnetic_field
        self.star_mass = star_mass
        self.star_mass_loss = star_mass_loss
        self.distance = distance

    def __repr__(self):
        return f"Exoplanet {self.name} with: a={self.semi_major_axis}, Bs={self.magnetic_field}, M={self.star_mass}," \
               f" Mdot={self.star_mass_loss}, D={self.distance} \n"


def max_freq(B):
    """
    Finds the maximum CMI emission frequency based on Jupiter's emission, as given in Ashtari et al. 2022. or Lynch et al. 2018
    :param B: Magnetic Field Strength at the surface of the exoplanet, must be given in Gauss.
    :return: Maximum CMI emission frequency of the exoplanet, given in Hz
    """
    # Bj = 14  # gauss
    # nuJ = 24 * 10 ** 6  # Hz
    # return nuJ * B / Bj
    return 2.8 * 10**6 * B

def density(p):
    """
    Returns the density of a planet in Jovian density
    :param p: density in g/cm^3
    :return: Density in p_j
    """
    return p / 1.326


def rotation(T, a):
    """
    Returns the planetary rotation rate in Jovian rotation rate. Assume tidally locked if orbit semi major axis is less
    than 0.1 AU.
    :param T: Planet orbital period in days
    :param a: Semi major axis in AU
    :return: Planetary rotation rate in Jovian rotation rate.
    """
    if a <= 0.1:
        return T / 0.414
    else:
        return 1


def convective_radius(exo):
    """
    Calculates the radius of the convective core of the exoplanet with the scaling provided in Curtis & Ness 1986.
    If the calculation exceeds planetary radius, planetary radius is returned. JOvain convective radius of 0.830 R_j
    is taken from Sharan 2022.
    :param exo: the exoplanet
    :return: Radius of the convective core
    """
    M = exo.mass
    p = exo.density
    r = exo.radius
    if p < 1.6 and M > 0.4:
        R = M ** 0.44  # Curtis & Ness (1986)
    else:
        R = M ** 0.75 * r ** (-0.96)  # Grießmeier (2004)
    return min(R, exo.radius) / 0.9


def magnetic_moment(p_c, w_p, r_c, sigma):
    """
    Finds the magnetic moment of the exoplanet using the density of its convective core, rotation rate, radius of its
    convective core, and the conductivity of its convective core. Results is in Jovian magnetic moments when all parameter
    values are given with respect to Jupiter.
    :param p_c: Density of the convective core.
    :param w_p: Planetary rotation rate.
    :param r_c: Radius of the convective core
    :param sigma: Conductivity of the convective core
    :return: Magnetic moment of the exoplanet.
    """
    return p_c ** (1 / 2) * w_p ** (1 / 2) * r_c ** 3 * sigma ** (-1 / 2)


def magnetic_field(mu, R):
    """
    Calculates the magnetic field strength at the surface of an exopplanet from its magnetic moment and radius. Given
    in Ashtari2022
    :param mu: Magnetic moment in A m^2
    :param R:Radius in Jupiter radii.
    :return: Magnetic field strength at the surface.
    """
    B_j = 4.17  # Gauss
    return mu / R ** 3 * B_j


def magnetopause(B, a):
    """
    Finds the magnetopause standoff distance of an exoplanet using its magnetic field strength
    and semi-major axis. The scaling equation in Ashtari et al. 2022 is used
    :param B: Magnetic field strength given in Gauss.
    :param a: Semi-major axis given in AU.
    :return: The magnetopause standoff distance of the exoplanet given in Jupiter radii.
    """
    Bj = 4.17  # gauss
    aj = 5.204  # AU
    Rj = 1  # Rj
    return 56 * Rj * ((B / Bj) / (aj / a)) ** (1 / 3)


def radio_lum(Rmp, a, M, Mdot):
    """
    Calculates the total output power of CMI emission of an exoplanet using its semi-major axis,
    The host star's mass and mass loss rate. The equation is given in Ashtari et al. 2022
    :param Rmp: Magnetopause standoff distance given in Jupiter Radii.
    :param a: Semi-major axis given in AU
    :param M: Host star's mass in solar masses.
    :param Mdot: Host star's mass loss rate in 10^-15 solar masses per year.
    :return: The total output power of the CMI emission, or the radio brightness of the exoplanet
            given in Watts.
    """
    Rmpj = 60  # Rj
    M_sun = 1  # M_sun
    M_sun_loss = 68  # e-15 solar masses
    Lj = 2.1 * 10 ** 11
    aj = 5.204  # AU
    return (Rmp / Rmpj) ** 2 * (aj / a) ** 2 * (M / M_sun) * (Mdot / M_sun_loss) * Lj


def radio_brightness(L, nu, D):
    """
    Calculates the beamed radio brightness of an exoplanet using its radio output power,
    emission frequency and observation distance.
    :param L: Total radio output power given in Watts.
    :param nu: Emission frequency of the exoplanet given in Hz.
    :param D: Distance to exoplanet given in meters.
    :return: Spectral Radio Flux density or Radio Brightness of the exoplanet in Janskies(Jy).
    """
    delta = 10.5
    return delta * L / (4 * np.pi * nu * D ** 2) * 10 ** 26


def complete(B, a, M, Mdot, D):
    """
    Finds the CMI radio brightness of the exoplanet directly using initial parameters.
    :param B: Magnetic field strength given in Gauss.
    :param a: Semi-major axis given in AU
    :param M: Host star's mass in solar masses.
    :param Mdot: Host star's mass loss rate in 10^-15 solar masses per year.
    :param D: Distance to exoplanet given in meters.
    :return: Spectral Radio Flux density or Radio Brightness of the exoplanet in Janskies(Jy).
    """
    nu = max_freq(B)
    Rmp = magnetopause(B, a)
    L = radio_lum(Rmp, a, M, Mdot)
    I = radio_brightness(L, nu, D)
    return I
