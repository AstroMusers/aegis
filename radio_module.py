import numpy as np


class Exoplanet:
    def __init__(self, name, semi_major_axis, magnetic_field, star_mass, star_mass_loss, distance):
        self.name = name
        self.semi_major_axis = semi_major_axis
        self.magnetic_field = magnetic_field
        self.star_mass = star_mass
        self.star_mass_loss = star_mass_loss
        self.distance = distance

    def __repr__(self):
        return f"Exoplanet {self.name} with: a={self.semi_major_axis}, Bs={self.magnetic_field}, M={self.star_mass}, Mdot={self.star_mass_loss}, D={self.distance} \n"


def max_freq(B):
    """
    Finds the maximum CMI emission frequency based on Jupiter's emission, as given in Ashtari et al. 2022.
    :param B: Magnetic Field Strength at the surface of the exoplanet, must be given in Gauss.
    :return: Maximum CMI emission frequency of the exoplanet, given in Hz
    """
    Bj = 14  # gauss
    nuJ = 24 * 10**6  # Hz
    return nuJ * B / Bj


def magnetopause(B, a):
    """
    Finds the magnetopause standoff distance of an exoplanet using its magnetic field strength
    and semi-major axis. The scaling equation in Ashtari et al. 2022 is used
    :param B: Magnetic field strength given in Gauss.
    :param a: Semi-major axis given in AU.
    :return: The magnetopause standoff distance of the exoplanet given in Jupiter radii.
    """
    Bj = 14  # gauss
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
