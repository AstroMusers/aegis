import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


class Exoplanet:
    def __init__(self, name, semi_major_axis, radius, mass, density, magnetic_field, star_mass, star_mass_loss,
                 distance, freq=0, intensity=0):
        self.name = name
        self.semi_major_axis = semi_major_axis
        self.radius = radius
        self.mass = mass
        self.density = density
        self.magnetic_field = magnetic_field
        self.star_mass = star_mass
        self.star_mass_loss = star_mass_loss
        self.distance = distance
        self.freq = freq
        self.intensity = intensity
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


def convective_radius(M, p, r):
    """
    Calculates the radius of the convective core of the exoplanet with the scaling provided in Curtis & Ness 1986.
    If the calculation exceeds planetary radius, planetary radius is returned. Jovain convective radius of 0.830 R_j
    is taken from Sharan 2022.
    :param M: Mass of the planet in Mj
    :param p: Density of the planet in g/cm3
    :param r: Radius of the planet in Rj

    :return: Radius of the convective core
    """
    # M = exo.mass
    # p = exo.density
    # r = exo.radius
    if p < 1.6 and M > 0.4:
        R = 0.049 * M ** 0.44  # Curtis & Ness (1986)
    else:
        R = 0.830 * M ** 0.75 * r ** (-0.96)  # Grießmeier (2004)
    return min(R, r) / 0.830


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


def freq_condition(nu, n):
    n /= 10**6
    nu_p = 8.98 * np.sqrt(n) * 10**(-3)
    return nu > nu_p


def keplerian(M, r):
    """

    :param M: Host star's mass in M_sun
    :param r: Semi-major axis in AU
    :return: Keplerian velocity in km/s
    """
    G = 8.87129 * 10 ** 2  # Gravitational constant in AU (M_sun)-1 (km/s)2
    return np.sqrt(G * M / r)


def v_eff(v_wind, v_kepler):
    """

    :param v_wind: Stellar wind speed
    :param v_kepler: Keplerian speed of the planet
    :return: effective speed of the planet in stellar wind, in units of the provided speeds.
    """
    return np.sqrt(v_wind**2 + v_kepler**2)


def star_surface_b(age):
    """

    :param age: Age of the star in Gyr
    :return: Surface magentic field of the star in Gauss
    """
    return 0.7 * age**0.655


def star_period(age):
    """

    :param age: Age of the STar in Gyr
    :return: Rotational period of the star in days.
    """
    tau = 2.56 * 10**(-2)
    return 0.67 * (1 + age / tau)**0.7


def B_r(B0, R0, r):
    """

    :param B0: Surface magentic field of the Star
    :param R0: Radius of the Star in solar radii
    :param r: distance from the star in AU
    :return: Radial component of the interplanetary magnetic field of the star.
    """
    R0 /= 215  # Conversion from solar radii to AU
    return B0 * (R0 / r)**2


def B_phi(B_r, P, r, v_eff):
    """

    :param B_r: Radial component of IMF
    :param P: Stellar rotation period in days
    :param r: Distance from star (sma) in AU
    :param v_eff: Effective speed of the planet in stellar wind.
    :return: Phi component of the IMF.
    """
    P *= 8.64 * 10**4  # Converted from days to seconds
    r *= 1.496 * 10**8  # Converted from AU to km
    Omega = 2 * np.pi / P
    return B_r * Omega * r / v_eff


def imf(B_r, B_phi):
    return np.sqrt(B_r**2 + B_phi**2)


def imf_complete(M, r, v_wind, age, R):
    v_k = keplerian(M, r)
    veff = v_eff(v_wind, v_k)
    B0 = star_surface_b(age)
    P = star_period(age)
    Br = B_r(B0, R, r)
    Bphi = B_phi(Br, P, r, veff)
    return imf(Br, Bphi)


def imf_perp(B_r, B_phi, v_k, v):
    """

    :param B_r: Radial component of IMF
    :param B_phi: Phi component of IMF
    :param v_k: Keplerian velocity of the planet
    :param v: Stellar wind velocity at planet distance
    :return:
    """
    return np.sqrt(B_r**2 + B_phi**2) * abs(np.sin(np.arctan(B_phi/B_r) - np.arctan(v_k/v)))


def imf_perp_complete(M, r, R, age, v_wind):
    """

    :param M: Mass of the Star
    :param r: Distance from the star
    :param R: Radius of the star
    :param age: Age of the star
    :param v_wind: Wind velocity of the star at planet distance
    :return: Perpendicular component of the Inteplanetary Magnetic Field at the distance of the planet
    """
    v_k = keplerian(M, r)
    veff = v_eff(v_wind, v_k)
    B0 = star_surface_b(age)
    P = star_period(age)
    Br = B_r(B0, R, r)
    Bphi = B_phi(Br, P, r, veff)
    B_perp = imf_perp(Br, Bphi, v_k, v_wind)
    return B_perp


# def number_density(Mdot, v_eff, r):
#     """
#
#     :param Mdot: Mass loss rate of the star in 10^-15 Msun/yr
#     :param v_eff: stellar wind speed in km/s
#     :param r: distance in AU
#     :return: number density in m^-3
#     """
#     # Mdot *= 2 * 10**16
#     # v_eff *= 10**3
#     # r *= 1.49 * 10**11
#     # m_p = 1.76 * 10**(-27)  # Mass of proton
#     # return Mdot / (4 * np.pi * m_p * v_eff * r**2)

def number_density(age):
    n0 = 1.04 * 10**11
    tau = 2.56 * 10**(-2)
    return n0 * (1 + age/tau)**(-1.86)



def Rm(B, R, n, T, v_eff, B_star):
    """

    :param B: Planetary surface magnetic flux ddensity in G
    :param R: Radius of the planet in
    :param n: Number density of the stellar wind in m^-3
    :param T: Temperature of the stellar wind in K
    :param v_eff: Speed of the stellar wind in km/s
    :param B_star: Magnetic field flux density of the star in G
    :return: Magnetopause distance of the planet.
    """
    k_b = 1.38 * 10**(-23)  # Boltzmann constant in J/K
    m_p = 1.67 * 10**(-27)  # Mass of proton
    R_magnet = 2.44**(1/3) * ((B**2 / (80 * np.pi * (2 * n * k_b * T + m_p * n * 10**6 * v_eff**2 + B_star**2/(80 * np.pi)))) ** (1/6)) * R
    if R_magnet > R:
        return R_magnet
    else:
        return R


def P_input_mag(imf_perp, v_eff, R_m, n):
    """
    I'm actually not sure about the units.
    :param imf_perp: Perpendicular component of IMF flux density in G
    :param v_eff: EFfective speed near the exoplanet in m/s
    :param R_m: Radius of planet Magnetosphere in m
    :param n: Number density of the wind in m^-3
    :return: Input Radio Power in Watts
    """
    m_p = 1.67 * 10**(-27)  # Mass of proton
    mu_0 = 4 * np.pi * 10**(-7)
    # return imf_perp**2 / (2*mu_0) * v_eff * R_m**2
    return imf_perp**2 / (80*np.pi) * v_eff * R_m**2


def P_input_kin(imf_perp, v_eff, R_m, n):
    """
    I'm actually not sure about the units.
    :param imf_perp: Perpendicular component of IMF flux density in G
    :param v_eff: EFfective speed near the exoplanet in m/s
    :param R_m: Radius of planet Magnetosphere in m
    :param n: Number density of the wind in m^-3
    :return: Input Radio Power in Watts
    """
    m_p = 1.67 * 10**(-27)  # Mass of proton
    mu_0 = 4 * np.pi * 10**(-7)
    # return imf_perp**2 / (2*mu_0) * v_eff * R_m**2
    return m_p * n * v_eff**3 * np.pi * R_m**2


def radio_power(P_input_mag, P_input_kin, nu, d):
    """

    :param P_input: Input power in Watts
    :param nu: Bandwidth of observation (Assumed to be maximum emission frequency)
    :param d: Distance from Earth to the Source.
    :return: Expected observed radio flux density in Jy.
    """
    epsilon_mag = 3.32 * 10**(-4)  # for magnetic power
    epsilon_kin = 7.86 * 10**(-7)  # Only if you combine incident kinetic and magnetic power
    return (epsilon_mag * P_input_mag + epsilon_kin * P_input_kin) / (1.6 * nu * d**2) * 10**26


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

# Plotting


mpl.use('Qt5Agg')

mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [10, 5]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)


def retro_noir(ax):
    ax.grid(alpha=0.2)
    ax.tick_params(direction="in", length=7, right=True, top=True, width=1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
