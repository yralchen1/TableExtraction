"""
ASTools — Atomic Spectroscopy Tools

Python translation of the VBA module ASTools.bas.
Contains functions for:
  - Physical constants (CODATA vintages)
  - Refractive index of air (multiple formulations)
  - Air↔vacuum wavelength conversions
  - Water vapor pressure
  - Transition rates, line strengths, oscillator strengths
  - Angular momentum coupling (6j-symbols, etc.)
  - NIST accuracy codes
  - Formatting / rounding utilities
  - Rowland wavelength correction
"""

import math
import numpy as np
from scipy.special import factorial2, gamma as scipy_gamma


# ════════════════════════════════════════════════════════════════
#  SECTION 1: Physical Constants (CODATA vintages)
# ════════════════════════════════════════════════════════════════

def C_cm_sec():
    """Speed of light in cm/s."""
    return 29979245800.0

def C_m_sec():
    """Speed of light in m/s."""
    return C_cm_sec() * 0.01

def Ry(year=2022):
    """Rydberg constant in cm⁻¹."""
    vals = {
        2022: 109737.31568157, 2018: 109737.3156816,
        2014: 109737.31568508, 2010: 109737.31568539,
        2006: 109737.31568527, 2002: 109737.31568525,
        1998: 109737.31568549, 1986: 109737.31534,
        1973: 109737.3177, 1969: 109737.312,
    }
    return vals.get(year, vals[2022])

def evcm(year=2022):
    """eV to cm⁻¹ conversion factor."""
    vals = {
        2022: 8065.54393734921, 2018: 8065.54393734921,
        2014: 8065.544005, 2010: 8065.54429,
        2006: 8065.54465, 2002: 8065.54445,
        1998: 8065.54477, 1986: 8065.541,
        1973: 8065.479, 1969: 8065.469,
    }
    return vals.get(year, vals[2022])

def cmev(year=2022):
    """cm⁻¹ to eV conversion factor."""
    return 1.0 / evcm(year)

def evcm_Moore_1970():
    return 1.0 / 0.000123981

def evcm_1952():
    return 1.0 / 0.00012395

def evcm_1932():
    return 1.0 / 0.00012336

def u_kg(year=2022):
    """Atomic mass unit in kg."""
    vals = {2022: 1.66053906892e-27, 2018: 1.6605390666e-27}
    return vals.get(year, vals[2022])

def u_eV(year=2022):
    """Atomic mass unit in eV."""
    vals = {
        2022: 931494103.72, 2018: 931494102.42,
        2014: 931494095.37, 2010: 931494060.92,
        2006: 931494027.62, 2002: 931494042.91,
        1998: 931494012.15, 1986: 931494320.1,
        1973: 931501600.0, 1969: 931481400.0,
    }
    return vals.get(year, vals[2022])

def m_e_u(year=2022):
    """Electron mass in atomic mass units."""
    vals = {
        2022: 5.485799090441e-04, 2018: 0.000548579909065,
        2014: 0.00054857990907, 2002: 0.00054857990945,
        1998: 0.000548579911, 1973: 0.00054858026,
    }
    return vals.get(year, vals[2022])

def alpha_fine_struc(year=2022):
    """Fine-structure constant."""
    vals = {2022: 0.0072973525643, 2018: 0.0072973525693}
    return vals.get(year, vals[2022])

def mu_B(year=2022):
    """Bohr magneton in J/T."""
    vals = {2022: 9.2740100657e-24, 2018: 9.2740100783e-24}
    return vals.get(year, vals[2022])

def Bohr_rad_m(year=2018):
    """Bohr radius in meters."""
    vals = {2018: 5.29177210903e-11, 2006: 5.2917720859e-11}
    return vals.get(year, vals[2018])

def K2eV():
    """Convert Kelvins to eV. CODATA 2018."""
    return 8.61733326214518e-05

def freq_au():
    """Atomic unit of frequency (s⁻¹)."""
    return 4.0 * math.pi * Ry() * C_cm_sec()

def torr2Pa(p):
    """Convert Torr to Pascal."""
    return p / 0.00750061683

def Pa2torr(Pa):
    """Convert Pascal to Torr."""
    return Pa * 0.00750061683

def C2K(t):
    """Celsius to Kelvin."""
    return t + 273.15


# ════════════════════════════════════════════════════════════════
#  SECTION 2: Refractive Index of Air
# ════════════════════════════════════════════════════════════════

def nair1e(sigma):
    """
    (n-1) of standard air by Edlén 1966.
    B. Edlén, Metrologia 2, 71–80 (1966).
    Dry air at 15°C, 760 mmHg, 0.033% CO₂.
    sigma: wavenumber in cm⁻¹.
    """
    sigma2 = sigma * sigma
    return (0.0000834213
            + 2406030.0 / (13000000000.0 - sigma2)
            + 15997.0 / (3890000000.0 - sigma2))

def nair1e53(sigma):
    """
    (n-1) of standard air by Edlén 1953.
    B. Edlén, JOSA 43, 339–344 (1953).
    sigma: wavenumber in cm⁻¹.
    """
    s = sigma * 0.0001
    s2 = s * s
    return (6432.8 + 2949810.0 / (146.0 - s2) + 25540.0 / (41.0 - s2)) * 1e-8

def nair1(sigma):
    """
    (n-1) by Peck & Reeder 1972.
    E.R. Peck and K. Reeder, JOSA 62, 958–962 (1972).
    sigma: wavenumber in cm⁻¹.
    """
    sigma2 = sigma * sigma
    return (0.0000806051
            + 2480990.0 / (13227400000.0 - sigma2)
            + 17455.7 / (3932957000.0 - sigma2))

def nair1MP(sigma):
    """
    (n-1) by Meggers & Peters 1919.
    W.F. Meggers and C.G. Peters, Astrophys. J. 50, 56–71 (1919).
    sigma: wavenumber in cm⁻¹.
    """
    lam = Lair(1e8 / sigma)
    lam2 = lam * lam
    return (2726.43 + 1228800000.0 / lam2 + 3.555e15 / (lam2 * lam2)) * 1e-7

def nair1B51(sigma):
    """
    (n-1) by Barrell 1951.
    H. Barrell, JOSA 41, 295–298 (1951).
    sigma: wavenumber in cm⁻¹.
    """
    lam = 10000.0 / sigma  # vacuum wavelength in μm
    lam2 = lam * lam
    return (272.729 + 1.4823 / lam2 + 0.02041 / (lam2 * lam2)) * 1e-6

def nair1B51Eq7(lambda_s):
    """
    (n-1) by Barrell 1951, Eq.(7).
    lambda_s: wavelength in standard air, in Å.
    """
    lam = 0.0001 * lambda_s  # convert to μm
    lam2 = lam * lam
    return (272.729 + 1.4814 / lam2 + 0.02039 / (lam2 * lam2)) * 1e-6

def nair1B51PT(sigma, p, t):
    """
    (n-1) by Barrell 1951 at pressure p [mmHg] and temperature t [°C].
    """
    n1s = nair1B51(sigma)
    alpha = 0.003661
    beta15 = 0.0000008135
    beta_t = 1e-6 * (1.049 - 0.0157 * t)
    return n1s * p * (1 + beta_t * p) * (1 + 15 * alpha) / (760 * (1 + 760 * beta15) * (1 + alpha * t))


# ════════════════════════════════════════════════════════════════
#  SECTION 2b: Water Vapor Pressure
# ════════════════════════════════════════════════════════════════

def wPs(t):
    """
    Saturated water vapor pressure in Pa at temperature t (°C).
    Stone-Zimmerman / modified Ciddor equations.
    http://emtoolbox.nist.gov/Wavelength/Documentation.asp
    """
    Tk = t + 273.15
    K1, K2, K3, K4 = 1167.05214528, -724213.167032, -17.0738469401, 12020.8247025
    K5, K6, K7, K8 = -3232555.03223, 14.9151086135, -4823.26573616, 405113.405421
    K9, K10 = -0.238555575678, 650.175348448

    if t >= 0:
        omega = Tk + K9 / (Tk - K10)
        omega2 = omega * omega
        a = omega2 + K1 * omega + K2
        b = K3 * omega2 + K4 * omega + K5
        c = K6 * omega2 + K7 * omega + K8
        X = -b + math.sqrt(b * b - 4 * a * c)
        p_sv = (2 * c / X)
        p_sv = p_sv * p_sv
        p_sv = 1e6 * p_sv * p_sv
    else:
        A1 = -13.928169
        A2 = 34.7078238
        theta = Tk / 273.16
        theta_sqrt = math.sqrt(theta)
        theta_pow_quart = math.sqrt(theta_sqrt)
        y = A1 * (1 - 1 / (theta * theta_sqrt)) + A2 * (1 - 1 / (theta * theta_pow_quart))
        p_sv = 611.657 * math.exp(y)
    return p_sv

def wPs_Buck(t):
    """
    Saturated water vapor pressure in Pa (Buck 1981/1996).
    """
    if t < 0:
        return 611.15 * math.exp((23.036 - t / 333.7) * (t / (279.82 + t)))
    else:
        return 611.21 * math.exp((18.678 - t / 234.5) * (t / (257.14 + t)))

def wPs_SW87(t):
    """
    Saturated water vapor pressure in Pa.
    Saul & Wagner 1987 / Wagner & Pruß 2002.
    """
    Tc = 647.096
    pc = 22064000.0
    Tk = t + 273.15
    tau = Tc / Tk
    theta = 1 - Tk / Tc
    theta_sqrt = math.sqrt(theta)
    theta3 = theta ** 3
    theta4 = theta3 * theta
    theta7 = theta4 * theta3
    A1, A2 = -7.85951783, 1.84408259
    A3, A4 = -11.7866497, 22.6807411
    A5, A6 = -15.9618719, 1.80122502
    return pc * math.exp(tau * (A1 * theta + A2 * theta * theta_sqrt
                                + A3 * theta3 + A4 * theta3 * theta_sqrt
                                + A5 * theta4 + A6 * theta7 * theta_sqrt))


# ════════════════════════════════════════════════════════════════
#  SECTION 2c: Ciddor / Stone-Zimmerman refractive index
# ════════════════════════════════════════════════════════════════

def nair1Ciddor(sigma, p, t, w, xCO2):
    """
    (n-1) of humid air — phase refractive index.
    Stone-Zimmerman / modified Ciddor.
    sigma: cm⁻¹, p: Torr, t: °C, w: % rel. humidity, xCO2: % CO₂.
    """
    return nair1Ciddor_gp(sigma, p, t, w, xCO2, group=0)

def nair1Ciddor_gp(sigma, p, t, w, xCO2, group=0):
    """
    (n-1) phase or group refractive index of humid air.
    group=0 → phase, group=1 → group.
    """
    kB = 1.3806505e-23
    Tk = t + 273.15
    Pa = torr2Pa(p)
    Pw = w / 100.0 * wPs(t)

    # Dry-air dispersion (Ciddor Eq. B1-B2)
    s2 = (sigma / 10000.0) ** 2
    if group == 0:
        n_as = 1e-8 * (5792105.0 / (238.0185 - s2) + 167917.0 / (57.362 - s2))
    else:
        n_as = 1e-8 * (5792105.0 * (238.0185 + s2) / (238.0185 - s2)**2
                       + 167917.0 * (57.362 + s2) / (57.362 - s2)**2)

    # Water vapor dispersion (Ciddor Eq. B3)
    cf = 1.022
    w0, w1, w2, w3 = 295.235, 2.6422, -0.032380, 0.004028
    if group == 0:
        r_vs = cf * 1e-8 * (w0 + w1 * s2 + w2 * s2**2 + w3 * s2**3)
    else:
        r_vs = cf * 1e-8 * (w0 + w1 * 3 * s2 + w2 * 5 * s2**2 + w3 * 7 * s2**3)

    # Modified ideal-gas constants
    R = 8.314472
    M_a = 1e-3 * (28.9635 + 12.011e-6 * (xCO2 * 10000 - 400))
    m_v = 0.018015

    # Compressibility
    a0 = 1.58123e-6
    a1 = -2.9331e-8
    a2 = 1.1043e-10
    b0 = 5.707e-6
    b1 = -2.051e-8
    c0 = 1.9898e-4
    c1 = -2.376e-6
    d = 1.83e-11
    e_val = -0.765e-8

    Z_a = 1 - (Pa / Tk) * (a0 + a1 * t + a2 * t**2 + (b0 + b1 * t) * (xCO2 / 100.0) + (c0 + c1 * t) * (Pw / Pa if Pa != 0 else 0)**2) + (Pa / Tk)**2 * (d + e_val * (Pw / Pa if Pa != 0 else 0)**2)

    # Reference densities
    P_R1 = 101325.0
    T_R1 = 288.15
    Z_a_ref = 1 - (P_R1 / T_R1) * (a0 + a1 * 15 + a2 * 225)
    rho_ax5 = P_R1 * M_a / (Z_a_ref * R * T_R1)

    T_vap = 293.15
    svp20 = wPs(20)
    Z_v_ref = 1 - (svp20 / T_vap) * (a0 + a1 * 20 + a2 * 400 + (c0 + c1 * 20))
    rho_v5 = svp20 * m_v / (Z_v_ref * R * T_vap)

    # Actual densities
    x_v = Pw / Pa if Pa != 0 else 0
    rho_a = (1 - x_v) * Pa * M_a / (Z_a * R * Tk)
    rho_v = x_v * Pa * m_v / (Z_a * R * Tk)

    r_ax5 = n_as * (1 + n_as) * (1 + 0.534e-6 * (xCO2 * 10000 - 450)) if group == 0 else n_as
    return (rho_a / rho_ax5) * r_ax5 + (rho_v / rho_v5) * r_vs


# ════════════════════════════════════════════════════════════════
#  SECTION 2d: Birch 1994
# ════════════════════════════════════════════════════════════════

def nair1_Birch94(sigma):
    """
    (n-1) by Birch & Downs 1994.
    K.P. Birch and M.J. Downs, Metrologia 31, 315 (1994).
    sigma: cm⁻¹. Dry air, 15°C, 0.033% CO₂.
    """
    s2 = (sigma / 10000.0) ** 2
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    K1 = 1 + nair1Ciddor(sigma, 760, 15, 0, 0.045)
    K2 = 1 + nair1Ciddor(sigma, 760, 15, 0, 0.033)
    k = K2 / K1
    return n * k - 1

def nair1_Birch94_group(sigma):
    """(n_group - 1) by Birch & Downs 1994."""
    s2 = (sigma / 10000.0) ** 2
    n = (1 + 0.0000834254
         + 0.02406147 * (130 + s2) / (130 - s2)**2
         + 0.00015998 * (38.9 + s2) / (38.9 - s2)**2)
    K1 = 1 + nair1Ciddor_gp(sigma, 760, 15, 0, 0.045, 1)
    K2 = 1 + nair1Ciddor_gp(sigma, 760, 15, 0, 0.033, 1)
    k = K2 / K1
    return n * k - 1

def nair1_Birch94_ext(sigma, p, t, w, xCO2):
    """(n-1) Birch 1994 extended to arbitrary conditions."""
    return nair1_Birch94_ext_gp(sigma, p, t, w, xCO2, 0)

def nair1_Birch94_ext_gp(sigma, p, t, w, xCO2, group=0):
    """
    (n-1) Birch 1994 extended.
    Birch & Downs, Metrologia 31, 315 (1994); Metrologia 30, 155 (1993).
    """
    s2 = (sigma / 10000.0) ** 2
    Pa_val = torr2Pa(p)
    w_Pa = wPs(t) * w / 100.0

    if group == 0:
        n1s = 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
        n1tp = Pa_val * n1s / 96095.43 * (1 + 1e-8 * (0.601 - 0.00972 * t) * Pa_val) / (1 + 0.003661 * t)
    else:
        n1s = (0.0000834254
               + 0.02406147 * (130 + s2) / (130 - s2)**2
               + 0.00015998 * (38.9 + s2) / (38.9 - s2)**2)
        n1tp = Pa_val * n1s / 96095.43 * (1 + 1e-8 * (0.601 - 0.00972 * t) * Pa_val) / (1 + 0.003661 * t)

    n1tpw = n1tp - w_Pa * (3.7345 - 0.0401 * s2) * 1e-10
    return n1tpw


# ════════════════════════════════════════════════════════════════
#  SECTION 3: Wavelength Conversions
# ════════════════════════════════════════════════════════════════

def _iterative_Lvac(La, nair1_func, *args):
    """Generic iterative air→vacuum wavelength conversion."""
    if La <= 1500.0:
        return La
    Lv0 = La * 1.000281383
    tol = 1e-18
    for _ in range(1000):
        Lv1 = La * (1 + nair1_func(1e8 / Lv0, *args))
        if abs(Lv1 - Lv0) <= tol:
            return Lv1
        Lv0 = Lv1
    return Lv1

def LvacE(La):
    """Vacuum wavelength from air wavelength (Å) — Edlén 1966."""
    return _iterative_Lvac(La, nair1e)

def LairE(Lvac):
    """Air wavelength from vacuum wavelength (Å) — Edlén 1966."""
    return Lvac / (1 + nair1e(1e8 / Lvac))

def Lvac(La):
    """Vacuum wavelength from air wavelength (Å) — Peck & Reeder 1972."""
    return _iterative_Lvac(La, nair1)

def Lair(Lvac_val):
    """Air wavelength from vacuum wavelength (Å) — Peck & Reeder 1972."""
    return Lvac_val / (1 + nair1(1e8 / Lvac_val))

def LvacE53(La):
    """Vacuum wavelength — Edlén 1953."""
    return _iterative_Lvac(La, nair1e53)

def LairE53(Lvac_val):
    """Air wavelength — Edlén 1953."""
    return Lvac_val / (1 + nair1e53(1e8 / Lvac_val))

def LvacMP(La):
    """Vacuum wavelength — Meggers & Peters 1919."""
    return _iterative_Lvac(La, nair1MP)

def LairMP(Lvac_val):
    """Air wavelength — Meggers & Peters 1919."""
    return Lvac_val / (1 + nair1MP(1e8 / Lvac_val))

def LvacB51(La):
    """Vacuum wavelength — Barrell 1951."""
    return _iterative_Lvac(La, nair1B51)

def LairB51(Lvac_val):
    """Air wavelength — Barrell 1951."""
    return Lvac_val / (1 + nair1B51(1e8 / Lvac_val))

def LvacB51PT(La, p, t):
    """Vacuum wavelength — Barrell 1951 at p [mmHg], t [°C]."""
    return _iterative_Lvac(La, nair1B51PT, p, t)

def LairB51PT(Lvac_val, p, t):
    """Air wavelength — Barrell 1951 at p [mmHg], t [°C]."""
    return Lvac_val / (1 + nair1B51PT(1e8 / Lvac_val, p, t))

def LvacCiddor(La):
    """Vacuum wavelength — Ciddor (standard conditions)."""
    return _iterative_Lvac(La, nair1Ciddor, 760, 15, 0, 0.033)

def LairCiddor(Lvac_val):
    """Air wavelength — Ciddor (standard conditions)."""
    return Lvac_val / (1 + nair1Ciddor(1e8 / Lvac_val, 760, 15, 0, 0.033))

def nair1MPext(sigma, p, t):
    """(n-1) Meggers & Peters at pressure p [mmHg], temperature t [°C]."""
    lam = Lair(1e8 / sigma)
    lam2 = lam * lam
    n1s = nair1MP(sigma)
    alpha = 0.00367
    lam3 = lam2 * lam
    return n1s * p / 760 * (1 + 15 * (alpha + 3e6 / lam3)) / (1 + t * (alpha + 3e6 / lam3))

def LvacMPext(La, p, t):
    """Vacuum wavelength — Meggers & Peters at p, t."""
    return _iterative_Lvac(La, nair1MPext, p, t)

def LairMPext(Lvac_val, p, t):
    """Air wavelength — Meggers & Peters at p, t."""
    return Lvac_val / (1 + nair1MPext(1e8 / Lvac_val, p, t))


# ════════════════════════════════════════════════════════════════
#  SECTION 4: Transition Rates / Line Strengths / Oscillator Strengths
# ════════════════════════════════════════════════════════════════

# Coefficients for E/M multipole transitions (2018 CODATA)
_TRANS_COEFFS = {
    "E1": 2.02612688864176e18,
    "M1": 26973500312730.8,
    "E2": 1.1199500334572e18,
    "M2": 14909714069266.2,
    "E3": 3.14441661858975e17,
    "M3": 4186111102930.18,
}

_F_COEFFS = {
    "E1": 303.755683527197,
    "M1": 4.04385039779375e-3,
    "E2": 167.902212756848,
    "M2": 2.23525506407994e-3,
    "E3": 47.1408984614135,
    "M3": 6.27579174097901e-4,
}

def _multipole_order(type_str):
    """Return the multipole order p from the type string."""
    if not type_str or type_str == "E1" or type_str == "M1":
        return 1
    elif type_str in ("E2", "M2"):
        return 2
    elif type_str in ("E3", "M3"):
        return 3
    return 0

def _wl_power(type_str):
    """Return the wavelength exponent (2p+1) for A-coefficient."""
    p = _multipole_order(type_str)
    return 2 * p + 1

def Np_factor(p):
    """
    Factor for multipole transition probabilities.
    Sobel'man (1972); Grant (2020) Eq. 13.1.14.
    p: multipole order (integer).
    """
    FD = factorial2(2 * p + 1, exact=True)
    return 2 * (2 * p + 1) * (p + 1) * (2 * math.pi) ** (2 * p + 1) / (p * FD * FD)

def AfromS(S, vac_wl, g_upp, type_str="E1"):
    """
    Transition probability A (s⁻¹) from line strength S (a.u.),
    vacuum wavelength vac_wl (Å), and statistical weight g_upp.
    type_str: 'E1','M1','E2','M2','E3','M3'.
    """
    key = type_str if type_str else "E1"
    coeff = _TRANS_COEFFS.get(key, 0)
    pw = _wl_power(key)
    return S * coeff / vac_wl ** pw / g_upp if coeff else 0.0

def gAfromS(S, vac_wl, type_str="E1"):
    """Weighted transition probability gA from line strength S."""
    key = type_str if type_str else "E1"
    coeff = _TRANS_COEFFS.get(key, 0)
    pw = _wl_power(key)
    return S * coeff / vac_wl ** pw if coeff else 0.0

def SfromA(A, vac_wl, g_upp, type_str="E1"):
    """Line strength S (a.u.) from A (s⁻¹)."""
    key = type_str if type_str else "E1"
    coeff = _TRANS_COEFFS.get(key, 0)
    pw = _wl_power(key)
    return A * g_upp * vac_wl ** pw / coeff if coeff else 0.0

def FfromA(A, vac_wl, g_low, g_upp):
    """Oscillator strength f from A (s⁻¹) — E1 only."""
    return 1.49919e-16 * A * g_upp * vac_wl ** 2 / g_low

def AfromF(f, vac_wl, g_low, g_upp):
    """A (s⁻¹) from oscillator strength f — E1 only."""
    return 0.66702e16 * f * g_low / (g_upp * vac_wl ** 2)

def FfromS(S, vac_wl, g_low, type_str="E1"):
    """Oscillator strength f from line strength S."""
    key = type_str if type_str else "E1"
    coeff = _F_COEFFS.get(key, 0)
    pw = 2 * _multipole_order(key) - 1
    return S * coeff / vac_wl ** pw / g_low if coeff else 0.0

def SfromF(f, vac_wl, g_low, type_str="E1"):
    """Line strength S from oscillator strength f."""
    key = type_str if type_str else "E1"
    coeff = _F_COEFFS.get(key, 0)
    pw = 2 * _multipole_order(key) - 1
    return (f * g_low / coeff) * vac_wl ** pw if coeff else 0.0

def SfromgF(gf, vac_wl, type_str="E1"):
    """Line strength S from weighted oscillator strength gf."""
    key = type_str if type_str else "E1"
    coeff = _F_COEFFS.get(key, 0)
    pw = 2 * _multipole_order(key) - 1
    return (gf / coeff) * vac_wl ** pw if coeff else 0.0

def AfromgF(gf, vac_wl, g_upp):
    """A (s⁻¹) from weighted oscillator strength gf — E1 only."""
    return 0.66702e16 * gf / (g_upp * vac_wl ** 2)


# ════════════════════════════════════════════════════════════════
#  SECTION 5: Angular Momentum
# ════════════════════════════════════════════════════════════════

_L_LETTERS = "SPDFGHIKLMNOQRTUVWXYZ"

def LQN(Lstr):
    """Orbital angular momentum quantum number from spectroscopic letter."""
    idx = _L_LETTERS.find(str(Lstr).upper())
    if idx < 0:
        raise ValueError(f"Unknown L designation: {Lstr}")
    return idx

def LfromLQN(L):
    """Spectroscopic letter from orbital quantum number L."""
    return _L_LETTERS[L]

def FCT(n):
    """
    Factorial for 6j-symbol calculations.
    Accepts integer n. Returns n! (with FSCALE=1).
    """
    return float(math.factorial(max(n, 0)))

def DELSQ(IA, IB, IC, ID):
    """Triangle coefficient for 6j-symbol. Arguments are 2*J values (integers)."""
    return FCT(ID - IA) * FCT(ID - IB) * FCT(ID - IC) / FCT(ID + 1)

def S6J(FJ1, FJ2, FJ3, FL1, FL2, FL3):
    """
    Wigner 6j-symbol { FJ1 FJ2 FJ3 }
                      { FL1 FL2 FL3 }.
    All arguments are angular momenta (float, can be half-integer).
    """
    s1 = FJ1 + FJ2 + FJ3
    s2 = FJ1 + FL2 + FL3
    s3 = FL1 + FJ2 + FL3

    IS1 = int(s1)
    IS2 = int(s2)
    IS3 = int(s3)

    # Check triangular conditions (sums must be integers)
    if (abs(s1 - IS1) > 0.01 or abs(s2 - IS2) > 0.01 or abs(s3 - IS3) > 0.01):
        return 0.0

    IS4 = int(FL1 + FL2 + FJ3)
    IS5 = int(FJ1 + FJ2 + FL1 + FL2)
    IS6 = int(FJ2 + FJ3 + FL2 + FL3)
    IS7 = int(FJ3 + FJ1 + FL3 + FL1)

    KN = max(IS1, IS2, IS3, IS4)
    KX = min(IS5, IS6, IS7)

    if KX < KN:
        return 0.0

    IT1 = int(2.000001 * FJ1)
    IT2 = int(2.000001 * FJ2)
    IT3 = int(2.000001 * FJ3)
    IT4 = int(2.000001 * FL1)
    IT5 = int(2.000001 * FL2)
    IT6 = int(2.000001 * FL3)

    COEF = (DELSQ(IT1, IT2, IT3, IS1) * DELSQ(IT1, IT5, IT6, IS2)
            * DELSQ(IT4, IT2, IT6, IS3) * DELSQ(IT4, IT5, IT3, IS4))
    COEF = math.sqrt(COEF)

    mySGN = -1 + 2 * (KN % 2)
    total = 0.0
    for k in range(KN, KX + 1):
        mySGN = -mySGN
        total += mySGN * FCT(k + 1) / (
            FCT(k - IS1) * FCT(k - IS2) * FCT(k - IS3) * FCT(k - IS4)
            * FCT(IS5 - k) * FCT(IS6 - k) * FCT(IS7 - k))
    return COEF * total

def dEhfs(a, b, I, J, F):
    """
    Hyperfine structure energy shift.
    a: magnetic dipole constant (cm⁻¹ or MHz).
    b: electric quadrupole constant.
    I: nuclear spin (float). J: total electronic angular momentum (float).
    F: total angular momentum (float).
    """
    c = F * (F + 1) - J * (J + 1) - I * (I + 1)
    if J == 0.5 or I == 0.5:
        dEquad = 0.0
    else:
        dEquad = b * (3 * c * (c + 1) / 4 - I * (I + 1) * J * (J + 1)) / (
            2 * I * (2 * I - 1) * J * (2 * J - 1))
    return a * c / 2 + dEquad

def Gtot(SL):
    """Total statistical weight of a term (e.g., '4P')."""
    s = int(SL[0])
    L = LQN(SL[1])
    Lnum = 2 * L + 1
    gmin = abs(Lnum - s)
    gmax = Lnum + s - 2
    num_J = (gmax - gmin) // 2
    G0 = gmin + 1
    g1 = G0
    for _ in range(1, num_J + 1):
        g1 += 2
        G0 += g1
    return G0

def numJ(SL):
    """Number of J-levels in a term (e.g., '4P')."""
    s = int(SL[0])
    L = LQN(SL[1])
    Lnum = 2 * L + 1
    gmin = abs(Lnum - s)
    gmax = Lnum + s - 2
    return (gmax - gmin) // 2 + 1

def termG_LS(Term):
    """Total statistical weight of an LS term."""
    s = int(Term[0])
    L = LQN(Term[1])
    Snum = (s - 1) / 2.0
    Jmin = abs(L - Snum)
    Jmax = L + Snum
    num_lev = int(Jmax - Jmin + 1)
    gnum = 0
    for i in range(num_lev):
        j = Jmin + i
        gnum += int(2 * j + 1)
    return gnum

def g_mult(Term):
    """Statistical weight of a multiplet (sum of 2J+1 over all J)."""
    s = int(Term[0])
    Sqn = (s - 1) / 2.0
    L = LQN(Term[1])
    Gm = 0
    j = abs(L - Sqn)
    while j <= L + Sqn:
        Gm += int(2 * j + 1)
        j += 1
    return Gm

def DLine(Term1, J1, Term2, J2):
    """
    Reduced dipole matrix element ratio D(J1,J2) for LS coupling.
    Term1, Term2: e.g. '2P'. J1, J2: float (half-integer OK).
    """
    s1 = Term1[0]
    s2 = Term2[0]
    if s1 != s2:
        return 0.0
    Sqn = (int(s1) - 1) / 2.0
    Lqn1 = LQN(Term1[1])
    Lqn2 = LQN(Term2[1])
    d = S6J(Lqn1, Sqn, J1, J2, 1, Lqn2)
    d *= math.sqrt((2 * J1 + 1) * (2 * J2 + 1))
    k = int(Lqn1 + Sqn + J2 + 1) % 2
    return -d if k != 0 else d

def DLine2(Term1, J1, Term2, J2):
    """Squared reduced dipole matrix element D²(J1,J2)."""
    s1 = Term1[0]
    s2 = Term2[0]
    if s1 != s2:
        return 0.0
    Sqn = (int(s1) - 1) / 2.0
    Lqn1 = LQN(Term1[1])
    Lqn2 = LQN(Term2[1])
    d = S6J(Lqn1, Sqn, J1, J2, 1, Lqn2)
    return d * d * (2 * J1 + 1) * (2 * J2 + 1)

def LineStrFromMult(Smul, Term1, J1, Term2, J2):
    """
    Line strength for a single J→J' transition from the multiplet strength Smul.
    """
    s1 = Term1[0]
    s2 = Term2[0]
    if s1 != s2:
        return 0.0
    Sqn = (int(s1) - 1) / 2.0
    Lqn1 = LQN(Term1[1])
    Lqn2 = LQN(Term2[1])
    # Sum D² over all J, J' to get total
    Stot = 0.0
    j1 = abs(Lqn1 - Sqn)
    while j1 <= Lqn1 + Sqn:
        j2 = abs(Lqn2 - Sqn)
        while j2 <= Lqn2 + Sqn:
            D2 = DLine2(Term1, j1, Term2, j2)
            Stot += D2
            j2 += 1
        j1 += 1
    if Stot == 0:
        return 0.0
    D2_this = DLine2(Term1, J1, Term2, J2)
    return Smul * D2_this / Stot


# ════════════════════════════════════════════════════════════════
#  SECTION 6: NIST Accuracy Codes
# ════════════════════════════════════════════════════════════════

_ACC_TO_PCNT = {
    "AAA": 0.3, "AA": 1, "A+": 2, "A": 3,
    "AA'": 4, "A'": 5,
    "B+": 7, "B+'": 9, "B": 10, "B'": 12,
    "C+": 18, "C": 25,
    "D+": 44, "D": 54, "E": 60,
}

_PCNT_THRESHOLDS = [
    (0.3, "AAA"), (1, "AA"), (2, "A+"), (3, "A"),
    (7, "B+"), (10, "B"), (18, "C+"), (25, "C"),
    (44, "D+"), (54, "D"),
]

def TPAcc(unc_pcnt):
    """NIST accuracy code from percentage uncertainty."""
    u1 = round(unc_pcnt, 1)
    for threshold, code in _PCNT_THRESHOLDS:
        if u1 <= threshold:
            return code
    return "E"

def TPUnc(Acc):
    """Percentage uncertainty from NIST accuracy code."""
    return _ACC_TO_PCNT.get(Acc, float('inf'))

def AccPcnt(Acc):
    """Alias for TPUnc."""
    return TPUnc(Acc)


# ════════════════════════════════════════════════════════════════
#  SECTION 7: Utilities
# ════════════════════════════════════════════════════════════════

def Num_digits(value):
    """
    Number of decimal digits after the point (positive) or
    trailing zeros before the point (negative).
    """
    s = str(value)
    dot = s.find(".")
    if dot >= 0:
        return len(s) - dot - 1
    else:
        L = len(s)
        k = L
        for j in range(L - 1, -1, -1):
            k = j + 1
            if s[j] != "0":
                break
        return k - L  # negative or zero

def NumHead(s):
    """Extract the leading numeric part of a string (including one decimal point)."""
    result = []
    has_dot = False
    for ch in s:
        if ch.isdigit():
            result.append(ch)
        elif ch == "." and not has_dot:
            has_dot = True
            result.append(ch)
        else:
            break
    return "".join(result)


# ════════════════════════════════════════════════════════════════
#  SECTION 8: Rowland Correction (data-driven)
# ════════════════════════════════════════════════════════════════

# Piecewise-constant correction table: (wl_lower, wl_upper, dwl)
# Polynomial segments are stored as (wl_lower, wl_upper, coefficients)
# where coefficients are [c6, c5, c4, c3, c2, c1, c0] for c5*wl^5 + ...
_ROWLAND_POLY_SEGMENTS = [
    (2000, 2931.7, [0, 0, 0, 0, -0.000000229357809, 0.001153419511, -1.54342646]),
    (2931.7, 3590, [0, -5.358313633E-15, 8.820237734E-11, -5.793201993E-07, 0.001897682093, -3.100083131, 2020.307718]),
    (3590, 3650.5, [0, 0, 0, 0, 0, 0, -0.141]),
    (3650, 3940.5, [0, 0, 0, 0, 0, 0, -0.14]),
    (3888, 3940.5, [0, 0, 0, 0, 0, 0, -0.14]),
    (3940.5, 4056, [8.83360197187633E-14, -2.11568654125636E-09, 2.11126530154593E-05, -0.112362706002819, 336.366854096498, -537021.84910787, 357230424.439763]),
    (4056, 4072, [0, 0, 0, 0, 0, 0, -0.155]),
    (4072, 4100, [0, -5.08067913066422E-09, 1.03769953086077E-04, -0.847774478073107, 3463.04338933707, -7073012.13982248, 5778428387.17827]),
    (4100, 4138, [0, -1.62592154881935E-09, 3.34809765199679E-05, -0.275775245790233, 1135.74567989595, -2338705.61027099, 1926321006.51177]),
    (4138, 4144, [0, 0, 0, 1.22737964115216E-04, -1.52471086564949, 6313.5645726071, -8714456.7177055]),
    (4144, 4329, [0, -1.30130464198497E-12, 2.75886195605215E-08, -2.33949807012531E-04, 0.991898000188946, -2102.62969692509, 1782788.49969217]),
    (4329, 4433, [0, 1.50387935395569e-11, -3.29731211519663e-07, 2.89171521793408e-03, -12.6797124574071, 27798.5229378821, -24377123.9784076]),
    (4433, 4565, [0, -4.26142369699592E-12, 9.58904496748238E-08, -8.63067238517732E-04, 3.88393982009225, -8738.94235832124, 7864908.92871308]),
    (4565, 4764, [0, 6.11528860366028E-13, -1.42797456974253E-08, 1.33362446233216E-04, -0.622682343536881, 1453.51593946411, -1357011.90601287]),
    (4764, 4921, [0, 9.28782160726723E-14, -2.23660028487117E-09, 2.1532269352303E-05, -0.103590885340315, 249.044354597247, -239352.727184225]),
]

# For the large piecewise-constant region (4921..7142+), store as sorted breakpoints
# Format: (wl_start, dwl)
_ROWLAND_CONST_SEGMENTS = [
    (4921, -0.171), (5050.5, -0.17), (5067.5, -0.169), (5080.5, -0.168),
    (5095.5, -0.167), (5110.5, -0.166), (5131.5, -0.165), (5160.5, -0.164),
    (5175.5, -0.163), (5188.5, -0.162), (5200.5, -0.161),
    # 5212.5–5223.5 is linear interpolation (handled specially)
    (5223.5, -0.155), (5225.5, -0.156), (5230.5, -0.157), (5234.5, -0.158),
    (5238.5, -0.159), (5242.5, -0.16), (5246.5, -0.161), (5250.5, -0.162),
    (5254.5, -0.163), (5258.5, -0.164), (5262.5, -0.165), (5266.5, -0.166),
    (5271.5, -0.167), (5275.5, -0.168), (5280.5, -0.169), (5286.5, -0.17),
    (5291.5, -0.171), (5297.5, -0.172), (5302.5, -0.173), (5307.5, -0.174),
    (5311.5, -0.175), (5316.5, -0.176), (5320.5, -0.177), (5324.5, -0.178),
    (5329.5, -0.179), (5333.5, -0.18), (5338.5, -0.181), (5342.5, -0.182),
    (5346.5, -0.183), (5350.5, -0.184), (5354.5, -0.185), (5358.5, -0.186),
    (5361.5, -0.187), (5365.5, -0.188), (5369.5, -0.189), (5373.5, -0.19),
    (5376.5, -0.191), (5379.5, -0.192), (5383.5, -0.193), (5386.5, -0.194),
    (5389.5, -0.195), (5393.5, -0.196), (5396.5, -0.197), (5399.5, -0.198),
    (5402.5, -0.199), (5405.5, -0.2), (5408.5, -0.201), (5412.5, -0.202),
    (5422.5, -0.203), (5438.5, -0.204), (5487.5, -0.205), (5496.5, -0.206),
    (5502.5, -0.207), (5506.5, -0.208), (5511.5, -0.209), (5515.5, -0.21),
    (5519.5, -0.211), (5524.5, -0.212), (5528.5, -0.213), (5532.5, -0.214),
    (5537.5, -0.215), (5543.5, -0.216), (5565.5, -0.217), (5604.5, -0.218),
    (5650.5, -0.217), (5660.5, -0.216), (5681.5, -0.215), (5700.5, -0.214),
    (5725.5, -0.213), (5737.5, -0.212), (5748.5, -0.211), (5762.5, -0.21),
    (5772.5, -0.209),
    # 5795.5–5805 linear interpolation
    (5805, -0.214), (5820, -0.213), (5871, -0.214), (5883, -0.215),
    (5895, -0.216), (5938, -0.217), (5973, -0.218), (5992, -0.217),
    (6008, -0.216), (6025, -0.215), (6045, -0.214), (6067, -0.213),
    (6085, -0.212), (6105, -0.211), (6124, -0.21), (6163, -0.211),
    (6177, -0.212), (6190, -0.211), (6197, -0.212), (6217, -0.213),
    (6236, -0.214), (6275, -0.215), (6295, -0.216), (6326, -0.215),
    (6338, -0.214), (6362, -0.213), (6369, -0.212), (6376, -0.211),
    (6445, -0.212), (6454, -0.213), (6460, -0.214), (6464, -0.215),
    (6468, -0.216), (6474, -0.217), (6478, -0.218), (6483, -0.219),
    (6487, -0.22), (6491, -0.221), (6496, -0.222), (6500, -0.223),
    (6505, -0.224), (6509, -0.225), (6513, -0.226), (6517, -0.227),
    (6520, -0.228), (6523, -0.229), (6525, -0.23), (6529, -0.231),
    (6532, -0.232), (6535, -0.233), (6538, -0.234), (6542, -0.235),
    (6548, -0.236), (6554, -0.237), (6567, -0.238), (6582, -0.239),
    (6593, -0.24), (6605, -0.241), (6618, -0.242), (6632, -0.243),
    (6644, -0.244), (6657, -0.245), (6670, -0.246), (6682, -0.247),
    (6696, -0.248), (6719, -0.249), (6739, -0.25), (6768, -0.251),
    (6795, -0.252), (6823, -0.253), (6848, -0.254),
    # 6868–6870 linear interpolation
    (6870, -0.242), (6873, -0.243), (6880, -0.244), (6887, -0.245),
    (6896, -0.246), (6908, -0.247), (6914, -0.248), (6921, -0.249),
    (6930, -0.25), (6938, -0.251), (6944, -0.252), (6953, -0.253),
    (6960, -0.254), (6965, -0.255), (6970, -0.256), (6975, -0.257),
    (6981, -0.258), (6986, -0.259), (6991, -0.26), (6997, -0.261),
    (7002, -0.262), (7007, -0.263), (7011, -0.264), (7017, -0.265),
    (7023, -0.266), (7029, -0.267), (7037, -0.268), (7050, -0.269),
    (7062, -0.27), (7074, -0.271), (7086, -0.272), (7097, -0.273),
    (7108, -0.274), (7119, -0.275), (7130, -0.276), (7142, -0.277),
]

def RowlandCorr(wl):
    """
    Rowland wavelength correction (Å) as a function of air wavelength.
    Returns the systematic correction dwl to be *added* to the Rowland wavelength.
    """
    if wl < 2000 or wl >= 7142:
        # Outside the polynomial range, use the terminal constant values
        if wl >= 7142:
            return -0.277
        return 0.0

    # Check polynomial segments first (wl < 4921)
    if wl < 4921:
        for (lo, hi, coeffs) in _ROWLAND_POLY_SEGMENTS:
            if lo <= wl < hi:
                if all(c == 0 for c in coeffs):
                    return 0.0
                return sum(c * wl ** (6 - i) for i, c in enumerate(coeffs))
        return 0.0

    # Special linear interpolation regions
    if 5212.5 <= wl < 5223.5:
        return -0.161 + 0.006 * (wl - 5212.5) / (5223.5 - 5212.5)
    if 5795.5 <= wl < 5805:
        return -0.209 - 0.005 * (wl - 5795.5) / (5805 - 5795.5)
    if 6868 <= wl < 6870:
        return -0.254 + 0.012 * (wl - 6868) / 2.0

    # Piecewise constant lookup
    import bisect
    breakpoints = [seg[0] for seg in _ROWLAND_CONST_SEGMENTS]
    idx = bisect.bisect_right(breakpoints, wl) - 1
    if 0 <= idx < len(_ROWLAND_CONST_SEGMENTS):
        return _ROWLAND_CONST_SEGMENTS[idx][1]
    return 0.0

