"""
Functions for computing plasma density profile and laser frequency for dephasingless Plasma Wakefield Photon Acceleration (PWPA)

Written by Ryan Sandberg
This version August 2023

Contains
---
get_wake_at_driver_end_when_nd_lt_1_2_np
get_wake_at_driver_end_when_nd_eq_1_2_np
get_wake_at_driver_end_when_nd_gt_1_2_np
get_gamma_max
get_tau_where_delta_n_is_0
get_pwpa_profile


Notes
---
For background, refer to following papers:
dePWPA demonstration [R. T. Sandberg and A. G. R. Thomas, "Photon acceleration
from optical to xuv," Phys. Rev. Lett. 130 (2023), 10.1103/Phys-
RevLett.130.085001.]

Original photon acceleration paper [S. C. Wilks, J. M. Dawson, W. B. Mori, T. Katsouleas,
and M. E. Jones, "Photon accelerator", Physical Review Letters 62 (1989)]

Wake solutions [J. B. Rosenzweig, "Nonlinear plasma dynamics in the plasma wake-field
accelerator", Physical Review Letters 58, 555-58 (1987).]

Sample usage
---

>>> import numpy as np
>>> from matplotlib import pyplot as plt

>>> import dephasingless_pwpa_plasma_profile
>>> pwpa = dephasingless_pwpa_plasma_profile

>>> nc0_np0 = 100
>>> nd0_np0 = 2
>>> kp0Ld0 = 0.25
>>> lambda_L0 = 8e-7
>>> lambda_p0 = lambda_L0 * np.sqrt(nc0_np0)
>>> kp0 = 2*np.pi / lambda_p0
>>> kp0z, n_n0_pwpa, k_kp0_pwpa = pwpa.get_pwpa_profile(nc0_np0, nd0_np0, kp0Ld0,z_final=80000,n_steps=160001)
"""

import numpy as np
from scipy.constants import c, e, m_e, epsilon_0
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import scipy.optimize as opt
optimize = opt
from scipy import special as sp


def get_wake_at_driver_end_when_nd_lt_1_2_np(nd0_np, kp0Ld, kp, bracket=[-20,20]):
    """
    Calculate wake parameters at drive end when drive density < 1/2 plasma density

    Parameters
    ---
    nd0_np : float, ratio of drive beam density to plasma density
    kp0Ld : float, drive beam length normalized to reference plasma wavenumber kp0
    kp : float, plasma wavenumber
    bracket: range of wave phase phi for numerical root finder

    Returns
    ---
    x_sub_drive_end: float, Rosenzweig's wake parameter x at end of drive beam
    sgn_E_drive_end: float, sign of E at the end of the drive beam
    """
    ksq = 2 * nd0_np
    def phi_eqn(phi):
        return (1-ksq) * kp*kp0Ld - 2 * (sp.ellipe(ksq) - sp.ellipeinc(phi, ksq))
    solns = opt.root_scalar(phi_eqn, bracket = bracket)
    phi_drive_end = solns.root

    x_sub_drive_end = 1 + ksq /(1-ksq) * np.cos(phi_drive_end)**2
    sgn_E_drive_end = np.sign(np.sin(2*phi_drive_end))
    return x_sub_drive_end, sgn_E_drive_end

def get_wake_at_driver_end_when_nd_eq_1_2_np(nd0_np, kp0Ld, kp, bracket=[1,100]):
    """
    Calculate wake parameters at drive end when drive density = 1/2 plasma density

    Parameters
    ---
    nd0_np : float, ratio of drive beam density to plasma density
    kp0Ld : float, drive beam length normalized to reference plasma wavenumber kp0
    kp : float, plasma wavenumber
    bracket: range of Rosenzweig's wake parameter x for numerical root finder

    Returns
    ---
    x_sub_drive_end: float, Rosenzweig's wake parameter x at end of drive beam
    sgn_E_drive_end: float, sign of E at the end of the drive beam
    """
    def x_eqn(x):
        return np.sqrt(x)*np.sqrt(x-1) + np.arcsinh(np.sqrt(x-1)) - kp*kp0Ld
    solns = opt.root_scalar(x_eqn)

    x_sub_drive_end = solns.root
    sgn_E_drive_end = 1.0
    return x_sub_drive_end, sgn_E_drive_end


def get_wake_at_driver_end_when_nd_gt_1_2_np(nd0_np, kp0Ld, kp, bracket=[-np.pi/2,0]):
    """
    Calculate wake parameters at drive end when drive density > 1/2 plasma density

    Parameters
    ---
    nd0_np : float, ratio of drive beam density to plasma density
    kp0Ld : float, drive beam length normalized to reference plasma wavenumber kp0
    kp : float, plasma wavenumber
    bracket: range of wave phase phi for numerical root finder

    Returns
    ---
    x_sub_drive_end: float, Rosenzweig's wake parameter x at end of drive beam
    sgn_E_drive_end: float, sign of E at the end of the drive beam
    """
    consta = np.sqrt(2*nd0_np)
    constb = np.sqrt(2*nd0_np - 1)
    kappa_sq = 1 / consta**2
    kappap_sq = constb**2 / consta**2

    def phi_eqn(phi):
        LHS = kappap_sq * sp.ellipkinc(phi,kappa_sq)
        LHS -= sp.ellipeinc(phi, kappa_sq)
        LHS += np.tan(phi) * np.sqrt(1 - kappa_sq*np.sin(phi)**2)
        LHS *= 2*consta
        LHS += kp * kp0Ld*constb**2
        return LHS

    solns = opt.root_scalar(phi_eqn, bracket=bracket)

    phi_drive_end = solns.root
    x_sub_drive_end = 1 + np.tan(phi_drive_end)**2
    sgn_E_drive_end = 1.0
    return x_sub_drive_end, sgn_E_drive_end

def get_gamma_max(nd0,kp0Ld,n_p):
    """
    Find wake amplitude behind beam driver

    Parameters
    ---
    nd0 : float, driver density
    kp0Ld : float, driver length (k_p0 Ld)
    n_p : float, plasma density

    Returns
    ---
    gamma_max : float, maximum gamma attained in wave behind dr
    gamma_drive_end : float, fluid gamma at the end of the driver
    x_sub_drive_end : float, Rosenzweig's `x` = potential at the end of the driver
    p_drive_end :
    E_drive_end :
    """

    nd0_np = nd0 / n_p
    kp = np.sqrt(n_p)
    if (nd0_np < 0.5):
        x_sub_drive_end, sgn_E_drive_end = get_wake_at_driver_end_when_nd_lt_1_2_np(nd0_np, kp0Ld, kp)
    elif nd0_np == 0.5:
        x_sub_drive_end, sgn_E_drive_end = get_wake_at_driver_end_when_nd_eq_1_2_np(nd0_np, kp0Ld, kp)
    elif nd0_np > 0.5:
        x_sub_drive_end, sgn_E_drive_end = get_wake_at_driver_end_when_nd_gt_1_2_np(nd0_np, kp0Ld, kp)

    else:
        raise TypeError('nd0_np must be nonnegative float')
    abs_dx_dtau = np.sqrt(2*(1-nd0_np) - 1./x_sub_drive_end + (2*nd0_np - 1) * x_sub_drive_end)
    E_drive_end = sgn_E_drive_end * kp*abs_dx_dtau
    beta_drive_end = (1 - x_sub_drive_end**2) / (1 + x_sub_drive_end**2)
    gamma_drive_end = (1 + x_sub_drive_end**2) / 2 / x_sub_drive_end
    p_drive_end = beta_drive_end * gamma_drive_end
    gamma_max = 1 - nd0_np + nd0_np * x_sub_drive_end
    return gamma_max,gamma_drive_end,p_drive_end,E_drive_end


v_gamma = np.vectorize(get_gamma_max)

def get_tau_where_delta_n_is_0(nd0_np0, kp0Ld0, n_p, wave_beta = 1.0):
    """
    Calculate the position `tau` in the wake where density perturbation dn goes to 0 -- the phase to match to for PWPA

    Parameters
    ---
    nd0_np : float, ratio of drive beam density to plasma density
    kp0Ld : float, drive beam length normalized to reference plasma wavenumber kp0
    n_p : float or array-like, plasma density
    wave_beta : float, the speed of the wake normalized to the speed of light `c`; default is 1, assuming ultra-relativistic driver

    Returns
    ---
    tau_dn0 : float or array-like, depending on input n_p; the position `tau` where density perturbation `delta n` is 0
    """
    kp = np.sqrt(n_p)
    gamma_max,gamma_drive_end,p_drive_end,E_drive_end = v_gamma(nd0_np0,kp0Ld0,n_p)
    q = -1
    qE_drive_end_kp = E_drive_end / kp
    beta_drive_end = p_drive_end / gamma_drive_end
    nb = n_p / (1-beta_drive_end)
    phi_wake = np.arctan2(qE_drive_end_kp * np.sqrt(1+gamma_drive_end), p_drive_end*np.sqrt(2))
    phi_wake_unwrap = np.unwrap(phi_wake * 2) / 2
    phi_wake_pi = phi_wake_unwrap / np.pi
    next_int = np.ceil(phi_wake_pi)
    next_odd_int = next_int + np.mod(next_int + 1,2)
    phi_dn0 = np.pi * (next_odd_int + 0.5)

    kappa_sq = (gamma_max-1.) / (gamma_max + 1.)
    kappa = np.sqrt(kappa_sq)
    kappa_p = np.sqrt(2 / (gamma_max + 1.))

    tau_b = wave_beta * 1/kp*(2 / kappa_p * sp.ellipeinc(phi_wake_unwrap, kappa_sq) - kappa_p * sp.ellipkinc(phi_wake_unwrap, kappa_sq))
    d_tau_b = -1/kp * 2 * kappa / kappa_p * np.sin(phi_wake_unwrap)
    tau_b= tau_b + d_tau_b

    tau_dn0 = wave_beta * 1/kp*(2 / kappa_p * sp.ellipeinc(phi_dn0, kappa_sq) - kappa_p * sp.ellipkinc(phi_dn0, kappa_sq))
    d_tau_dn0 = -1/kp * 2 * kappa / kappa_p * np.sin(phi_dn0)
    tau_dn0 = tau_dn0 + d_tau_dn0 - tau_b+kp0Ld0
    return tau_dn0


def get_pwpa_profile(nc0_np0,
                     nd0_np0,
                     kp0Ld0,
                     n_n0_min = 0.01,
                     n_n0_max = 1.1,
                     kp0z_initial = 0,
                     kp0z_final = 50000,
                     n_steps = 100001):
    """
    Calculate the Plasma wakefield photon acceleration (PWPA) plasma density profile and theoretical laser frequency

    Parameters
    ---
    nc0_np0: float, ratio of initial laser critical density to plasma density or ratio of laser frequency squared to plasma frequency squared
    nd0_np0: float, ratio of initial drive beam density to initial plasma density
    kp0Ld0: float, drive beam length normalized to initial plasma wavenumber
    n_n0_min: minimum possible plasma density ratio that can be computed
    n_n0_max: maximum possible plasma density ratio that can be computed
    kp0z_initial: float, start time/position of PWPA profile
    kp0z_final: float, final time/position of PWPA profile
    n_steps: number of integration steps in numerical solving of the PWPA coupled differential equations

    Returns
    ---
    kp0z : numpy array

    """
    n_p_interp = np.linspace(n_n0_min,n_n0_max)
    tau_dn0_list = get_tau_where_delta_n_is_0(nd0_np0, kp0Ld0, n_p_interp)
    dn = n_p_interp[1]-n_p_interp[0]
    dtaudn0_dn = (tau_dn0_list[2:] - tau_dn0_list[:-2]) / 2 / dn

    dtau_dn0_dn_fun = interp1d(n_p_interp[1:-1], dtaudn0_dn)

    def pa_coupled_nonlinear(y,z, nd0_np0, kp0Ld0):
        n_p, nc = y

        dxi_dn = -dtau_dn0_dn_fun(n_p)
        dnp_dz = -n_p/2/nc / dxi_dn

        gamma_max,_,_,_ = v_gamma(nd0_np0, kp0Ld0, n_p)
        dnc_dz = np.sqrt(2)*n_p**1.5 * np.sqrt(gamma_max - 1)
        return [dnp_dz, dnc_dz]


    kp0z = np.linspace(kp0z_initial, kp0z_final, n_steps)
    sol = odeint(pa_coupled_nonlinear, [1,nc0_np0], kp0z, args=(nd0_np0, kp0Ld0))
    n_n0_pwpa = sol[:,0]
    k_kp0_pwpa_sq = sol[:,1]
    return kp0z, n_n0_pwpa, k_kp0_pwpa_sq