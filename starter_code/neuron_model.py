"""Hodgkin-Huxley Cortical Pyramidal Axon Model for TMS Stimulation.

This shared module provides a calibrated, unconditionally stable pure-NumPy
Hodgkin-Huxley membrane simulator (using Rush-Larsen exponential integration)
for Step 3 of the EEC2204 Term Project. Electromagnetic field calculations
(Steps 1-2) supply the time-varying local electric field E_eff(t) [V/m] along
the neuron's polarization axis, which drives the transmembrane depolarizing
current density:

    I_stim(t) = g_couple * E_eff(t)   [uA/cm^2]

where g_couple = 5.2 (uA/cm^2)/(V/m) represents the effective sub-millimeter
axon-bend / terminal polarization coupling factor and phi_T = 3.0 models
Q10-scaled cortical ion-channel kinetics. Under a 2.5 kHz damped-sinusoidal
TMS pulse, the all-or-none Action Potential threshold occurs near E_th ~ 96 V/m,
matching the canonical ~100 V/m cortical stimulation threshold.
"""

from __future__ import annotations

import math
import numpy as np

# Default Hodgkin-Huxley biophysical parameters
DEFAULT_CM = 1.0          # Membrane capacitance [uF/cm^2]
DEFAULT_G_NA = 120.0      # Peak Na+ conductance [mS/cm^2]
DEFAULT_G_K = 36.0        # Peak K+ conductance [mS/cm^2]
DEFAULT_G_L = 0.3         # Leak conductance [mS/cm^2]
DEFAULT_E_NA = 50.0       # Na+ reversal potential [mV]
DEFAULT_E_K = -77.0       # K+ reversal potential [mV]
DEFAULT_V_REST = -65.0    # Resting membrane potential [mV]
DEFAULT_PHI_T = 3.0       # Q10 temperature scaling factor for cortical kinetics
DEFAULT_G_COUPLE = 5.2    # E-field to transmembrane current coupling [(uA/cm^2)/(V/m)]
MAX_DT_MS = 0.002         # Maximum time step [ms] (2 us)


def _clamp(val: float, lo: float, hi: float) -> float:
    return lo if val < lo else (hi if val > hi else val)


def _hh_rates(V: float) -> tuple[float, float, float, float, float, float]:
    """Compute Hodgkin-Huxley gating transition rates alpha and beta [1/ms]."""
    dV_m = V + 40.0
    if abs(dV_m) < 1e-6:
        alpha_m = 1.0
    else:
        alpha_m = 0.1 * dV_m / (1.0 - math.exp(_clamp(-dV_m / 10.0, -50.0, 50.0)))
    beta_m = 4.0 * math.exp(_clamp(-(V + 65.0) / 18.0, -50.0, 50.0))

    alpha_h = 0.07 * math.exp(_clamp(-(V + 65.0) / 20.0, -50.0, 50.0))
    beta_h = 1.0 / (1.0 + math.exp(_clamp(-(V + 35.0) / 10.0, -50.0, 50.0)))

    dV_n = V + 55.0
    if abs(dV_n) < 1e-6:
        alpha_n = 0.1
    else:
        alpha_n = 0.01 * dV_n / (1.0 - math.exp(_clamp(-dV_n / 10.0, -50.0, 50.0)))
    beta_n = 0.125 * math.exp(_clamp(-(V + 65.0) / 80.0, -50.0, 50.0))

    return alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n


def _initial_steady_state(v_rest: float) -> tuple[float, float, float, float]:
    """Return resting gating variables (m0, h0, n0) and equilibrium leak potential E_L."""
    am, bm, ah, bh, an, bn = _hh_rates(v_rest)
    m0 = am / (am + bm)
    h0 = ah / (ah + bh)
    n0 = an / (an + bn)
    i_na0 = DEFAULT_G_NA * (m0 ** 3) * h0 * (v_rest - DEFAULT_E_NA)
    i_k0 = DEFAULT_G_K * (n0 ** 4) * (v_rest - DEFAULT_E_K)
    e_l = v_rest + (i_na0 + i_k0) / DEFAULT_G_L
    return m0, h0, n0, e_l


def simulate_hh_neuron(
    t_ms: np.ndarray,
    E_eff_Vm: np.ndarray,
    g_couple: float = DEFAULT_G_COUPLE,
    phi_T: float = DEFAULT_PHI_T,
    v_rest: float = DEFAULT_V_REST,
    ap_threshold_mV: float = 0.0,
) -> dict[str, np.ndarray | float | bool | None]:
    """Integrate the Hodgkin-Huxley equations under a TMS electric-field waveform.

    Parameters
    ----------
    t_ms : np.ndarray
        Monotonically increasing 1D array of time samples in milliseconds [ms].
    E_eff_Vm : np.ndarray
        Effective depolarizing electric field along the neuron axis [V/m],
        sampled at `t_ms`. Positive values depolarize the membrane.
    g_couple : float
        Field-to-current coupling coefficient [(uA/cm^2)/(V/m)].
    phi_T : float
        Temperature rate-scaling factor (dimensionless).
    v_rest : float
        Resting membrane potential [mV].
    ap_threshold_mV : float
        Peak membrane voltage threshold to classify an Action Potential [mV].

    Returns
    -------
    dict
        Dictionary containing `t_ms`, `V_m`, `m`, `h`, `n`, `I_stim`, `I_Na`,
        `I_K`, `I_L`, `peak_V_m`, `min_V_m`, `ap_triggered`, and `ap_latency_ms`.
    """
    t_in = np.asarray(t_ms, dtype=float)
    e_in = np.asarray(E_eff_Vm, dtype=float)
    if t_in.ndim != 1 or e_in.shape != t_in.shape or len(t_in) < 2:
        raise ValueError("t_ms and E_eff_Vm must be 1D arrays of identical length >= 2.")
    if np.any(np.diff(t_in) <= 0.0):
        raise ValueError("t_ms must be strictly increasing.")

    # Refine onto a uniform grid with dt <= MAX_DT_MS for high temporal accuracy
    t_span = t_in[-1] - t_in[0]
    n_sub = max(int(math.ceil(t_span / MAX_DT_MS)), len(t_in) - 1)
    t_fine = np.linspace(t_in[0], t_in[-1], n_sub + 1)
    dt = float(t_fine[1] - t_fine[0])
    e_fine = np.interp(t_fine, t_in, e_in)
    i_stim_list = (g_couple * e_fine).tolist()

    m0, h0, n0, e_l = _initial_steady_state(v_rest)
    v_arr = np.empty_like(t_fine)
    m_arr = np.empty_like(t_fine)
    h_arr = np.empty_like(t_fine)
    n_arr = np.empty_like(t_fine)

    V, m, h, n = float(v_rest), m0, h0, n0
    v_arr[0], m_arr[0], h_arr[0], n_arr[0] = V, m, h, n

    # Rush-Larsen exponential Euler: unconditionally stable for any pulse intensity
    for k in range(n_sub):
        i_mid = 0.5 * (i_stim_list[k] + i_stim_list[k + 1])
        g_na = DEFAULT_G_NA * (m ** 3) * h
        g_k = DEFAULT_G_K * (n ** 4)
        g_tot = g_na + g_k + DEFAULT_G_L
        v_inf = (i_mid + g_na * DEFAULT_E_NA + g_k * DEFAULT_E_K + DEFAULT_G_L * e_l) / g_tot
        V = v_inf + (V - v_inf) * math.exp(-dt * g_tot / DEFAULT_CM)

        am, bm, ah, bh, an, bn = _hh_rates(V)
        r_m = am + bm
        m_inf = am / r_m
        m = m_inf + (m - m_inf) * math.exp(-phi_T * r_m * dt)

        r_h = ah + bh
        h_inf = ah / r_h
        h = h_inf + (h - h_inf) * math.exp(-phi_T * r_h * dt)

        r_n = an + bn
        n_inf = an / r_n
        n = n_inf + (n - n_inf) * math.exp(-phi_T * r_n * dt)

        v_arr[k + 1] = V
        m_arr[k + 1] = m
        h_arr[k + 1] = h
        n_arr[k + 1] = n

    # Sample trajectories back onto the caller's t_ms grid
    V_out = np.interp(t_in, t_fine, v_arr)
    m_out = np.interp(t_in, t_fine, m_arr)
    h_out = np.interp(t_in, t_fine, h_arr)
    n_out = np.interp(t_in, t_fine, n_arr)
    I_stim_out = g_couple * e_in
    I_Na_out = DEFAULT_G_NA * (m_out ** 3) * h_out * (V_out - DEFAULT_E_NA)
    I_K_out = DEFAULT_G_K * (n_out ** 4) * (V_out - DEFAULT_E_K)
    I_L_out = DEFAULT_G_L * (V_out - e_l)

    peak_V_m = float(np.max(v_arr))
    min_V_m = float(np.min(v_arr))
    ap_triggered = bool(peak_V_m >= ap_threshold_mV)
    if ap_triggered:
        crossing_idx = int(np.argmax(v_arr >= ap_threshold_mV))
        ap_latency_ms: float | None = float(t_fine[crossing_idx])
    else:
        ap_latency_ms = None

    return {
        "t_ms": t_in,
        "V_m": V_out,
        "m": m_out,
        "h": h_out,
        "n": n_out,
        "I_stim": I_stim_out,
        "I_Na": I_Na_out,
        "I_K": I_K_out,
        "I_L": I_L_out,
        "peak_V_m": peak_V_m,
        "min_V_m": min_V_m,
        "ap_triggered": ap_triggered,
        "ap_latency_ms": ap_latency_ms,
    }


def find_electric_field_threshold(
    t_ms: np.ndarray,
    unit_E_waveform: np.ndarray,
    g_couple: float = DEFAULT_G_COUPLE,
    phi_T: float = DEFAULT_PHI_T,
    v_rest: float = DEFAULT_V_REST,
    ap_threshold_mV: float = 0.0,
    low_Vm: float = 1.0,
    high_Vm: float = 300.0,
    tol_Vm: float = 0.05,
) -> float:
    """Find the minimum peak electric field amplitude [V/m] triggering an AP.

    Parameters
    ----------
    t_ms : np.ndarray
        Time grid in milliseconds [ms].
    unit_E_waveform : np.ndarray
        Normalized electric-field pulse shape (max(|unit_E_waveform|) == 1.0).
    """
    lo = float(low_Vm)
    hi = float(high_Vm)
    if simulate_hh_neuron(t_ms, lo * unit_E_waveform, g_couple, phi_T, v_rest, ap_threshold_mV)["ap_triggered"]:
        return lo
    while not simulate_hh_neuron(t_ms, hi * unit_E_waveform, g_couple, phi_T, v_rest, ap_threshold_mV)["ap_triggered"]:
        hi *= 2.0
        if hi > 1e5:
            raise RuntimeError("Neuron did not fire even at 100,000 V/m.")

    while (hi - lo) > tol_Vm:
        mid = 0.5 * (lo + hi)
        res = simulate_hh_neuron(t_ms, mid * unit_E_waveform, g_couple, phi_T, v_rest, ap_threshold_mV)
        if res["ap_triggered"]:
            hi = mid
        else:
            lo = mid
    return float(0.5 * (lo + hi))


if __name__ == "__main__":
    t_test = np.linspace(0.0, 8.0, 4001)
    omega_ms = 2.0 * math.pi * 2.5
    alpha_ms = 3.0
    unit_wave = np.exp(-alpha_ms * t_test) * (np.cos(omega_ms * t_test) - (alpha_ms / omega_ms) * np.sin(omega_ms * t_test))
    eth = find_electric_field_threshold(t_test, unit_wave)
    sub = simulate_hh_neuron(t_test, 0.8 * eth * unit_wave)
    sup = simulate_hh_neuron(t_test, 1.2 * eth * unit_wave)
    print(f"HH self-test | E_th = {eth:.2f} V/m")
    print(f"  0.8*E_th: peak V_m = {sub['peak_V_m']:+.2f} mV, AP = {sub['ap_triggered']}")
    print(f"  1.2*E_th: peak V_m = {sup['peak_V_m']:+.2f} mV, AP = {sup['ap_triggered']}, latency = {sup['ap_latency_ms']:.3f} ms")
