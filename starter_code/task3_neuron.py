"""Task 3 — 대뇌 피질 뉴런 활성화(cortical neuron activation): Task 3.1.

피질의 유도 전기장을 제공된 호지킨-헉슬리(Hodgkin–Huxley) 뉴런 모델
(neuron_model.py — 수정 금지)에 입력합니다. Task 1 함수가 필요합니다.

    python task3_neuron.py

를 실행하면 Task 3.1 결과가 출력됩니다.
"""

from __future__ import annotations

import numpy as np

from neuron_model import find_electric_field_threshold, simulate_hh_neuron
from tms_config import (
    DIDT_RISE,
    I_PEAK_NOMINAL,
    R_OFF_TARGET,
    R_TARGET,
    T_MS,
    T_RISE_S,
    T_S,
    U_NEURON,
)
from task1_coil_field import (
    build_circular_coil,
    build_figure8_coil,
    compute_inductance,
    compute_unit_A_and_B,
    trapezoidal_pulse_current_and_didt,
)


# ==============================================================================
# Task 3.1 — HH 뉴런 결합(coupling), 발화 선택성(selectivity), 임계 구동 조건(critical threshold)
# ==============================================================================
def compute_E_eff_waveform(
    pts: np.ndarray,
    dl: np.ndarray,
    r_obs: np.ndarray,
    di_dt: np.ndarray,
) -> np.ndarray:
    """유도 전기장의 축삭 방향(axon direction) 성분 E_eff(t)를 계산합니다.

    E_eff(t) = E(r_obs, t) · U_NEURON   (E는 Task 2.1에서 유도한 유도 전기장 식으로 계산)

    인자(Parameters)
    ----------------
    pts, dl : ndarray, shape (K, 3)
        코일 선분 중점과 현 벡터 [m].
    r_obs : ndarray, shape (3,)
        뉴런 위치 [m].
    di_dt : ndarray, shape (T,)
        전류의 시간 미분 dI/dt [A/s].

    반환값(Returns)
    ---------------
    e_eff : ndarray, shape (T,)
        부호를 가진 축삭 방향 전기장 E_eff(t) [V/m]. simulate_hh_neuron(T_MS, e_eff)에 그대로 입력합니다.
    """
    # === TODO [Task 3.1a] ===
    raise NotImplementedError("TODO [Task 3.1a]: compute_E_eff_waveform()을 구현하십시오.")


def compute_fig8_thresholds(
    e_th: float,
    e_target_peak: float,
    l_fig8: float,
) -> dict:
    """발화 역치(firing threshold)에 해당하는 8자형 코일의 임계 전류·전압·에너지를 계산합니다 (T_rise = 0.06 ms 고정).

    인자(Parameters)
    ----------------
    e_th : float
        뉴런 발화 역치 E_th [V/m] (find_electric_field_threshold의 결과).
    e_target_peak : float
        명목 6 kA에서 표적점의 상승 구간 |E_eff| [V/m].
    l_fig8 : float
        8자형 코일 인덕턴스 [H].

    반환값(Returns)
    ---------------
    thresholds : dict
        키: "I_peak_th_kA" (임계 첨두 전류 [kA]), "V_coil_th_kV" (임계 구동 전압 [kV]),
        "W_m_th_J" (임계 저장 에너지 [J]).
    """
    # === TODO [Task 3.1b] ===
    # E는 I_peak에 비례하므로 I_peak_th = I_PEAK_NOMINAL · e_th / e_target_peak.
    # 이어서 V_coil_th = L · I_peak_th / T_RISE_S, W_m_th = 0.5 · L · I_peak_th².
    raise NotImplementedError("TODO [Task 3.1b]: compute_fig8_thresholds()를 구현하십시오.")


def find_max_standoff(
    e_th: float,
    h_lo: float = 5e-3,
    h_hi: float = 60e-3,
    n_iter: int = 40,
) -> float:
    """명목 6 kA에서 표적 뉴런을 여전히 발화시키는 최대 이격 거리 h_max를 이분법(bisection)으로 찾습니다.

    인자(Parameters)
    ----------------
    e_th : float
        뉴런 발화 역치 E_th [V/m].
    h_lo, h_hi : float
        탐색 구간 [m]. h_lo에서는 발화하고 h_hi에서는 발화하지 않아야 합니다.
    n_iter : int
        이분법 반복 횟수.

    반환값(Returns)
    ---------------
    h_max : float
        상승 구간에서 E_eff(R_TARGET) ≥ e_th를 만족하는 최대 이격 거리 [m].
    """
    # === TODO [Task 3.1c] ===
    raise NotImplementedError("TODO [Task 3.1c]: find_max_standoff()를 구현하십시오.")


# ==============================================================================
# Task 3 검증 실행(이 파일을 직접 실행할 때)
# ==============================================================================
def main() -> None:
    circ_pts, circ_dl = build_circular_coil()
    fig8_pts, fig8_dl = build_figure8_coil()
    l_fig8 = compute_inductance(fig8_pts, fig8_dl)
    _, di_dt = trapezoidal_pulse_current_and_didt(T_S)

    # --- Task 3.1a: 표적 vs. 비표적 뉴런 응답 ---
    cases = {
        "Figure-8 Target    ": (fig8_pts, fig8_dl, R_TARGET),
        "Figure-8 Off-Target": (fig8_pts, fig8_dl, R_OFF_TARGET),
        "Circular Target    ": (circ_pts, circ_dl, R_TARGET),
        "Circular Off-Target": (circ_pts, circ_dl, R_OFF_TARGET),
    }
    e_fig8_target_peak = 0.0
    for name, (pts, dl, r) in cases.items():
        e_eff = compute_E_eff_waveform(pts, dl, r, di_dt)
        hh = simulate_hh_neuron(T_MS, e_eff)
        print(f"[Task 3.1a] {name} | E_peak = {np.max(np.abs(e_eff)):8.2f} V/m, "
              f"Vm_max = {hh['peak_V_m']:+.2f} mV, AP = {hh['ap_triggered']}")
        if name.startswith("Figure-8 Target"):
            e_fig8_target_peak = float(np.max(np.abs(e_eff)))

    # --- Task 3.1b/c: 임계 구동 조건 ---
    e_th = find_electric_field_threshold(T_MS, di_dt / DIDT_RISE)
    th = compute_fig8_thresholds(e_th, e_fig8_target_peak, l_fig8)
    h_max = find_max_standoff(e_th)
    print(f"[Task 3.1b] E_th = {e_th:.2f} V/m | I_peak_th = {th['I_peak_th_kA']:.3f} kA, "
          f"V_coil_th = {th['V_coil_th_kV']:.3f} kV, Wm_th = {th['W_m_th_J']:.2f} J")
    print(f"[Task 3.1c] h_max = {h_max*1e3:.2f} mm (at I_peak = {I_PEAK_NOMINAL/1e3:.1f} kA)")


if __name__ == "__main__":
    main()
