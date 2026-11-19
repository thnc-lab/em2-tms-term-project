"""Task 1 — 코일 자기장과 인덕턴스(coil magnetic field & inductance): Task 1.1, Task 1.2.

아래의 `# === TODO [Task X.Y] ===` 블록을 모두 구현한 뒤

    python task1_coil_field.py

를 실행하면 Task 1.1 / 1.2 검증 값이 출력됩니다. Task 2와 Task 3은 이 파일의
함수를 불러 쓰므로 Task 1을 먼저 완성하십시오.
"""

from __future__ import annotations

import numpy as np

from tms_config import (
    DIDT_RISE,
    H_NOMINAL,
    I_PEAK_NOMINAL,
    K_SEG,
    MU0,
    N_TURNS,
    R_CIRC,
    R_FIG8,
    R_OFF_TARGET,
    R_TARGET,
    S_FIG8,
    T_FALL_S,
    T_HOLD_S,
    T_RISE_S,
    T_S,
    WIRE_RADIUS,
)


# ==============================================================================
# Task 1.1 — 코일 형상, A·B 계산, 중심축 검증(on-axis validation)
# ==============================================================================
def make_circular_loop(
    cx: float,
    cy: float,
    z0: float,
    radius: float,
    direction_sign: float = 1.0,
    k_seg: int = K_SEG,
) -> tuple[np.ndarray, np.ndarray]:
    """수평 원형 루프(circular loop)를 k_seg개의 직선 선분(segment)으로 분할합니다.

    인자(Parameters)
    ----------------
    cx, cy, z0 : float
        루프 중심 좌표 [m]. 루프는 z = z0 평면에 놓입니다.
    radius : float
        루프 반경 [m].
    direction_sign : float
        감는 방향. +1 = 반시계방향(CCW), -1 = 시계방향(CW) (+z 쪽에서 내려다볼 때).
    k_seg : int
        선분 개수 K.

    반환값(Returns)
    ---------------
    pts : ndarray, shape (k_seg, 3)
        각 선분의 중점(midpoint) 위치 [m].
    dl : ndarray, shape (k_seg, 3)
        각 선분의 방향을 가진 현 벡터(chord vector) [m]. 크기는 현의 길이이고
        방향은 전류 방향입니다.
    """
    # === TODO [Task 1.1a] ===
    # 힌트: 중점 각도 phi = (np.arange(k_seg) + 0.5) * (2*pi / k_seg),
    #       현의 길이 = 2 * radius * sin(0.5 * dphi), 방향은 각 점의 접선(tangent) 방향.
    raise NotImplementedError("TODO [Task 1.1a]: make_circular_loop()를 구현하십시오.")


def build_circular_coil(
    radius: float = R_CIRC,
    h_standoff: float = H_NOMINAL,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    k_seg: int = K_SEG,
) -> tuple[np.ndarray, np.ndarray]:
    """(제공) 중심 (offset_x, offset_y, h_standoff)의 반시계방향 원형 코일을 만듭니다.

    인자(Parameters)
    ----------------
    radius : float
        코일 반경 [m].
    h_standoff : float
        코일 높이 = 코일–두피 이격 거리(standoff distance) h [m].
    offset_x, offset_y : float
        코일 중심의 수평 위치 [m].
    k_seg : int
        선분 개수 K.

    반환값(Returns)
    ---------------
    pts, dl : ndarray, shape (k_seg, 3)
        make_circular_loop와 같은 형식.
    """
    return make_circular_loop(offset_x, offset_y, h_standoff, radius, direction_sign=+1.0, k_seg=k_seg)


def build_figure8_coil(
    wing_radius: float = R_FIG8,
    wing_offset: float | None = None,
    h_standoff: float = H_NOMINAL,
    offset_x: float = 0.0,
    offset_y: float = 0.0,
    k_seg: int = K_SEG,
) -> tuple[np.ndarray, np.ndarray]:
    """원형 루프 두 개를 좌우에 배치해 8자형 코일(figure-8 coil)을 만듭니다.

    왼쪽 날개(wing) 중심은 (offset_x - wing_offset, offset_y, h_standoff),
    오른쪽 날개 중심은 (offset_x + wing_offset, offset_y, h_standoff)입니다.

    인자(Parameters)
    ----------------
    wing_radius : float
        날개 반경(wing radius) a [m].
    wing_offset : float or None
        날개 중심의 x 방향 오프셋(offset) s [m]. None이면 wing_radius + 1 mm.
    h_standoff : float
        코일 높이 h [m].
    offset_x, offset_y : float
        코일 전체의 수평 이동량 [m] (수평 정렬 오차(lateral offset) 스윕에 사용).
    k_seg : int
        날개당 선분 개수 K.

    반환값(Returns)
    ---------------
    pts, dl : ndarray, shape (2*k_seg, 3)
        왼쪽 날개, 오른쪽 날개 순서로 쌓은 선분 중점과 현 벡터 [m].
    """
    if wing_offset is None:
        wing_offset = wing_radius + 1e-3
    # === TODO [Task 1.1b] ===
    # 중앙 접합부(junction, x = offset_x)에서 두 날개의 전류가 모두 +a_y 방향이 되도록
    # 각 날개의 direction_sign을 정한 뒤, np.vstack으로 두 날개를 합치십시오.
    raise NotImplementedError("TODO [Task 1.1b]: build_figure8_coil()을 구현하십시오.")


def compute_unit_A_and_B(
    pts: np.ndarray,
    dl: np.ndarray,
    obs_points: np.ndarray,
    n_turns: int = N_TURNS,
    batch_size: int = 2500,
) -> tuple[np.ndarray, np.ndarray]:
    """단위 전류(unit current, I = 1 A)당 벡터 퍼텐셜(vector potential)과 자기선속밀도를 계산합니다.

    인자(Parameters)
    ----------------
    pts, dl : ndarray, shape (K, 3)
        코일 선분 중점과 현 벡터 [m] (make_circular_loop / build_*_coil의 출력).
    obs_points : ndarray, shape (M, 3) 또는 (3,)
        관측점(observation point) 위치 [m].
    n_turns : int
        권선 수 N (모든 선분에 곱해짐).
    batch_size : int
        메모리 절약을 위해 한 번에 계산할 관측점 개수.

    반환값(Returns)
    ---------------
    A_tilde : ndarray, shape (M, 3)
        단위 전류당 벡터 퍼텐셜 Ã [Wb/(m·A)]. 실제 A = I · Ã.
    B_tilde : ndarray, shape (M, 3)
        단위 전류당 자기선속밀도(magnetic flux density) B̃ [T/A]. 실제 B = I · B̃.
    """
    obs = np.atleast_2d(np.asarray(obs_points, dtype=float))
    # === TODO [Task 1.1c] ===
    #   Ã(r) = (MU0 * n_turns / (4*pi)) * Σ_k dl_k / R_k
    #   B̃(r) = (MU0 * n_turns / (4*pi)) * Σ_k (dl_k × (r - r'_k)) / R_k^3,   R_k = |r - r'_k|
    # 팁: obs를 batch_size개씩 나누고 NumPy 브로드캐스팅(broadcasting)을 사용하면 빠릅니다.
    raise NotImplementedError("TODO [Task 1.1c]: compute_unit_A_and_B()를 구현하십시오.")


def analytical_circular_on_axis_Bz(
    z_arr: np.ndarray,
    radius: float = R_CIRC,
    h_standoff: float = H_NOMINAL,
    n_turns: int = N_TURNS,
    current: float = 1.0,
) -> np.ndarray:
    """중심 (0, 0, h_standoff)인 N회 권선 원형 루프의 중심축(on-axis) B_z(0, 0, z) 해석해를 계산합니다.

    인자(Parameters)
    ----------------
    z_arr : ndarray
        중심축 위의 z 좌표 [m].
    radius : float
        루프 반경 a [m].
    h_standoff : float
        루프 높이 h [m].
    n_turns : int
        권선 수 N.
    current : float
        코일 전류 I [A].

    반환값(Returns)
    ---------------
    Bz : ndarray, z_arr와 같은 shape
        중심축 자기선속밀도 B_z [T].
    """
    # === TODO [Task 1.1d] ===
    # 보고서에서 손으로 유도한 중심축 공식을 구현하십시오.
    raise NotImplementedError("TODO [Task 1.1d]: analytical_circular_on_axis_Bz()를 구현하십시오.")


def verify_curl_A(
    pts: np.ndarray,
    dl: np.ndarray,
    r0: np.ndarray,
    delta: float = 1e-5,
) -> tuple[np.ndarray, np.ndarray, float]:
    """점 r0에서 중앙 유한 차분(central finite difference)으로 ∇×Ã를 구해 B̃와 비교합니다.

    인자(Parameters)
    ----------------
    pts, dl : ndarray, shape (K, 3)
        코일 선분 중점과 현 벡터 [m].
    r0 : ndarray, shape (3,)
        검증할 점의 위치 [m].
    delta : float
        유한 차분 간격 δ [m].

    반환값(Returns)
    ---------------
    B_tilde : ndarray, shape (3,)
        비오-사바르(Biot–Savart)로 구한 B̃(r0) [T/A].
    curl_A : ndarray, shape (3,)
        유한 차분으로 구한 ∇×Ã(r0) [T/A].
    rel_error : float
        상대오차(relative error) |curl_A − B_tilde| / |B_tilde|.
    """
    # === TODO [Task 1.1e] ===
    # r0 ± delta (x, y, z 방향, 총 6점)에서 compute_unit_A_and_B로 Ã를 구하고,
    # ∂A_i/∂x_j ≈ (A_i(+) − A_i(−)) / (2·delta)로 회전(curl)의 세 성분을 조립하십시오.
    raise NotImplementedError("TODO [Task 1.1e]: verify_curl_A()를 구현하십시오.")


# ==============================================================================
# Task 1.2 — 전류 펄스, 인덕턴스, 구동 전압(driving voltage)
# ==============================================================================
def trapezoidal_pulse_current_and_didt(
    t_s: np.ndarray,
    i_peak: float = I_PEAK_NOMINAL,
    t_rise: float = T_RISE_S,
    t_hold: float = T_HOLD_S,
    t_fall: float = T_FALL_S,
) -> tuple[np.ndarray, np.ndarray]:
    """사다리꼴 전류 펄스(trapezoidal current pulse) I(t)와 시간 미분 dI/dt를 만듭니다.

    인자(Parameters)
    ----------------
    t_s : ndarray
        시간 배열 [s].
    i_peak : float
        첨두 전류(peak current) I_peak [A].
    t_rise, t_hold, t_fall : float
        상승(rise), 평탄 유지(hold), 하강(fall) 구간 길이 [s].

    반환값(Returns)
    ---------------
    i_t : ndarray, t_s와 같은 shape
        코일 전류 I(t) [A].
    di_dt : ndarray, t_s와 같은 shape
        전류의 시간 미분 dI/dt [A/s].
    """
    # === TODO [Task 1.2a] ===
    #   상승 [0, t_rise)                : I = (i_peak / t_rise)·t,   dI/dt = +i_peak / t_rise
    #   유지 [t_rise, t_rise + t_hold)  : I = i_peak,                dI/dt = 0
    #   하강 [t_rise + t_hold, t_end)   : I가 0까지 선형 감소,       dI/dt = −i_peak / t_fall
    #   펄스 밖                         : I = 0,                     dI/dt = 0
    raise NotImplementedError("TODO [Task 1.2a]: trapezoidal_pulse_current_and_didt()를 구현하십시오.")


def compute_inductance(
    pts1: np.ndarray,
    dl1: np.ndarray,
    pts2: np.ndarray | None = None,
    dl2: np.ndarray | None = None,
    n_turns: int = N_TURNS,
    wire_radius: float = WIRE_RADIUS,
) -> float:
    """쇄교 자속(flux linkage) 이중합으로 자기 인덕턴스(self-inductance) 또는 상호 인덕턴스(mutual inductance)를 계산합니다.

    인자(Parameters)
    ----------------
    pts1, dl1 : ndarray, shape (K1, 3)
        첫 번째 루프(또는 코일)의 선분 중점과 현 벡터 [m].
    pts2, dl2 : ndarray, shape (K2, 3) 또는 None
        두 번째 루프의 선분. None이면 pts1, dl1을 그대로 써서 자기 인덕턴스 L을 계산하고,
        주어지면 상호 인덕턴스 M_12를 계산합니다.
    n_turns : int
        권선 수 N (결과에 N^2이 곱해짐).
    wire_radius : float
        도선 반경 r_w [m] (정칙화(regularization) 거리 d_eff 계산에 사용).

    반환값(Returns)
    ---------------
    L_or_M : float
        인덕턴스 [H].
    """
    # === TODO [Task 1.2b] ===
    #   L 또는 M_12 = (MU0 * n_turns^2 / (4*pi)) * Σ_{i,j} (dl1_i · dl2_j) / R_ij
    #   R_ij = sqrt(|r1_i − r2_j|^2 + d_eff^2),  d_eff = wire_radius * exp(−0.25)  (세선 GMD 정칙화)
    # L_1 + L_2 + 2·M_12 = L_fig8가 정확히 성립하도록 상호 인덕턴스에도 같은 d_eff를 사용하십시오.
    raise NotImplementedError("TODO [Task 1.2b]: compute_inductance()를 구현하십시오.")


# ==============================================================================
# Task 1 검증 실행(이 파일을 직접 실행할 때)
# ==============================================================================
def main() -> None:
    circ_pts, circ_dl = build_circular_coil()
    fig8_pts, fig8_dl = build_figure8_coil()

    # --- Task 1.1: 중심축 B_z 해석해 vs. 수치해 ---
    z_axis = np.linspace(-60e-3, 0.0, 121)
    obs_axis = np.column_stack([np.zeros_like(z_axis), np.zeros_like(z_axis), z_axis])
    _, b_axis = compute_unit_A_and_B(circ_pts, circ_dl, obs_axis)
    bz_num = b_axis[:, 2] * I_PEAK_NOMINAL
    bz_ana = analytical_circular_on_axis_Bz(z_axis, current=I_PEAK_NOMINAL)
    max_rel_err = float(np.max(np.abs(bz_num - bz_ana) / np.abs(bz_ana)))
    print(f"[Task 1.1] On-axis Bz max rel error = {max_rel_err:.3e} (goal < 1e-4)")

    # --- Task 1.1: B = ∇×A 검증 ---
    _, _, err_c = verify_curl_A(circ_pts, circ_dl, R_OFF_TARGET)
    _, _, err_f = verify_curl_A(fig8_pts, fig8_dl, R_TARGET)
    print(f"[Task 1.1] curl(A) rel error: Circular @ (45,0,-20) mm = {err_c:.3e}, Figure-8 @ (0,0,-20) mm = {err_f:.3e}")

    # --- Task 1.2: 인덕턴스 분해와 구동 전압 ---
    i_t, _ = trapezoidal_pulse_current_and_didt(T_S)
    i_peak = float(np.max(np.abs(i_t)))

    wing1_pts, wing1_dl = make_circular_loop(-S_FIG8, 0.0, H_NOMINAL, R_FIG8, direction_sign=+1.0)
    wing2_pts, wing2_dl = make_circular_loop(+S_FIG8, 0.0, H_NOMINAL, R_FIG8, direction_sign=-1.0)
    l_circ = compute_inductance(circ_pts, circ_dl)
    l_w1 = compute_inductance(wing1_pts, wing1_dl)
    l_w2 = compute_inductance(wing2_pts, wing2_dl)
    m_12 = compute_inductance(wing1_pts, wing1_dl, wing2_pts, wing2_dl)
    l_fig8 = compute_inductance(fig8_pts, fig8_dl)

    print(f"[Task 1.2] L_circ = {l_circ*1e6:.3f} uH | V_coil = +/-{l_circ*DIDT_RISE*1e-3:.3f} kV | Wm = {0.5*l_circ*i_peak**2:.2f} J")
    print(f"           Figure-8: L1 = {l_w1*1e6:.3f} uH, L2 = {l_w2*1e6:.3f} uH, M12 = {m_12*1e6:+.3f} uH")
    print(f"           L1 + L2 + 2*M12 = {(l_w1 + l_w2 + 2*m_12)*1e6:.6f} uH vs direct L_fig8 = {l_fig8*1e6:.6f} uH")
    print(f"           V_coil = +/-{l_fig8*DIDT_RISE*1e-3:.3f} kV | Wm = {0.5*l_fig8*i_peak**2:.2f} J")


if __name__ == "__main__":
    main()
