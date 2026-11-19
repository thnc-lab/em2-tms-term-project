"""EEC2204 TMS 텀 프로젝트 공통 설정(shared configuration) — 제공 파일, 수정 금지.

모든 Task 파일(task1_coil_field.py, task2_tissue_induction.py, task3_neuron.py)은
물리 상수(physical constant), 코일·펄스 매개변수(parameter), 4층 두부 모델(4-layer head model)을
이 파일에서 불러 씁니다. 길이는 모두 미터[m], 시간은 초[s] 단위입니다.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

_mpl_dir = Path(tempfile.gettempdir()) / "em2_tms_mpl_cache"
_mpl_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_dir))

import numpy as np

# ==============================================================================
# 물리 상수 및 시스템 매개변수(physical constants & system parameters)
# ==============================================================================
MU0 = 4.0 * np.pi * 1e-7          # 자유공간 투자율(free-space permeability) [H/m]

# 코일 형상(coil geometry)
N_TURNS = 10                      # 코일(또는 날개)당 권선 수(number of turns) N
WIRE_RADIUS = 1.5e-3              # 도선 반경(wire radius) r_w [m]
K_SEG = 180                       # 원형 루프당 선분 분할 수(number of segments) K
R_CIRC = 45e-3                    # 원형 코일 반경 a_circ [m]
R_FIG8 = 35e-3                    # 8자형 코일 날개 반경(wing radius) a_fig8 [m]
S_FIG8 = 36e-3                    # 8자형 날개 중심의 x 방향 오프셋(offset) s [m]
H_NOMINAL = 5e-3                  # 명목 코일–두피 이격 거리(standoff distance) h [m]

# 사다리꼴 전류 펄스(trapezoidal current pulse), 총 길이 0.18 ms
I_PEAK_NOMINAL = 6000.0           # 명목 첨두 전류(peak current) I_peak [A]
T_RISE_S = 0.06e-3                # 상승 구간(rise) 길이 [s] -> dI/dt = +100 A/us
T_HOLD_S = 0.06e-3                # 평탄 유지 구간(hold) 길이 [s] -> dI/dt = 0
T_FALL_S = 0.06e-3                # 하강 구간(fall) 길이 [s] -> dI/dt = -100 A/us
DIDT_RISE = I_PEAK_NOMINAL / T_RISE_S  # 상승 구간의 dI/dt [A/s]

# 공통 시간 격자(time grid): 0~6 ms, 2 us 간격 (뉴런 솔버(solver)와 동일)
T_MS = np.linspace(0.0, 6.0, 3001)  # [ms]
T_S = T_MS * 1e-3                   # [s]

# 4층 평면 두부 모델(4-layer planar head model): 두피 표면 z = 0, 아래로 갈수록 z < 0
HEAD_LAYERS = [
    {"name": "Skin (Scalp)", "z_top": 0.0,     "z_bot": -5e-3,  "sigma": 0.33,  "color": "#f4d3b8"},
    {"name": "Skull",        "z_top": -5e-3,   "z_bot": -12e-3, "sigma": 0.010, "color": "#e5e5e5"},
    {"name": "CSF",          "z_top": -12e-3,  "z_bot": -15e-3, "sigma": 1.79,  "color": "#cce8f4"},
    {"name": "Brain",        "z_top": -15e-3,  "z_bot": -65e-3, "sigma": 0.33,  "color": "#f2c6c6"},
]

Z_TARGET = -20e-3                           # 피질 평가 평면(cortical evaluation plane) z [m]
R_TARGET = np.array([0.0, 0.0, Z_TARGET])   # 표적 뉴런(target neuron) 위치 [m]
R_OFF_TARGET = np.array([40e-3, 0.0, Z_TARGET])  # 비표적 뉴런(off-target neuron) 위치 [m]
U_NEURON = np.array([0.0, -1.0, 0.0])       # 축삭 방향 단위 벡터(axon unit vector), -a_y

E_GOAL_TRADEOFF = 150.0           # 트레이드오프(trade-off) 비교용 공통 목표 전기장 E_goal [V/m]

FIG_DIR = Path(__file__).resolve().parent / "figures"   # 그림 저장 폴더
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ==============================================================================
# 제공 보조 함수(provided helper)
# ==============================================================================
def get_layer_properties(z: float) -> tuple[str, float]:
    """깊이 z에 해당하는 조직 층(tissue layer)의 이름과 전기전도도(conductivity)를 반환합니다.

    인자(Parameters)
    ----------------
    z : float
        깊이 좌표 [m] (z <= 0).

    반환값(Returns)
    ---------------
    name : str
        층 이름 ("Skin (Scalp)", "Skull", "CSF", "Brain").
    sigma : float
        전기전도도(electrical conductivity) σ [S/m].
    """
    for layer in HEAD_LAYERS:
        if z >= layer["z_bot"] or layer["name"] == "Brain":
            return layer["name"], float(layer["sigma"])
    return HEAD_LAYERS[-1]["name"], float(HEAD_LAYERS[-1]["sigma"])
