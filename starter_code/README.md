# EEC2204 전자기학 2 — 텀 프로젝트 시작 코드(starter code)

## 파일 구성

| 파일 | 내용 |
|---|---|
| `tms_config.py` | 물리 상수(physical constant), 코일·펄스 매개변수(parameter), 4층 두부 모델(4-layer head model), `get_layer_properties` (**수정 금지**) |
| `neuron_model.py` | 호지킨-헉슬리(Hodgkin–Huxley) 뉴런 시뮬레이터(simulator) `simulate_hh_neuron`, `find_electric_field_threshold` (**수정 금지**) |
| `task1_coil_field.py` | Task 1 — 코일 자기장과 인덕턴스 |
| `task2_tissue_induction.py` | Task 2 — 조직 레이어에서 유도 전기장 (Task 1 함수 사용) |
| `task3_neuron.py` | Task 3 — 대뇌 피질 뉴런 활성화 (Task 1 함수 사용) |
| `tms_simulation_starter.ipynb` | `task1~3`의 함수를 불러와(import) 검증을 실행하고 보고서용 그림을 그리는 주피터 노트북(Jupyter Notebook) (선택 사용) |

과제 설명은 저장소 루트의 `README.md`(과제 명세서)에, 각 함수의 인자(argument)와 입출력(input/output) 형식은 함수의 docstring에 있음. `# === TODO [Task X.Y] ===` 블록을 구현한 뒤 실행하면 해당 Task의 검증 값이 출력됨.

## 실행 환경

- Python 3.10 이상, `numpy`, `matplotlib` (`scipy` 불필요)
- 노트북을 사용할 경우 JupyterLab, VS Code, Google Colab 중 하나

먼저 환경을 확인함:

```bash
python neuron_model.py      # 정상이면 "HH self-test | E_th = 96.35 V/m" 출력
```

## 작성 방법

함수 구현은 항상 `task1_coil_field.py`, `task2_tissue_induction.py`, `task3_neuron.py`의 `# === TODO [Task X.Y] ===` 블록에서 함. Task 2와 Task 3은 `task1_coil_field.py`의 함수를 불러 쓰므로 Task 1을 먼저 완성함. 검증 실행은 다음 두 방법 중 편한 것을 사용함.

**(1) 터미널에서 `.py` 파일 직접 실행**

```bash
python task1_coil_field.py        # Task 1.1, 1.2
python task2_tissue_induction.py  # Task 2.2, 2.3 (Task 2.1은 보고서만)
python task3_neuron.py            # Task 3.1
```

**(2) 주피터 노트북(`tms_simulation_starter.ipynb`)에서 실행**

- 첫 코드 셀이 `task1~3_*.py`의 함수를 불러옴(import). 노트북에는 함수 구현이 없으며, 검증 셀과 보고서용 그림 셀만 있음.
- `%autoreload 2`가 켜져 있으므로 `.py` 파일을 수정·저장하면 커널을 재시작하지 않아도 다음 셀 실행부터 반영됨.
- 노트북은 `starter_code/` 폴더에서 열어야 함. Google Colab에서는 `tms_config.py`, `neuron_model.py`, `task1~3_*.py`를 함께 업로드함.
- 마지막의 "보고서용 그림" 셀부터 보고서에 필요한 그림을 자유롭게 작성함.

## 제출

- 완성한 `task1~3_*.py` 파일 3개 (필수)
- 노트북을 사용했다면 `tms_simulation_starter.ipynb` (선택)
- 보고서(PDF)와 그림 파일