# VibraStruct
## Structural Dynamics and Vibration Analysis

Version 1.1 — reviewed numerical implementation.

A reproducible Python project exploring natural frequencies, mode shapes, resonance, damping, transient vibration and stiffness sensitivity.

### Open the browser dashboard (Windows)
1. Extract the ZIP and open the inner `VibraStruct` folder containing `app.py`.
2. Click the File Explorer address bar, type `cmd`, and press Enter.
3. Run `python -m pip install -r requirements.txt`.
4. Run `python -m streamlit run app.py`.
5. Your browser should open automatically. Otherwise open http://localhost:8501 while the command window stays open.

For subsequent launches, double-click `OPEN_VIBRASTRUCT.bat` after dependencies are installed. If port 8501 is busy, use `python -m streamlit run app.py --server.port 8502` and open http://localhost:8502. Close the server with Ctrl+C. A local address works on the laptop running the server; it is not a public portfolio URL.

The dashboard supports single-mass and coupled three-mass models, adjustable parameters, harmonic and transient response, mode shapes and CSV/PNG/JSON downloads. Press **Update analysis** after editing the sidebar inputs.

### Publish with Streamlit Community Cloud
This repository keeps the project in the `VibraStruct` subfolder. In Streamlit Community Cloud, select branch `main` and `VibraStruct/app.py` as the main file. Dependencies are in `VibraStruct/requirements.txt`; dashboard styling is in `.streamlit/config.toml` at the repository root. Live dashboard: https://vibrastruct.streamlit.app/.

### Run the original reproducible study
1. Extract the ZIP, then open a terminal inside `VibraStruct`.
2. Create an environment: `python -m venv .venv`.
3. Activate it on Windows: `.venv\Scripts\activate`; on macOS/Linux: `source .venv/bin/activate`.
4. Install dependencies: `python -m pip install -r requirements.txt`.
5. Run the analysis: `python run_analysis.py` (includes automatic verification checks).
6. Verify it: `python -m unittest discover -s tests -v`.

Python 3.10 or newer is required. The supplied results can be read without installing anything. `requirements-tested.txt` records the exact scientific package versions used for the included results; use `python -m pip install -r requirements-tested.txt` to request those versions. This environment was tested on Python 3.12; other supported Python/package combinations have not all been tested. Internet is needed to install packages unless already available. The code itself runs offline. No paid software or confidential data is required.

### Contents
- `vibrastruct/core.py`: reusable numerical methods, SI units, input validation.
- `run_analysis.py`: complete reproducible study, six figures, data tables, numerical summary and generated report.
- `VibraStruct_Walkthrough.ipynb`: optional notebook; install Jupyter separately if desired.
- `results/`: completed plots, CSV data, summary and test log.
- `docs/TECHNICAL_REPORT.md`: methods, actual numerical findings, assumptions and limitations.
- `docs/LEARNING_GUIDE.md`: explanations and questions for discussing the project.
- `tests/`: analytical and physical verification checks.

### What the study investigates
The SDOF case compares damping ratios of 1%, 3% and 10% under a 100 N harmonic force. The coupled system represents an idealised support chain with three translating masses, one grounded end and one free end. It is analysed under separate free-vibration and harmonic-force cases. A parameter sweep scales all spring stiffnesses from 50% to 150% of baseline.

This is a lumped-parameter dynamics project. SolidFEM addresses a separate structural finite-element topic. VibraStruct does not claim to be a real machine model, a plant diagnosis or an experimentally validated design.

### Changing the study
Edit the clearly named baseline values in `run_analysis.py`, rerun the script and the tests, and inspect the regenerated tables. All generated results overwrite files in `results/` and the technical report. Keep a separate copy if retaining an alternative scenario. Keep masses in kg, springs in N/m, time in seconds and force in N. Angular frequency is rad/s; plots use Hz.

### Numerical checks and API details
The revised package includes 20 unit tests, plus script-level verification recorded in `results/verification.json`. The transient response is compared against an independent matrix-exponential solution; a refined integration measures sensitivity to solver settings. Harmonic response is checked against modal superposition. Actual numerical errors are in `results/summary.json` and the technical report.

The toolkit requires symmetric positive-definite mass and stiffness matrices. Free rigid-body modes, unconstrained structures and buckling/negative stiffness are outside its scope. `response(..., rtol=..., atol=..., max_step=...)` accepts explicit solver controls. Initial conditions apply at the first supplied time, even if it is not zero. `harmonic_response` returns complex displacement amplitudes for the exp(+i omega t) convention; its frequency input is rad/s. `receptance` explicitly rejects an exactly singular undamped SDOF resonance. Near-singular undamped models remain numerically sensitive; use meaningful damping and convergence checks when studying peaks. Modal sign is fixed deterministically for repeatable figures; it has no physical significance.

### Data and reproducibility
All input parameters are synthetic; no site-specific measurements are used. Results are reproducible from the supplied code. See the technical report for model assumptions, verification and limitations. Software is supplied under the MIT license.
