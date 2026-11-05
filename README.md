# VibraStruct

Interactive Python dashboard for structural dynamics and vibration analysis.

**[Open the live dashboard](https://vibrastruct.streamlit.app/)**

Explore single-mass and coupled three-mass systems: natural frequencies, mode shapes, harmonic response, free vibration, damping and CSV/PNG/JSON downloads. All example inputs are synthetic.

## Run locally

Requires Python 3.10 or newer. From this repository's root:

```bash
cd VibraStruct
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Open the local URL shown in the terminal. Keep the terminal running while using the dashboard.

## Reproduce and verify

From the `VibraStruct` project directory:

```bash
python run_analysis.py
python -m unittest discover -s tests -v
```

The numerical suite contains 20 tests. The reproducible study additionally checks eigen residuals, modal orthogonality, energy behaviour, independent transient accuracy and harmonic modal agreement.

## Project files

| Path | Contents |
|---|---|
| [VibraStruct/app.py](VibraStruct/app.py) | Streamlit dashboard |
| [VibraStruct/vibrastruct/](VibraStruct/vibrastruct/) | Numerical toolkit: `__init__.py` and `core.py` |
| [VibraStruct/docs/TECHNICAL_REPORT.md](VibraStruct/docs/TECHNICAL_REPORT.md) | Methods, results, verification and limitations |
| [VibraStruct/results/](VibraStruct/results/) | Generated figures and data |
| [VibraStruct/tests/](VibraStruct/tests/) | Analytical and numerical checks |
| [VibraStruct/README.md](VibraStruct/README.md) | Full installation and usage guide |

## Deployment

Streamlit entry point: `VibraStruct/app.py`, branch `main`. Dependencies are in `VibraStruct/requirements.txt`. The repository-root `.streamlit/config.toml` defines dashboard styling.

## Scope and license

Linear lumped-parameter models with synthetic inputs. Numerical verification does not constitute experimental validation. The project does not determine safe operating conditions for real machinery. See the technical report for assumptions and limitations.

Licensed under the [MIT license](LICENSE).
