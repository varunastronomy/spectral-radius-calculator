<p align="center"><img src="assets/banner.svg" alt="Spectral Radius Calculator cover" width="100%"></p>

# Spectral Radius Calculator

[![Quality checks](https://github.com/varunastronomy/spectral-radius-calculator/actions/workflows/quality-checks.yml/badge.svg)](https://github.com/varunastronomy/spectral-radius-calculator/actions/workflows/quality-checks.yml)

Convert diskbb and bbodyrad normalizations into apparent and colour-corrected radius diagnostics with propagated uncertainties.

> **Scientific boundary:** These radii depend on distance, inclination, hardening factors, boundary corrections, geometry, and continuum decomposition; they are not direct geometric measurements.

## Copyright, permission, and mandatory citation

> **COPYRIGHT © 2026 M. VARUN, CHRIST (DEEMED TO BE UNIVERSITY), BENGALURU, INDIA. ALL RIGHTS RESERVED.**
>
> **NO PERMISSION IS GRANTED TO USE, COPY, MODIFY, REDISTRIBUTE, REPUBLISH, OR INCORPORATE THIS SOFTWARE INTO ANOTHER PROJECT WITHOUT PRIOR WRITTEN PERMISSION FROM THE AUTHOR. APPROVED USE REQUIRES CLEAR ATTRIBUTION AND CITATION.**

## Installation and use

```bash
python -m pip install -r requirements.txt
python spectral_radius_calculator.py --help
```

Some interactive scripts do not implement `--help`; run them from a prepared analysis directory and follow their prompts. Read the source-level defaults before execution. Mission-specific programs additionally require the relevant mission software, calibration, network access, and data layout.

## Provenance and validation

The scientific source is copied from `radius_calculation.py` in the private research archive. The archived original remains unchanged. Repository checks validate Python syntax only. No real observation, mission pipeline, fit, download, plot suite, or scientific result was executed or validated for publication. Outputs must be inspected before scientific use.

No open-source license is granted. Citation metadata are in [`CITATION.cff`](CITATION.cff).

## Author

**M. Varun**<br>X-ray astronomy researcher and scientific Python developer<br>CHRIST (Deemed to be University), Bengaluru, India<br>[GitHub: varunastronomy](https://github.com/varunastronomy)
