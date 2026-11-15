#!/usr/bin/env bash
# Regenerate everything from Supplementary Data 1, in order. Run from the repo root
# with the project's environment active (python -m pip install -r requirements.txt).
#   bash scripts/run_all.sh             conventions, fit pass, QC, validation, tests
#   bash scripts/run_all.sh --with-sim  also the synthetic validation (~12 min)
# Override the interpreter with PYTHON=/path/to/python.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="./.venv/bin/python"

"$PY" --version
"$PY" -c "import hashlib, yaml; c = yaml.safe_load(open('config.yaml'))['paper']; \
h = hashlib.sha256(open(c['supplement_file'], 'rb').read()).hexdigest(); \
assert h == c['supplement_sha256'], f'supplement sha256 {h} != config'; print('supplement sha256 ok')"

echo "== conventions (src/validate.py)";      "$PY" src/validate.py
echo "== fit pass (src/fit_spectra.py)";      "$PY" src/fit_spectra.py
echo "== QC (src/qc.py)";                     "$PY" src/qc.py
echo "== validation (src/validation.py)";     "$PY" src/validation.py
echo "== Fig 1e diagnostic (src/fig1e_diagnostic.py)"; "$PY" src/fig1e_diagnostic.py
if [[ "${1:-}" == "--with-sim" ]]; then
  echo "== synthetic validation (src/simulate.py)"; "$PY" src/simulate.py
fi
echo "== tests";                              "$PY" -m pytest -q
