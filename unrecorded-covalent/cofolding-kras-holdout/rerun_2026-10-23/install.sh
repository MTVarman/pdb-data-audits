#!/bin/bash
set -e
cd /content
pip install -q gemmi rdkit requests numpy 2>&1 | tail -1
echo "--- protenix venv"
python3 -m venv --system-site-packages /content/px
/content/px/bin/pip install -q --upgrade pip 2>&1 | tail -1
/content/px/bin/pip install -q --upgrade protenix --index-url https://pypi.org/simple 2>&1 | tail -2
echo "--- boltz venv"
python3 -m venv --system-site-packages /content/bz
/content/bz/bin/pip install -q --upgrade pip 2>&1 | tail -1
/content/bz/bin/pip install -q -U boltz 2>&1 | tail -2
echo "--- check"
LAYERNORM_TYPE=torch /content/px/bin/protenix --help 2>&1 | head -3
/content/bz/bin/boltz predict --help 2>&1 | head -2
echo INSTALL_OK
