#!/usr/bin/env bash
# Build nucleotides.npy, the data used by examples/trinucleotides.R.
#
# For each trinucleotide, downloads the 0.5 A conformer library
# (github.com/sjdv1982/nucleotide-library) and its PDB template
# (github.com/sjdv1982/nucleotide-fragment-templates), runs
# npy2pdb2angles.py (next to this script) to get the dihedral angles, keeps
# the middle residue of every conformer and saves all trinucleotides as a
# dict {trinucleotide: DataFrame} in a pickled .npy.
#
# Usage (from examples/data/):
#   ./build_nucleotides.sh [output.npy]     # default: ./nucleotides.npy
#
# Requires: curl, python3 with numpy and pandas.
# Environment: PYTHON (default python3).
#
# This file is AI-generated with review.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT="${1:-$SCRIPT_DIR/nucleotides.npy}"
[[ "$OUTPUT" = /* ]] || OUTPUT="$PWD/$OUTPUT"   # absolute: we cd to a temp dir
PYTHON="${PYTHON:-python3}"

LIBRARY_URL="https://raw.githubusercontent.com/sjdv1982/nucleotide-library/main/library"
TEMPLATE_URL="https://raw.githubusercontent.com/sjdv1982/nucleotide-fragment-templates/main"

TRINUCS=(AAA AAC CAC CAA CCC CCA ACC ACA)

if ! "$PYTHON" -c "import numpy, pandas" 2>/dev/null; then
    echo "Error: $PYTHON needs numpy and pandas (pip install numpy pandas)" >&2
    exit 1
fi

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

### Download ##################################################################

echo "Working in $WORKDIR"
for t in "${TRINUCS[@]}"; do
    echo "Downloading $t library and template..."
    curl -fsSL -O "$LIBRARY_URL/trinuc-$t-0.5.npy"
    curl -fsSL -O "$TEMPLATE_URL/$t.pdb"
done

### Angles ####################################################################

echo "Computing angles..."
for t in "${TRINUCS[@]}"; do
    "$PYTHON" "$SCRIPT_DIR/npy2pdb2angles.py" "trinuc-$t-0.5.npy" "$t.pdb" -o "$t"
    echo "  $t done"
done

### Assemble ##################################################################

echo "Writing $OUTPUT"
mkdir -p "$(dirname "$OUTPUT")"
"$PYTHON" - "$OUTPUT" "${TRINUCS[@]}" <<'EOF'
import sys
import numpy as np
import pandas as pd

output, trinucs = sys.argv[1], sys.argv[2:]
angles = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta",
          "nu0", "nu1", "nu2", "nu3", "nu4", "eta", "theta", "chi"]

data = {}
for t in trinucs:
    # keep_default_na=False so a blank chain id stays ' ' instead of NaN
    df = pd.read_csv(f"angles_{t}.csv", keep_default_na=False, na_values=[""])
    middle = df[df["res_id"] == 2]   # middle residue: all angles are defined
    data[t] = middle[["chain"] + angles]
    print(f"  {t}: {len(middle)} conformers")

np.save(output, data, allow_pickle=True)
EOF

echo "Done."
