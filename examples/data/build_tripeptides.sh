#!/usr/bin/env bash
# Download tripeptides/, the data used by examples/tripeptides.R.
#
# Downloads the tripeptide angle archive from moma.laas.fr (~580 MB) and
# extracts only its data_coil/ folder (one <AA>_angles.RData per central amino
# acid) into tripeptides/, replacing any existing copy.
#
# Usage (from examples/data/):
#   ./build_tripeptides.sh
#
# Requires: curl, tar.
#
# This file is AI-generated with review.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT="$SCRIPT_DIR/tripeptides"

URL="https://moma.laas.fr/static/data/tripeptide_angles_data.tar"
MEMBER="tripeptide_angles_data/data_coil"
N_FILES=20

# Temp dir next to the output so the final mv is a rename on the same disk
WORKDIR="$(mktemp -d "$SCRIPT_DIR/.tripeptides.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

echo "Downloading $URL (~580 MB) and extracting data_coil/..."
curl -fSL "$URL" | tar -x -C "$WORKDIR" --strip-components=2 "$MEMBER"

n=$(find "$WORKDIR" -name '*_angles.RData' | wc -l)
if [ "$n" -ne "$N_FILES" ]; then
    echo "Error: expected $N_FILES .RData files, got $n" >&2
    exit 1
fi

rm -rf "$OUTPUT"
mv "$WORKDIR" "$OUTPUT"
chmod 755 "$OUTPUT"   # mktemp creates it as 700
echo "Saved $n files to $OUTPUT"
