#!/bin/bash
# Batch-generate Hirshfeld-I anion reference densities (Route 1, Psi4).
# Usage: batch_anions.sh <LOT> <outdir>   e.g. batch_anions.sh PBE rho_pbe
set -u
LOT="${1:-PBE}"
OUT="${2:-rho_${LOT,,}}"
GEN="$(dirname "$0")/gen_anion_rho.py"
mkdir -p "$OUT"
SUMMARY="$OUT/SUMMARY.txt"
: > "$SUMMARY"

# q=-1 for all main-group H-Kr (incl. noble gases, flagged as unphysical)
Q1="h he li be b c n o f ne na mg al si p s cl ar k ca ga ge as se br kr"
# q=-2 only where chemically conceivable (groups 13-16 p-block)
Q2="b c n o al si p s ga ge as se"

for el in $Q1; do
  python "$GEN" "$el" -1 --lot "$LOT" --outdir "$OUT" 2>>"$OUT/errors.log" | tee -a "$SUMMARY"
done
for el in $Q2; do
  python "$GEN" "$el" -2 --lot "$LOT" --outdir "$OUT" 2>>"$OUT/errors.log" | tee -a "$SUMMARY"
done

echo "=== done $LOT -> $OUT ===" | tee -a "$SUMMARY"
