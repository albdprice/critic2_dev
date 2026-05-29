#!/bin/bash
# Batch-generate Route 2 (confined ld1.x) Hirshfeld-I anion references.
# Same method as critic2's shipped neutral/cation dat/wfc set, so all
# charge states are consistent. Usage: batch_anions_ld1.sh <outdir>
set -u
OUT="${1:-rho_ld1}"
GEN="$(dirname "$0")/gen_anion_rho_ld1.py"
mkdir -p "$OUT"
SUMMARY="$OUT/SUMMARY.txt"; : > "$SUMMARY"

Q1="h he li be b c n o f ne na mg al si p s cl ar k ca ga ge as se br kr"
Q2="b c n o al si p s ga ge as se"

for el in $Q1; do
  python3 "$GEN" "$el" -1 --outdir "$OUT" 2>>"$OUT/errors.log" | tee -a "$SUMMARY" \
    || echo "FAILED: $el q-1" | tee -a "$SUMMARY"
done
for el in $Q2; do
  python3 "$GEN" "$el" -2 --outdir "$OUT" 2>>"$OUT/errors.log" | tee -a "$SUMMARY" \
    || echo "FAILED: $el q-2" | tee -a "$SUMMARY"
done
# tidy ld1.x scratch
rm -f "$OUT"/ld1.wfc "$OUT"/*.ld1in "$OUT"/ld1.* 2>/dev/null
echo "=== done -> $OUT ===" | tee -a "$SUMMARY"
