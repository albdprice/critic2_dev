#!/bin/bash
# Batch-generate Route 2 (confined ld1.x) Hirshfeld-I anion references for
# the whole periodic table. Same numerical scalar-relativistic PBE method as
# critic2's shipped neutral/cation dat/wfc set. Standardized confinement box
# rmax = alpha * R99(neutral), with step-down fallback for unbound cases.
# Usage: batch_anions_ld1.sh <outdir> [alpha]
set -u
set -o pipefail   # so a generator failure is seen through the tee pipe
# Note: stock ld1.x caps the nuclear charge at Z=103; Rf-Og (Z=104-118)
# require an ld1.x patched to raise the max Z (see gen.m header), and are
# skipped here. Superheavy anions are not chemically meaningful anyway.
OUT="${1:-rho_ld1_all}"
ALPHA="${2:-3.6}"
GEN="$(dirname "$0")/gen_anion_rho_ld1.py"
mkdir -p "$OUT"
SUMMARY="$OUT/SUMMARY.txt"; : > "$SUMMARY"

# q = -1 for every element H..Og
Q1="h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v cr mn fe co ni \
cu zn ga ge as se br kr rb sr y zr nb mo tc ru rh pd ag cd in sn sb te i xe cs \
ba la ce pr nd pm sm eu gd tb dy ho er tm yb lu hf ta w re os ir pt au hg tl pb \
bi po at rn fr ra ac th pa u np pu am cm bk cf es fm md no lr rf db sg bh hs mt \
ds rg cn nh fl mc lv ts og"
# q = -2 for the p-block (groups 13-16), where doubly-charged anions are conceivable
Q2="b c n o al si p s ga ge as se in sn sb te tl pb bi po"

for el in $Q1; do
  python3 "$GEN" "$el" -1 --alpha "$ALPHA" --outdir "$OUT" 2>>"$OUT/errors.log" \
    | tee -a "$SUMMARY" || echo "FAILED: $el q-1" | tee -a "$SUMMARY"
done
for el in $Q2; do
  python3 "$GEN" "$el" -2 --alpha "$ALPHA" --outdir "$OUT" 2>>"$OUT/errors.log" \
    | tee -a "$SUMMARY" || echo "FAILED: $el q-2" | tee -a "$SUMMARY"
done
rm -f "$OUT"/ld1.wfc "$OUT"/*.ld1in "$OUT"/ld1.* 2>/dev/null
echo "=== done -> $OUT (alpha=$ALPHA) ===" | tee -a "$SUMMARY"
