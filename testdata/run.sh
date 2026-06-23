#!/usr/bin/env bash
# Regression test: convert comprehensive.xml and diff UniPEFF's PEFF output against the
# committed golden files, for both default and -AnnotationIdentifiers.
# `dotnet` may not be on PATH; override with: DOTNET=/path/to/dotnet ./run.sh
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
DOTNET="${DOTNET:-dotnet}"

"$DOTNET" build "$ROOT/UniPEFF.csproj" -nologo -clp:NoSummary >/dev/null
DLL="$ROOT/bin/Debug/net8.0/UniPEFF.dll"

cd "$HERE"   # the tool reads ptmlist.txt from the current working directory
trap 'rm -f actual_A.peff actual_B.peff actual_compact.peff' EXIT

"$DOTNET" "$DLL" -in comprehensive.xml -out actual_A.peff >/dev/null
"$DOTNET" "$DLL" -in comprehensive.xml -out actual_B.peff -AnnotationIdentifiers >/dev/null

diff -u expected_A.peff actual_A.peff
diff -u expected_B.peff actual_B.peff
echo "PASS: PEFF output matches goldens."

# Whitespace independence: compact.xml has no pretty-print whitespace between elements. The
# name/depth-driven parser MUST still extract the same annotations (the old fixed-Read()-count
# parser would have crashed or mis-parsed this).
"$DOTNET" "$DLL" -in compact.xml -out actual_compact.peff >/dev/null
for want in '>sp:P00003' '\GName=CGENE' '\ModResPsi=(3|MOD:00046|O-phospho-L-serine)' '\VariantSimple=(5|V)' '\Processed=(1|4|PEFF:0001021|signal peptide)'; do
  grep -qF -- "$want" actual_compact.peff || { echo "FAIL: compact.xml output missing: $want"; exit 1; }
done
echo "PASS: compact (whitespace-stripped) XML parses identically."

# Spec-conformance assertions, independent of the goldens:
#  - the DB description block MUST contain DbVersion (PEFF 1.0 section 3.3.2)
#  - PEFF permits only ASCII characters
for f in expected_A.peff expected_B.peff; do
  grep -q '^# DbVersion=' "$f" || { echo "FAIL: $f lacks the mandatory '# DbVersion=' header"; exit 1; }
  if LC_ALL=C grep -q '[^ -~]' "$f"; then echo "FAIL: $f contains non-ASCII or control characters"; exit 1; fi
  # multi-DB blocks: the TrEMBL entry MUST sit under its own >tr: block, and there must be
  # three '# //' separators (file-description block + sp block + tr block).
  grep -q '^>tr:Q67890' "$f" || { echo "FAIL: $f: TrEMBL entry is not under a >tr: block"; exit 1; }
  test "$(grep -c '^# //' "$f")" -eq 3 || { echo "FAIL: $f: expected 3 '# //' separators (file + sp + tr blocks)"; exit 1; }
  # isoforms: the alternative-products comment + splice-variant feature MUST be ignored (isoform
  # expansion is not implemented), so no isoform-id (e.g. P12345-2) descriptor line may appear.
  if grep -qE '^>[a-z]+:[A-Za-z0-9]+-[0-9]+ ' "$f"; then echo "FAIL: $f emits an isoform-id descriptor line"; exit 1; fi
done
echo "PASS: conformance checks (DbVersion, ASCII, per-dataset blocks, isoforms ignored)."
