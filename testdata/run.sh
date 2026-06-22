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
trap 'rm -f actual_A.peff actual_B.peff' EXIT

"$DOTNET" "$DLL" -in comprehensive.xml -out actual_A.peff >/dev/null
"$DOTNET" "$DLL" -in comprehensive.xml -out actual_B.peff -AnnotationIdentifiers >/dev/null

diff -u expected_A.peff actual_A.peff
diff -u expected_B.peff actual_B.peff
echo "PASS: PEFF output matches goldens."

# Spec-conformance assertions, independent of the goldens:
#  - the DB description block MUST contain DbVersion (PEFF 1.0 section 3.3.2)
#  - PEFF permits only ASCII characters
for f in expected_A.peff expected_B.peff; do
  grep -q '^# DbVersion=' "$f" || { echo "FAIL: $f lacks the mandatory '# DbVersion=' header"; exit 1; }
  if LC_ALL=C grep -q '[^ -~]' "$f"; then echo "FAIL: $f contains non-ASCII or control characters"; exit 1; fi
done
echo "PASS: conformance checks (mandatory DbVersion present, ASCII-only)."
