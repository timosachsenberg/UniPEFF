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
