# UniPEFF test fixture

A self-contained regression test for the UniProtKB-XML → PEFF conversion. There is no
unit-test framework; `run.sh` diffs the tool's output against committed golden files.

## Run

```bash
./run.sh
# dotnet not on PATH?  ->  DOTNET=$HOME/.dotnet/dotnet ./run.sh
```

It builds the project, converts `comprehensive.xml` in both modes, diffs against
`expected_A.peff` (default) and `expected_B.peff` (`-AnnotationIdentifiers`), and runs a few
PEFF-conformance assertions (mandatory `DbVersion` header present, ASCII-only output). Exit 0
= match. The tool reads `ptmlist.txt` from the working directory, so `run.sh` runs from this
folder. If a code change *intentionally* alters output, regenerate the goldens with the two
`dotnet … -out expected_*.peff` commands.

## What `comprehensive.xml` covers

Entry **P12345** (Swiss-Prot) is a deliberate "kitchen sink":

- **Entry fields:** organism (`\NcbiTaxId`/`\TaxName`, with a `<name type="common">` before the
  taxonomy reference to exercise the bounded organism scan), sequence/entry versions
  (`\SV`/`\EV`), protein existence (`\PE`); two accessions (the first is the primary id, the
  second becomes `\AltAC`); mnemonic (`\ID`); a fullName mixing paired and unpaired parens plus
  `|` and `\` (escaping: only `|`, `\` and *unpaired* parens are backslash-escaped; balanced
  parens are left intact; `\PName` is a bare scalar value, **not** a `(parenthesized)` list); and
  a gene whose **synonym precedes the primary name** — a bounded `<gene>` subtree scan still picks
  the primary `\GName`, independent of child order.
- **Molecular processing:** all five types → their PEFF CVs — initiator methionine (single
  position), signal peptide, transit peptide, propeptide, chain; plus an unknown-position
  processing feature that MUST be omitted (`\Processed` positions count from 1; "?" is
  ModRes-only).
- **Modifications:** PSI-MOD (`\ModResPsi`), Unimod-only (`\ModResUnimod`), generic
  no-accession (`\ModRes`), with names resolved from the bundled `psi-mod.obo`/`unimod.obo`
  (PEFF requires the OBO `name:`, e.g. `O-phospho-L-serine`, not the UniProt synonym); the
  `(Microbial infection)` and `;`-suffix description strips; a description absent from
  `ptmlist.txt` (warns, produces no output, does not crash); and an **unknown-position**
  modified residue and glycosylation site (`<position status="unknown"/>`) emitted as `?`
  (the writer renders position 0 as `?`, legal only inside ModRes) rather than crashing.
- **Glycosylation / lipidation / cross-links → modifications:** a glycosylation site with an
  embedded paired-paren description (`N-linked (GlcNAc...)`, left unescaped) and one with no
  description (name falls back to the feature type); a lipidation matching `ptmlist.txt`
  (→ `\ModResPsi`) and one that does not (→ generic `\ModRes`); a single-position cross-link
  and an intra-chain begin/end cross-link (→ two `\ModRes` residues).
- **Disulfides:** two-ended (45–80), unknown begin (→ position 0 → `?`, …90) and single-ended
  cross-chain (95). Verifies the half-cystine merge feeds `\ModResPsi` position-sorted and,
  in Option B, that `\DisulfideBond` references the correct half-cystine ids — including the
  unknown-position one by **id, not position** — while the unpaired cross-chain cystine is
  `\ModResPsi`-only.
- **Variants:** simple (25), complex replacement (30–31), single-residue deletion (33), range
  deletion (36–38); a single-residue substitution UniProt expressed as a range (42–42|M) that
  is **demoted** to `\VariantSimple`; a single-position deletion expressed with `<original>` but
  **no `<variation>`** (→ `\VariantComplex (44|44|)`, not a crash); a **multi-residue insertion**
  at one position (`A→APT` → `\VariantComplex (46|46|APT)`, never a truncated `\VariantSimple
  (46|A)`); and out-of-bounds variants (a simple at 250, a complex 200–201) that are **omitted
  with a warning** (positions must lie within the length-100 sequence).
- **Skipped entries:** a trailing `<entry>` with **no accession**, and another with an accession
  but **no sequence**, are both skipped (a PEFF entry is `>Prefix:DbUniqueId` + a sequence block),
  so only two entries are written and `NumberOfEntries` counts only those.
- **Still-unhandled feature types** (`sequence conflict`, `non-standard amino acid`) are parsed
  without error and produce no output.

`compact.xml` is the same kind of entry serialized with **no pretty-print whitespace** between
elements; `run.sh` asserts the name/depth-driven parser extracts identical annotations from it
(the previous fixed-`Read()`-count parser was whitespace-layout-dependent and would mis-parse it).

Entry **Q67890** (TrEMBL) is annotation-free but carries metadata (organism, versions,
`\PE=4` predicted, a second accession → `\AltAC`): it checks multi-entry handling, the
annotation-free `>` line, **per-dataset database blocks** (being TrEMBL it is emitted under
its own `# Prefix=tr` / `>tr:` block, separate from the Swiss-Prot block — and a `-prefix`
override instead collapses everything into one block), and non-ASCII text in the protein name
(sanitized to ASCII, as PEFF requires).

## Integration test (optional, real data)

The fixture above is a controlled unit test. For an end-to-end check at scale, convert the
whole reviewed human proteome. The tool reads `ptmlist.txt`, `psi-mod.obo` and `unimod.obo`
from the working directory; the OBO files supply the canonical modification names PEFF
requires — **without them, mod names fall back to the UniProt `ptmlist.txt` IDs** (a warning
is printed and the output is not strictly OBO-named). Download all three plus the proteome:

```bash
curl -O 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt'
curl -fsSL 'https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/master/PSI-MOD.obo' -o psi-mod.obo
curl -O 'http://www.unimod.org/obo/unimod.obo'
curl -G 'https://rest.uniprot.org/uniprotkb/stream' \
  --data-urlencode 'query=organism_id:9606 AND reviewed:true' \
  --data-urlencode 'format=xml' --data-urlencode 'compressed=true' \
  -o human_sprot.xml.gz
dotnet "$ROOT/bin/Debug/net8.0/UniPEFF.dll" -in human_sprot.xml.gz -out human.peff                       # default
dotnet "$ROOT/bin/Debug/net8.0/UniPEFF.dll" -in human_sprot.xml.gz -out human_B.peff -AnnotationIdentifiers
```

Last verified run: 20,431 entries (~169 MB gz in → ~19 MB PEFF out, ~12 s), no crash and no
missing-PTM warnings. Spot check — insulin `P01308` emits its six half-cystines as `\ModResPsi`
and, with `-AnnotationIdentifiers`, the correct disulfide connectivity `31-96, 43-109, 95-100`.
