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

- **Entry fields:** two accessions (only the first is kept), mnemonic, a fullName mixing
  paired and unpaired parens plus `|` and `\` (escaping: only `|`, `\` and *unpaired* parens
  are backslash-escaped; balanced parens are left intact), and a gene with a synonym (only
  the primary name is kept).
- **Molecular processing:** all five types → their PEFF CVs — initiator methionine (single
  position), signal peptide, transit peptide, propeptide, chain; plus an unknown-position
  processing feature that MUST be omitted (`\Processed` positions count from 1; "?" is
  ModRes-only).
- **Modifications:** PSI-MOD (`\ModResPsi`), Unimod-only (`\ModResUnimod`), generic
  no-accession (`\ModRes`), with names resolved from the bundled `psi-mod.obo`/`unimod.obo`
  (PEFF requires the OBO `name:`, e.g. `O-phospho-L-serine`, not the UniProt synonym); the
  `(Microbial infection)` and `;`-suffix description strips; and a description absent from
  `ptmlist.txt` (warns, produces no output, does not crash).
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
  deletion (36–38).
- **Still-unhandled feature types** (`sequence conflict`, `non-standard amino acid`) are
  parsed without error and produce no output.

Entry **Q67890** (TrEMBL) is a minimal, annotation-free entry: it checks multi-entry
handling, the annotation-free `>` line, a single `# Prefix=` for the whole file (the TrEMBL
entry is written as `>sp:` to stay consistent with the header), and non-ASCII text in the
protein name (sanitized to ASCII, as PEFF requires).

## Integration test (optional, real data)

The fixture above is a controlled unit test. For an end-to-end check at scale, convert the
whole reviewed human proteome. From a working directory holding the **full** `ptmlist.txt`:

```bash
curl -O 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt'
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
