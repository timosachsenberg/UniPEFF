# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

UniPEFF is a .NET 8 command-line tool (proteomics/bioinformatics) that reads a UniProtKB protein database in XML and is intended to generate a PEFF file (PSI Extended FASTA Format — annotated protein sequences used in mass spectrometry). Author: David L. Tabb, UMC Groningen.

Current state: the tool parses UniProt XML and **writes PEFF output**. Two output modes (see `-AnnotationIdentifiers`) cover the de-facto bulk-export convention and the spec-canonical annotation-identifier form. It still also dumps the last entry to the console via `DebugPrint()`. Version string is `"20260413 alpha"`.

## Build & run

```bash
dotnet build                     # or: dotnet build UniPEFF.sln
dotnet run -- -in path/to/uniprot_sprot.xml.gz
```

CLI arguments (parsed in `Program.Main`):
- `-in <file>` — input UniProt XML, either `.xml` or gzipped `.xml.gz` (extension determines decompression)
- `-out <file>` — PEFF output path (default: the input path with its `.xml`/`.xml.gz` extension replaced by `.peff`)
- `-prefix <p>` — force a single `# Prefix=`/`>prefix:accession` for the whole file. By default no override is given and entries are grouped by their dataset into one DB-description block per prefix (`sp` for Swiss-Prot, `tr` for TrEMBL); passing `-prefix` collapses everything into one block with that prefix.
- `-dbversion <v>` — value for the mandatory `# DbVersion=` header (default: `unknown`).
- `-AnnotationIdentifiers` — emit PEFF "Option B": `# HasAnnotationIdentifiers=true`, a sequential `id:` on every annotation tuple, and `\DisulfideBond` referencing the two half-cystine ids. Default ("Option C") emits half cystines as position-based `\ModResPsi` and no `\DisulfideBond`, matching UniProt/neXtProt bulk exports.
- `-OmitMolecularProcessing`, `-OmitAminoAcidModifications`, `-OmitSequenceVariations` — disable each annotation category

The entire program is a single file: `Program.cs`. There is no lint config and no unit-test framework, but `testdata/` holds a comprehensive fixture and a regression runner — see Testing below.

## Testing

`testdata/` is a self-contained regression test (see `testdata/README.md` for the full layout):
- `comprehensive.xml` — one Swiss-Prot "kitchen-sink" entry exercising every parser path (all five molecular-processing types, every modification routing incl. cross-link/glycosylation/lipidation, the `(Microbial infection)`/`;` description strips and a not-in-`ptmlist` warning, all three disulfide shapes, all four variant shapes, out-of-bounds variant guards, entry metadata, and an `alternative products` comment + `splice variant` feature that must be **ignored**) plus a minimal TrEMBL entry, and every writer path (escaping, wrapping, both modes).
- `compact.xml` — the same kind of input with no pretty-print whitespace, to prove the depth-driven parser is layout-independent.
- `ptmlist.txt` — minimal PTM CV covering only the fixture's mods. `psi-mod.obo` / `unimod.obo` — minimal OBO name maps for the canonical-name lookup.
- `expected_A.peff` / `expected_B.peff` — golden outputs (default and `-AnnotationIdentifiers`).
- `run.sh` — builds, runs both modes and diffs against the goldens, parses `compact.xml` and checks key annotations survive, and asserts spec-conformance independently of the goldens (`# DbVersion=` present, ASCII-only, one `# //`-separated DB block per dataset, TrEMBL under `>tr:`, and no isoform-id descriptor leaks). Run with `DOTNET=/path/to/dotnet ./run.sh` if `dotnet` is not on PATH. If a change intentionally alters output, regenerate the goldens with the two commands inside `run.sh`.

## Required runtime input

When amino-acid modifications are enabled (the default), the program reads **`ptmlist.txt` from the current working directory** at startup. This is UniProt's PTM controlled vocabulary, downloaded from:
`https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt`
Modified-residue features whose description is not found in this file are skipped with a warning. One PTM ("Half cystine", for disulfides) is hard-coded because it is absent from `ptmlist.txt`.

It also tries to read **`psi-mod.obo`** and **`unimod.obo`** from the working directory to resolve the strictly PEFF-required OBO `name:` for each `\ModResPsi`/`\ModResUnimod` accession. These files are *optional*: if absent (or missing an entry) the writer falls back to the UniProt `ptmlist.txt` name and prints a warning — the output is still produced but the modification names are not the canonical OBO names.

## Architecture

`Program.cs` holds everything. Three layers:

1. **`Program.Main`** — arg parsing, loads `ptmlist.txt` into a `PTMList`, opens the input (gzip or plain), hands an `XmlReader` to the parser, then calls the writer.
2. **`PEFFmodel.FromUniProtXML`** — the parser. A single `while (InputStream.Read())` loop over a forward-only `XmlReader`, dispatching on element name. `<entry>` starts a new `PEFFentry`; each `<feature>` subtree is slurped once by `ReadFeatureContent` into a `FeatureContent` struct (position/begin/end/original/variation), and the feature is then sorted into the entry's annotation lists by its `type` attribute and mapped to PEFF CV accessions. Only a `<sequence>` carrying a `length` attribute is taken as the canonical sequence, so isoform `<sequence>` elements inside `alternative products` comments are ignored. A final pass calls `MergeDisulfideResiduesIntoModifiedResidues` on each entry.
3. **`PEFFmodel.WritePeff` / `PEFFentry.WritePeffEntry`** — the writer. A `# PEFF 1.0` file-description block, then one DB-description block per dataset prefix (entries grouped `sp`/`tr`, each with its own `# NumberOfEntries=`; `-prefix` collapses this to a single block). Each entry is one `>prefix:accession` descriptor line (tags omitted when absent) with the sequence wrapped at 60 columns. `EscapePeff` backslash-escapes the four PEFF-reserved characters `\ | ( )`; unknown positions (stored as 0) become `?` via `Pos`. Half cystines reach the `\ModResPsi` list through the merge pass at the end of `FromUniProtXML`, which is also what makes the Option-B `\DisulfideBond` id references resolvable.

### Two conventions that pervade the code

- **Hand-rolled singly-linked lists, not `List<T>`.** Every model class (`PEFFentry`, `MolecularProcessing`, `ModifiedResidueUniProt`, `LinkedModifications`, `SequenceVariantSimple`, `SequenceVariantComplex`, `PTMList`) is a node with a `Next` pointer. Each list has a **sentinel head node** whose own fields are unused; real data starts at `head.Next`. The parser keeps a `*Runner` cursor per list and appends via `Runner.Next = new(...); Runner = Runner.Next`. `DebugPrint()` and `PTMList.Find()` both skip the head by starting at `this.Next`.

- **Depth-driven XML reading, not positional Read() counts.** An earlier version advanced by counting `InputStream.Read()` calls against an exact element layout; that was replaced. The parser now bounds every subtree scan by `XmlReader.Depth` (`int d = InputStream.Depth; while (Read() && !(NodeType == EndElement && Depth == d)) …`), used for `<feature>` (via `ReadFeatureContent`), `<gene>`, and `<organism>`. This is why `compact.xml` (no inter-element whitespace) parses identically. When editing the parser, keep working in terms of element name + depth rather than reintroducing fixed read counts.

### Model relationships

- `PEFFmodel` → linked list of `PEFFentry` (one per protein).
- Each `PEFFentry` owns six annotation lists plus scalar fields: primary accession, name, gene, base sequence, plus the metadata that becomes `\NcbiTaxId`/`\TaxName`/`\SV`/`\EV`/`\PE`/`\ID` and the dataset (→ `sp`/`tr` prefix). The first `<accession>` is the primary id used in `>prefix:accession`; any further accessions are kept for `\AltAC`. Only the first `<name>`/`<fullName>` is kept.
- Disulfide bonds with two endpoints create two `ModifiedResidueUniProt` (Half cystine) nodes linked by a `LinkedModifications` node; single-endpoint disulfides (cross-chain) create only one residue node.

### Known incomplete areas

- Feature `type` values with no `case` (e.g. `non-standard amino acid`, `sequence conflict`, `splice variant`) fall through unhandled and are omitted. (`cross-link`, `glycosylation site`, and `lipid moiety-binding region` are now handled — see the modification switch in `FromUniProtXML`.)
- **Isoforms are not expanded.** `alternative products` comments and their isoform `<sequence>`/`splice variant` features are deliberately ignored; only the canonical sequence is emitted, and no `accession-N` isoform descriptor is produced. A `run.sh` assertion guards against isoform-id leakage.
- Carbamidomethyl and other sample-prep modifications are intentionally **not** emitted — they are search-time parameters, not native sequence annotations.
