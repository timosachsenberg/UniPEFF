using System;
using System.IO;
using System.IO.Compression;
using System.Xml;
using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text.Json;
using System.Text;

namespace UniPEFF
{
    class PEFFmodel
    {
        string Source;
        PEFFentry Entries = new PEFFentry();
        // Strip UniProt's "(Microbial infection) " prefix and any trailing "; ..." qualifier from a
        // PTM description, matching the long-standing logic in the "modified residue" case.
        static string CleanPtmDescription(string Description)
        {
            if (string.IsNullOrEmpty(Description)) return Description;
            // Strip the "(Microbial infection) " qualifier by matching the whole literal (incl. its
            // trailing space), so the substring length can never overrun a bare "(Microbial infection)".
            const string MicrobialInfection = "(Microbial infection) ";
            if (Description.StartsWith(MicrobialInfection)) Description = Description.Substring(MicrobialInfection.Length);
            int SemiPosition = Description.IndexOf(';');
            if (SemiPosition > -1) Description = Description.Substring(0, SemiPosition);
            return Description;
        }

        // The parsed content of a <feature> subtree: the original/variation residues (sequence
        // variants) and the location. A null Begin/End/Position means UniProt gave the coordinate
        // as status="unknown" (or omitted it); callers decide whether that is legal for their key
        // (ModRes accepts "?"; Variant/Processed positions MUST be known and in range).
        class FeatureContent
        {
            public string Original;
            public string Variation;
            public bool HasPosition;   // a single <position> was present
            public bool HasRange;      // a <begin>/<end> pair was present
            public int? Position;
            public int? Begin;
            public int? End;
        }

        // Parse a UniProt "position" attribute safely: returns null for a missing attribute or a
        // status="unknown" coordinate (where the attribute is absent). Never throws.
        static int? ParsePosition(string Attribute)
            => (Attribute != null && Int32.TryParse(Attribute, out int Value)) ? Value : (int?)null;

        // PEFF VariantSimple newAminoAcid MUST be an amino-acid letter (ambiguity codes included) or '*'.
        public static bool IsResidueCode(char C)
            => (C >= 'A' && C <= 'Z') || (C >= 'a' && C <= 'z') || C == '*';

        // Read a <feature> element's subtree by NAME and DEPTH (not by counting Read() calls),
        // collecting the location and any original/variation residues. The reader must be on the
        // <feature> start element; on return it is on the matching </feature> end element. Robust to
        // whitespace, element ordering, optional elements, and status="unknown".
        static FeatureContent ReadFeatureContent(XmlReader InputStream)
        {
            var Content = new FeatureContent();
            if (InputStream.IsEmptyElement) return Content;
            int FeatureDepth = InputStream.Depth;
            while (InputStream.Read() && !(InputStream.NodeType == XmlNodeType.EndElement && InputStream.Depth == FeatureDepth))
            {
                if (InputStream.NodeType != XmlNodeType.Element) continue;
                switch (InputStream.Name)
                {
                    case "original":
                        if (!InputStream.IsEmptyElement && InputStream.Read()) Content.Original = InputStream.Value;
                        break;
                    case "variation":
                        if (Content.Variation == null && !InputStream.IsEmptyElement && InputStream.Read()) Content.Variation = InputStream.Value;
                        break;
                    case "position":
                        Content.HasPosition = true;
                        Content.Position = ParsePosition(InputStream["position"]);
                        break;
                    case "begin":
                        Content.HasRange = true;
                        Content.Begin = ParsePosition(InputStream["position"]);
                        break;
                    case "end":
                        Content.HasRange = true;
                        Content.End = ParsePosition(InputStream["position"]);
                        break;
                }
            }
            return Content;
        }

        public static PEFFmodel FromUniProtXML(XmlReader InputStream, bool RecordMolecularProcessing, bool RecordAminoAcidModifications, bool RecordSequenceVariations, PTMList PTMCV)
        {
            Console.WriteLine("Reading XML...");
            PEFFmodel ThisPEFF = new PEFFmodel();
            PEFFentry PERunner = ThisPEFF.Entries;
            MolecularProcessing    MPRunner = null;
            ModifiedResidueUniProt MRURunner = null;
            ModifiedResidueUniProt DSRunner = null;
            SequenceVariantSimple  SVSRunner = null;
            SequenceVariantComplex SVCRunner = null;
            LinkedModifications    DSLMRunner = null;
            AltAccession           AAARunner = null;

            // Append helpers for the sentinel-headed lists: each advances the captured runner cursor.
            MolecularProcessing AddProcessing(string CV, string Type, int Begin, int End)
            {
                MPRunner.Next = new MolecularProcessing();
                MPRunner = MPRunner.Next;
                MPRunner.CV = CV; MPRunner.Type = Type; MPRunner.Begin = Begin; MPRunner.End = End;
                return MPRunner;
            }
            ModifiedResidueUniProt AddModifiedResidue(int Position, PTMList Modification)
            {
                MRURunner.Next = new ModifiedResidueUniProt();
                MRURunner = MRURunner.Next;
                MRURunner.Position = Position; MRURunner.Modification = Modification;
                return MRURunner;
            }
            ModifiedResidueUniProt AddDisulfideResidue(int Position, PTMList Modification)
            {
                DSRunner.Next = new ModifiedResidueUniProt();
                DSRunner = DSRunner.Next;
                DSRunner.Position = Position; DSRunner.Modification = Modification;
                return DSRunner;
            }

            while (InputStream.Read())
            {
                var ThisNodeType = InputStream.NodeType;
                if (ThisNodeType == XmlNodeType.Element)
                {
                    if (InputStream.Name == "entry")
                    {
                        //This is the start of a new entry
                        PERunner.Next = new PEFFentry();
                        PERunner  = PERunner.Next;
                        MPRunner  = PERunner.MolecularProcessings;
                        MRURunner = PERunner.ModifiedResiduesUniProt;
                        DSRunner  = PERunner.ModifiedResiduesDisulfide;
                        DSLMRunner = PERunner.Disulfides;
                        SVSRunner = PERunner.SequenceVariantsSimple;
                        SVCRunner = PERunner.SequenceVariantsComplex;
                        AAARunner = PERunner.AltAccessions;
                        PERunner.DataSet = InputStream["dataset"];
                        PERunner.EntryVersion = InputStream["version"];
                    }
                    else if (InputStream.Name == "accession")
                    {
                        // The first accession is the primary id (in >Prefix:Accession); the rest -> \AltAC.
                        if (InputStream.Read())
                        {
                            if (PERunner.Accession == null) PERunner.Accession = InputStream.Value;
                            else
                            {
                                AAARunner.Next = new AltAccession();
                                AAARunner = AAARunner.Next;
                                AAARunner.Accession = InputStream.Value;
                            }
                        }
                    }
                    else if (InputStream.Name == "name")
                    {
                        if (PERunner.Name == null)
                        {
                            if (InputStream.Read())
                            {
                                PERunner.Name = InputStream.Value;
                            }
                        }
                    }
                    else if (InputStream.Name == "fullName")
                    {
                        if (PERunner.FullName == null)
                        {
                            if (InputStream.Read())
                            {
                                PERunner.FullName = InputStream.Value;
                            }
                        }
                    }
                    else if (InputStream.Name == "gene")
                    {
                        // Bounded scan of the <gene> subtree: capture <name type="primary"> wherever it
                        // appears (a synonym/ORF/ordered-locus name may precede it), stopping at </gene>.
                        if (!InputStream.IsEmptyElement)
                        {
                            int GeneDepth = InputStream.Depth;
                            while (InputStream.Read() && !(InputStream.NodeType == XmlNodeType.EndElement && InputStream.Depth == GeneDepth))
                            {
                                if (InputStream.NodeType == XmlNodeType.Element && InputStream.Name == "name"
                                    && InputStream["type"] == "primary" && PERunner.PrimaryGene == null)
                                {
                                    if (InputStream.Read()) PERunner.PrimaryGene = InputStream.Value;
                                }
                            }
                        }
                    }
                    else if ((InputStream.Name == "sequence") && (InputStream["length"] != null))
                    {
                        PERunner.SequenceVersion = InputStream["version"];
                        if (InputStream.Read()) PERunner.BaseSequence = InputStream.Value;
                    }
                    else if (InputStream.Name == "organism")
                    {
                        // Bounded scan of the <organism> subtree: capture the scientific <name> (in any
                        // order) and the NCBI Taxonomy id, stopping at </organism> so we never walk into
                        // the entry body, and so the organism's own <name> elements cannot leak back to
                        // the bare "name" branch above.
                        int OrganismDepth = InputStream.Depth;
                        while (InputStream.Read() && !(InputStream.NodeType == XmlNodeType.EndElement && InputStream.Depth == OrganismDepth))
                        {
                            if (InputStream.NodeType != XmlNodeType.Element) continue;
                            if (InputStream.Name == "name" && InputStream["type"] == "scientific" && PERunner.TaxName == null)
                            {
                                if (InputStream.Read()) PERunner.TaxName = InputStream.Value;
                            }
                            else if (InputStream.Name == "dbReference" && InputStream["type"] == "NCBI Taxonomy")
                            {
                                PERunner.NcbiTaxId = InputStream["id"];
                            }
                        }
                    }
                    else if (InputStream.Name == "proteinExistence")
                    {
                        PERunner.ProteinExistence = PEFFentry.ProteinExistenceCode(InputStream["type"]);
                    }
                    if (InputStream.Name == "feature")
                    {
                        string FeatureType = InputStream["type"];
                        string FeatureDescription = InputStream["description"];
                        // Read the whole feature subtree once, by name/depth, so every category below
                        // works from structured content instead of fragile Read()-count navigation.
                        FeatureContent FC = ReadFeatureContent(InputStream);
                        // A single-residue site: prefer <position>, else a range's <begin>; an unknown
                        // coordinate becomes 0, which the writer renders as "?" (legal only in ModRes).
                        int ModPosition = FC.Position ?? FC.Begin ?? 0;

                        if (RecordMolecularProcessing)
                        {
                            switch (FeatureType)
                            {
                                case "chain":
                                    AddProcessing("PEFF:0001020", "mature protein", FC.Begin ?? 0, FC.End ?? 0);
                                    break;
                                case "initiator methionine":
                                    AddProcessing("PEFF:0001035", "initiator methionine", FC.Position ?? 0, FC.Position ?? 0);
                                    break;
                                case "propeptide":
                                    AddProcessing("PEFF:0001034", "propeptide", FC.Begin ?? 0, FC.End ?? 0);
                                    break;
                                case "signal peptide":
                                    AddProcessing("PEFF:0001021", "signal peptide", FC.Begin ?? 0, FC.End ?? 0);
                                    break;
                                case "transit peptide":
                                    AddProcessing("PEFF:0001022", "transit peptide", FC.Begin ?? 0, FC.End ?? 0);
                                    break;
                            }
                        }
                        if (RecordAminoAcidModifications)
                        {
                            switch (FeatureType)
                            {
                                case "cross-link":
                                    // Isopeptide cross-link (ubiquitin/SUMO etc.). PEFF has no generic
                                    // cross-link-bond key, so annotate the modified residue(s) by position
                                    // as a generic \ModRes -- a two-residue bond is not a residue PTM, so it
                                    // is never routed into \ModResPsi/\ModResUnimod. Endpoints are not paired.
                                    {
                                        string Desc = CleanPtmDescription(FeatureDescription);
                                        if (string.IsNullOrEmpty(Desc)) Desc = "cross-link";
                                        PTMList XLink = PTMList.Synthetic(Desc);
                                        if (FC.HasRange)
                                        {
                                            AddModifiedResidue(FC.Begin ?? 0, XLink);
                                            AddModifiedResidue(FC.End ?? 0, XLink);
                                        }
                                        else AddModifiedResidue(ModPosition, XLink);
                                    }
                                    break;
                                case "disulfide bond":
                                    // <begin>/<end> => intra-chain bond (two linked half cystines); a single
                                    // <position> => inter-chain (one endpoint here, not linked). Half cystine
                                    // is hard-coded because it is absent from ptmlist.txt.
                                    {
                                        PTMList HalfCystine = PTMCV.Find("Half cystine");
                                        if (FC.HasRange)
                                        {
                                            var First  = AddDisulfideResidue(FC.Begin ?? 0, HalfCystine);
                                            var Second = AddDisulfideResidue(FC.End ?? 0, HalfCystine);
                                            DSLMRunner.Next = new LinkedModifications();
                                            DSLMRunner = DSLMRunner.Next;
                                            DSLMRunner.Residue1 = First;
                                            DSLMRunner.Residue2 = Second;
                                        }
                                        else AddDisulfideResidue(ModPosition, HalfCystine);
                                    }
                                    break;
                                case "glycosylation site":
                                    // Glycans are not in PSI-MOD/Unimod -> generic \ModRes with the description.
                                    {
                                        string Desc = CleanPtmDescription(FeatureDescription);
                                        if (string.IsNullOrEmpty(Desc)) Desc = "glycosylation site";
                                        AddModifiedResidue(ModPosition, PTMList.Synthetic(Desc));
                                    }
                                    break;
                                case "lipid moiety-binding region":
                                    // Lipidation -> \ModResPsi/\ModResUnimod when it is a known ptmlist PTM,
                                    // else a generic \ModRes.
                                    {
                                        string Desc = CleanPtmDescription(FeatureDescription);
                                        if (string.IsNullOrEmpty(Desc)) Desc = "lipid moiety-binding region";
                                        AddModifiedResidue(ModPosition, PTMCV.Find(Desc) ?? PTMList.Synthetic(Desc));
                                    }
                                    break;
                                case "modified residue":
                                    {
                                        string Desc = CleanPtmDescription(FeatureDescription);
                                        PTMList CVHit = PTMCV.Find(Desc);
                                        if (CVHit != null) AddModifiedResidue(ModPosition, CVHit);
                                        else
                                        {
                                            Console.WriteLine("Failed to find this modified residue in ptmlist.txt:\t" + Desc);
                                            Console.WriteLine("Please update your copy of ptmlist.txt from here:");
                                            Console.WriteLine("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt");
                                        }
                                    }
                                    break;
                            }
                        }
                        if (RecordSequenceVariations && FeatureType == "sequence variant")
                        {
                            // newSequence is the <variation> (empty = a deletion); the location is a single
                            // <position> or a <begin>/<end> range. A single-<position> change is VariantSimple
                            // only when the new residue is exactly one valid amino-acid letter; an insertion,
                            // a multi-residue change, or a deletion at one point becomes VariantComplex. An
                            // unknown position is skipped ("?" is ModRes-only, illegal for variants).
                            string NewSeq = FC.Variation ?? "";
                            if (FC.HasRange)
                            {
                                if (FC.Begin == null || FC.End == null)
                                    Console.Error.WriteLine("\tWarning: " + (PERunner.Accession ?? "?") + " sequence variant with unknown range; omitted.");
                                else
                                {
                                    SVCRunner.Next = new SequenceVariantComplex();
                                    SVCRunner = SVCRunner.Next;
                                    SVCRunner.Begin = FC.Begin.Value;
                                    SVCRunner.End = FC.End.Value;
                                    SVCRunner.NewSeq = NewSeq;
                                }
                            }
                            else if (FC.HasPosition)
                            {
                                if (FC.Position == null)
                                    Console.Error.WriteLine("\tWarning: " + (PERunner.Accession ?? "?") + " sequence variant with unknown position; omitted.");
                                else if (NewSeq.Length == 1 && IsResidueCode(NewSeq[0]))
                                {
                                    SVSRunner.Next = new SequenceVariantSimple();
                                    SVSRunner = SVSRunner.Next;
                                    SVSRunner.NewAA = NewSeq[0];
                                    SVSRunner.Position = FC.Position.Value;
                                }
                                else
                                {
                                    SVCRunner.Next = new SequenceVariantComplex();
                                    SVCRunner = SVCRunner.Next;
                                    SVCRunner.Begin = FC.Position.Value;
                                    SVCRunner.End = FC.Position.Value;
                                    SVCRunner.NewSeq = NewSeq;
                                }
                            }
                        }
                    }
                }
            }
            // Both source lists are complete only now; fold disulfide cystines into each
            // entry's modified-residue list (position-sorted) before handing back the model.
            for (var EntryRunner = ThisPEFF.Entries.Next; EntryRunner != null; EntryRunner = EntryRunner.Next)
            {
                EntryRunner.MergeDisulfideResiduesIntoModifiedResidues();
            }
            PERunner.DebugPrint();
            return ThisPEFF;
        }

        public void WritePeff(TextWriter Writer, PeffHeader Header, string PrefixOverride, bool AnnotationIdentifiers, OboNameMap PsiModNames, OboNameMap UnimodNames)
        {
            Header.HasAnnotationIdentifiers = AnnotationIdentifiers;

            // A PEFF entry is a description line PLUS a sequence block, and MUST start with
            // >Prefix:DbUniqueId, so entries lacking an accession OR a sequence are skipped (with a
            // warning); if that leaves nothing, write no (invalid, empty) PEFF.
            int Writables = 0;
            for (var ER = Entries.Next; ER != null; ER = ER.Next)
            {
                if (ER.IsWritable()) Writables++;
                else if (string.IsNullOrEmpty(ER.Accession)) Console.Error.WriteLine("\tWarning: skipping an entry with no accession.");
                else Console.Error.WriteLine("\tWarning: skipping entry " + ER.Accession + " with no sequence.");
            }
            if (Writables == 0)
            {
                Console.Error.WriteLine("\tError: no writable entries; no PEFF written (a PEFF file must contain >= 1 sequence entry).");
                return;
            }

            WriteFileDescriptionBlock(Writer);
            int Fallbacks = 0;

            // Group entries by database prefix (sp/tr, from the dataset) and emit one database-
            // description block per group -- each with its own per-block NumberOfEntries -- followed
            // by that group's entries under their own prefix. A -prefix override forces one block
            // over all entries. (Keying on the prefix, not the raw dataset, keeps null/Swiss-Prot/
            // unknown datasets in a single "sp" block rather than producing duplicate prefixes.)
            bool Single = !string.IsNullOrEmpty(PrefixOverride);
            var Prefixes = new List<string>();
            if (Single) Prefixes.Add(PrefixOverride);
            else
            {
                for (var ER = Entries.Next; ER != null; ER = ER.Next)
                {
                    if (!ER.IsWritable()) continue;
                    string P = PEFFentry.PrefixForDataset(ER.DataSet);
                    if (!Prefixes.Contains(P)) Prefixes.Add(P);
                }
            }
            // Header section: ALL database-description blocks come first -- the spec requires the
            // whole header section to precede the sequence-entry section (the entries then follow).
            foreach (string P in Prefixes)
            {
                int Count = 0;
                for (var ER = Entries.Next; ER != null; ER = ER.Next)
                    if (ER.IsWritable() && (Single || PEFFentry.PrefixForDataset(ER.DataSet) == P)) Count++;
                Header.Prefix = P;
                Header.NumberOfEntries = Count;
                WriteDbDescriptionBlock(Writer, Header);
            }
            // Sequence-entry section: entries grouped by their database prefix.
            foreach (string P in Prefixes)
            {
                for (var ER = Entries.Next; ER != null; ER = ER.Next)
                    if (ER.IsWritable() && (Single || PEFFentry.PrefixForDataset(ER.DataSet) == P))
                        Fallbacks += ER.WritePeffEntry(Writer, P, AnnotationIdentifiers, PsiModNames, UnimodNames);
            }
            if (Fallbacks > 0)
                Console.Error.WriteLine("\tWARNING: " + Fallbacks + " ModResPsi/ModResUnimod name(s) had no OBO 'name:' (psi-mod.obo/unimod.obo absent or incomplete); emitted the UniProt ptmlist name instead, which is NOT the strictly PEFF-required OBO name. Place the OBO files in the working directory for conformant output.");
        }

        static void WriteFileDescriptionBlock(TextWriter Writer)
        {
            Writer.Write("# PEFF 1.0\n");
            Writer.Write("# GeneralComment=Generated by UniPEFF\n");
            Writer.Write("# //\n");
        }

        static void WriteDbDescriptionBlock(TextWriter Writer, PeffHeader Header)
        {
            Writer.Write("# DbName=" + Header.DbName + "\n");
            Writer.Write("# Prefix=" + Header.Prefix + "\n");
            if (!string.IsNullOrEmpty(Header.DbDescription)) Writer.Write("# DbDescription=" + Header.DbDescription + "\n");
            Writer.Write("# Decoy=false\n");
            Writer.Write("# DbSource=" + Header.DbSource + "\n");      // mandatory per PEFF 1.0 sec 3.3.2
            Writer.Write("# DbVersion=" + Header.DbVersion + "\n");    // mandatory per PEFF 1.0 sec 3.3.2
            if (!string.IsNullOrEmpty(Header.DbDate)) Writer.Write("# DbDate=" + Header.DbDate + "\n");
            Writer.Write("# NumberOfEntries=" + Header.NumberOfEntries + "\n");
            Writer.Write("# SequenceType=AA\n");
            if (Header.HasAnnotationIdentifiers) Writer.Write("# HasAnnotationIdentifiers=true\n");
            Writer.Write("# //\n");
        }
    }
    class AltAccession
    {
        public string Accession;
        public AltAccession Next = null;
        public void DebugPrint()
        {
            AltAccession Runner = this.Next;
            while (Runner != null) { Console.WriteLine("AltAccession\t" + Runner.Accession); Runner = Runner.Next; }
        }
    }
    class PEFFentry
    {
        public string DataSet;
        public string Accession;
        public string Name;
        public string FullName;
        public string PrimaryGene;
        public string BaseSequence;
        public string NcbiTaxId;          // NCBI taxonomy id, e.g. "9606"
        public string TaxName;            // scientific organism name, e.g. "Homo sapiens"
        public string SequenceVersion;    // <sequence version=> -> \SV
        public string EntryVersion;       // <entry version=>    -> \EV
        public string ProteinExistence;   // PE digit "1".."5"   -> \PE
        public AltAccession AltAccessions = new AltAccession();   // secondary accessions -> \AltAC
        public MolecularProcessing MolecularProcessings = new MolecularProcessing();
        public ModifiedResidueUniProt ModifiedResiduesUniProt = new ModifiedResidueUniProt();
        public ModifiedResidueUniProt ModifiedResiduesDisulfide = new ModifiedResidueUniProt();
        public LinkedModifications Disulfides = new LinkedModifications();
        public SequenceVariantSimple SequenceVariantsSimple = new SequenceVariantSimple();
        public SequenceVariantComplex SequenceVariantsComplex = new SequenceVariantComplex();

        public PEFFentry Next = null;

        // Fold the disulfide-bond half cystines into the main modified-residue list,
        // position-sorted, so PEFF output can emit a single ordered ModRes annotation.
        public void MergeDisulfideResiduesIntoModifiedResidues()
        {
            ModifiedResiduesUniProt.AbsorbAndSortFrom(ModifiedResiduesDisulfide);
        }

        public static string PrefixForDataset(string DataSet)
        {
            if (DataSet == "TrEMBL") return "tr";
            return "sp";   // Swiss-Prot and the safe default
        }

        // A PEFF entry needs both a primary id (for >Prefix:DbUniqueId) and a sequence block.
        public bool IsWritable() => !string.IsNullOrEmpty(Accession) && !string.IsNullOrEmpty(BaseSequence);

        // Map UniProt's <proteinExistence type=> to the PEFF \PE digit (1-5); null if unrecognised.
        public static string ProteinExistenceCode(string Type) => Type switch
        {
            "evidence at protein level"    => "1",
            "evidence at transcript level" => "2",
            "inferred from homology"       => "3",
            "predicted"                    => "4",
            "uncertain"                    => "5",
            _ => null,
        };

        // Make text PEFF-safe. PEFF requires ASCII, so non-ASCII is transliterated/stripped first
        // (ToAscii). Then reserved characters are escaped -- but per the spec only "|", "\" and
        // UNPAIRED parentheses are escaped; balanced parentheses are left intact so readers can
        // handle embedded parens (e.g. "N-linked (GlcNAc...)").
        static string EscapePeff(string Value)
        {
            if (string.IsNullOrEmpty(Value)) return Value ?? "";
            Value = ToAscii(Value);
            // Mark unpaired parentheses; only those (plus | and \) get a backslash.
            var Unpaired = new HashSet<int>();
            var Open = new Stack<int>();
            for (int i = 0; i < Value.Length; i++)
            {
                if (Value[i] == '(') Open.Push(i);
                else if (Value[i] == ')') { if (Open.Count > 0) Open.Pop(); else Unpaired.Add(i); }
            }
            foreach (int i in Open) Unpaired.Add(i);

            var Builder = new StringBuilder(Value.Length);
            for (int i = 0; i < Value.Length; i++)
            {
                char C = Value[i];
                if (C == '\\' || C == '|' || ((C == '(' || C == ')') && Unpaired.Contains(i))) Builder.Append('\\');
                Builder.Append(C);
            }
            return Builder.ToString();
        }

        // PEFF allows only ASCII characters. Normalize (so accents such as e-acute collapse to
        // "e"), drop combining marks, and replace any remaining non-ASCII (e.g. Greek letters)
        // with "?" (a character explicitly permitted in PEFF item text).
        static string ToAscii(string Value)
        {
            bool AlreadyAscii = true;
            foreach (char C in Value) if (C > '\x7F') { AlreadyAscii = false; break; }
            if (AlreadyAscii) return Value;
            string Nfkd = Value.Normalize(NormalizationForm.FormKD);
            var Builder = new StringBuilder(Nfkd.Length);
            foreach (char C in Nfkd)
            {
                if (CharUnicodeInfo.GetUnicodeCategory(C) == UnicodeCategory.NonSpacingMark) continue;
                Builder.Append(C <= '\x7F' ? C : '?');
            }
            return Builder.ToString();
        }

        // PEFF encodes an unknown sequence position (our parser leaves it 0) as "?".
        static string Pos(int Position) => Position == 0 ? "?" : Position.ToString();

        // Write one PEFF entry: the ">" descriptor line (tags omitted when we have no data for
        // them) followed by the sequence wrapped at a fixed width. When AnnotationIdentifiers is
        // set, every annotation tuple carries a sequential "id:" prefix and disulfide bonds are
        // emitted as \DisulfideBond referencing the two half-cystine annotation ids (PEFF
        // "Option B"). Otherwise half cystines are written as plain ModResPsi and no bond
        // connectivity is emitted (PEFF "Option C", matching UniProt/neXtProt bulk exports).
        public int WritePeffEntry(TextWriter Writer, string Prefix, bool AnnotationIdentifiers, OboNameMap PsiModNames, OboNameMap UnimodNames)
        {
            const int SequenceWidth = 60;
            int NextId = 0;
            int Fallbacks = 0;
            var ResidueId = new Dictionary<ModifiedResidueUniProt, int>();

            void Tag(string Key, string ValueText) => Writer.Write(" \\" + Key + "=" + ValueText);
            string IdPrefix() => AnnotationIdentifiers ? (NextId++) + ":" : "";

            Writer.Write(">" + Prefix + ":" + (Accession ?? ""));
            if (!string.IsNullOrEmpty(FullName))         Tag("PName", EscapePeff(FullName));   // single-value key: bare, not a (list)
            if (!string.IsNullOrEmpty(PrimaryGene))      Tag("GName", EscapePeff(PrimaryGene));
            if (!string.IsNullOrEmpty(NcbiTaxId))        Tag("NcbiTaxId", EscapePeff(NcbiTaxId));
            if (!string.IsNullOrEmpty(TaxName))          Tag("TaxName", EscapePeff(TaxName));
            if (!string.IsNullOrEmpty(BaseSequence))     Tag("Length", BaseSequence.Length.ToString());
            if (!string.IsNullOrEmpty(SequenceVersion))  Tag("SV", EscapePeff(SequenceVersion));
            if (!string.IsNullOrEmpty(EntryVersion))     Tag("EV", EscapePeff(EntryVersion));
            if (!string.IsNullOrEmpty(ProteinExistence)) Tag("PE", ProteinExistence);
            // \DbUniqueId is intentionally NOT emitted: the PSI-MS CV (PEFF:0001001) states it
            // "shall not be used in the PEFF 1.0 serialization as it is redundant with the primary
            // identifier following the >".
            if (!string.IsNullOrEmpty(Name))             Tag("ID", EscapePeff(Name));
            {
                // Secondary accessions as one \AltAC=(ac)(ac) list (a key MUST NOT appear twice).
                var AltBuilder = new StringBuilder();
                for (var A = AltAccessions.Next; A != null; A = A.Next)
                    AltBuilder.Append('(').Append(EscapePeff(A.Accession)).Append(')');
                if (AltBuilder.Length > 0) Tag("AltAC", AltBuilder.ToString());
            }

            // Modified residues split into the three PEFF namespaces, emitted (and numbered) in
            // the order PSI-MOD, Unimod, generic. Each node's id is recorded so \DisulfideBond
            // can reference the half cystines by id.
            void WriteMods(string Key, Func<ModifiedResidueUniProt, bool> Belongs, Func<ModifiedResidueUniProt, string> AccessionOf, OboNameMap OboNames)
            {
                var Builder = new StringBuilder();
                for (var M = ModifiedResiduesUniProt.Next; M != null; M = M.Next)
                {
                    if (!Belongs(M)) continue;
                    Builder.Append('(');
                    if (AnnotationIdentifiers) { ResidueId[M] = NextId; Builder.Append(NextId + ":"); NextId++; }
                    Builder.Append(Pos(M.Position));
                    string Accession = AccessionOf(M);
                    Builder.Append('|').Append(EscapePeff(Accession));
                    // PEFF requires the OBO "name:" field for ModResPsi/ModResUnimod (spec sec 3.4),
                    // not the UniProt ptmlist ID. Resolve it from the loaded CV; on a miss fall back
                    // to the UniProt ID (an empty name is illegal).
                    string DisplayName = M.Modification != null ? (M.Modification.ID ?? "") : "";
                    if (OboNames != null)
                    {
                        string OboName = OboNames.Find(Accession);
                        if (OboName != null) DisplayName = OboName;
                        else Fallbacks++;   // count every miss, incl. when the OBO file is entirely absent
                    }
                    if (!string.IsNullOrEmpty(DisplayName)) Builder.Append('|').Append(EscapePeff(DisplayName));
                    Builder.Append(')');
                }
                if (Builder.Length > 0) Tag(Key, Builder.ToString());
            }
            WriteMods("ModResPsi",
                      M => M.Modification != null && !string.IsNullOrEmpty(M.Modification.PSIModAccession),
                      M => M.Modification.PSIModAccession, PsiModNames);
            WriteMods("ModResUnimod",
                      M => M.Modification != null && string.IsNullOrEmpty(M.Modification.PSIModAccession) && !string.IsNullOrEmpty(M.Modification.UnimodAccession),
                      M => "UNIMOD:" + M.Modification.UnimodAccession, UnimodNames);
            WriteMods("ModRes",
                      M => M.Modification == null || (string.IsNullOrEmpty(M.Modification.PSIModAccession) && string.IsNullOrEmpty(M.Modification.UnimodAccession)),
                      M => "", null);

            // Simple (single-residue) variants: the real list, plus any single-residue substitution
            // that UniProt expressed as a range (PEFF requires those to be VariantSimple, not Complex).
            {
                var Builder = new StringBuilder();
                for (var V = SequenceVariantsSimple.Next; V != null; V = V.Next)
                {
                    if (V.Position == 0) continue;   // PEFF VariantSimple position MUST be > 0; "?" is ModRes-only
                    if (BaseSequence != null && V.Position > BaseSequence.Length)
                    {
                        Console.Error.WriteLine("\tWarning: " + Accession + " VariantSimple position " + V.Position + " exceeds sequence length " + BaseSequence.Length + "; omitted.");
                        continue;   // spec: position MUST be <= protein length
                    }
                    Builder.Append('(').Append(IdPrefix()).Append(V.Position).Append('|').Append(EscapePeff(V.NewAA.ToString())).Append(')');
                }
                for (var V = SequenceVariantsComplex.Next; V != null; V = V.Next)
                {
                    // (b|b|X) single-residue substitution -> VariantSimple; a one-char deletion (b|b|) stays complex.
                    if (V.Begin != 0 && V.Begin == V.End && V.NewSeq != null && V.NewSeq.Length == 1 && PEFFmodel.IsResidueCode(V.NewSeq[0])
                        && !(BaseSequence != null && V.Begin > BaseSequence.Length))
                        Builder.Append('(').Append(IdPrefix()).Append(V.Begin).Append('|').Append(EscapePeff(V.NewSeq)).Append(')');
                }
                if (Builder.Length > 0) Tag("VariantSimple", Builder.ToString());
            }
            // Complex (multi-residue / deletion) variants, excluding the single-residue substitutions
            // demoted to VariantSimple above.
            {
                var Builder = new StringBuilder();
                for (var V = SequenceVariantsComplex.Next; V != null; V = V.Next)
                {
                    if (V.Begin == 0 || V.End == 0) continue;   // positions count from 1; "?" is ModRes-only
                    if (V.Begin > V.End || (BaseSequence != null && (V.Begin > BaseSequence.Length || V.End > BaseSequence.Length)))
                    {
                        Console.Error.WriteLine("\tWarning: " + Accession + " VariantComplex " + V.Begin + "-" + V.End + " is out of bounds (length " + BaseSequence?.Length + "); omitted.");
                        continue;   // spec: positions count from 1 and must lie within the sequence
                    }
                    if (V.Begin == V.End && V.NewSeq != null && V.NewSeq.Length == 1 && PEFFmodel.IsResidueCode(V.NewSeq[0])) continue;   // demoted to VariantSimple
                    Builder.Append('(').Append(IdPrefix()).Append(V.Begin).Append('|').Append(V.End).Append('|').Append(EscapePeff(V.NewSeq ?? "")).Append(')');
                }
                if (Builder.Length > 0) Tag("VariantComplex", Builder.ToString());
            }
            // Molecular processing -> \Processed
            {
                var Builder = new StringBuilder();
                for (var P = MolecularProcessings.Next; P != null; P = P.Next)
                {
                    if (P.Begin == 0 || P.End == 0) continue;   // \Processed positions count from 1; "?" is ModRes-only
                    if (P.Begin > P.End || (BaseSequence != null && (P.Begin > BaseSequence.Length || P.End > BaseSequence.Length)))
                    {
                        Console.Error.WriteLine("\tWarning: " + Accession + " Processed " + P.Begin + "-" + P.End + " is out of bounds (length " + BaseSequence?.Length + "); omitted.");
                        continue;   // spec: positions count from 1 and must lie within the sequence
                    }
                    Builder.Append('(').Append(IdPrefix()).Append(P.Begin).Append('|').Append(P.End).Append('|').Append(EscapePeff(P.CV ?? "")).Append('|').Append(EscapePeff(P.Type ?? "")).Append(')');
                }
                if (Builder.Length > 0) Tag("Processed", Builder.ToString());
            }
            // Disulfide connectivity (Option B only): reference the two half-cystine ids
            if (AnnotationIdentifiers)
            {
                var Builder = new StringBuilder();
                for (var D = Disulfides.Next; D != null; D = D.Next)
                {
                    if (D.Residue1 == null || D.Residue2 == null) continue;
                    if (!ResidueId.TryGetValue(D.Residue1, out int Id1)) continue;
                    if (!ResidueId.TryGetValue(D.Residue2, out int Id2)) continue;
                    Builder.Append('(').Append(NextId++ + ":").Append(Id1).Append(',').Append(Id2).Append(')');
                }
                if (Builder.Length > 0) Tag("DisulfideBond", Builder.ToString());
            }

            Writer.Write("\n");
            string Sequence = BaseSequence ?? "";
            for (int Offset = 0; Offset < Sequence.Length; Offset += SequenceWidth)
                Writer.Write(Sequence.Substring(Offset, Math.Min(SequenceWidth, Sequence.Length - Offset)) + "\n");
            return Fallbacks;
        }

        public void DebugPrint ()
        {
            Console.WriteLine();
            Console.WriteLine(this.DataSet);
            Console.WriteLine(this.Accession);
            Console.WriteLine(this.Name);
            Console.WriteLine(this.FullName);
            Console.WriteLine(this.PrimaryGene);
            Console.WriteLine(this.BaseSequence);
            Console.WriteLine("-----Molecular Processings");
            this.MolecularProcessings.DebugPrint();
            Console.WriteLine("-----Sequence Variants Simple");
            this.SequenceVariantsSimple.DebugPrint();
            Console.WriteLine("-----Sequence Variants Complex");
            this.SequenceVariantsComplex.DebugPrint();
            Console.WriteLine("-----Modified Residues");
            this.ModifiedResiduesUniProt.DebugPrint();
            Console.WriteLine("-----Half cystines");
            this.ModifiedResiduesDisulfide.DebugPrint();
            Console.WriteLine("-----Disulfides");
            this.Disulfides.DebugPrint();

        }
    }
    class MolecularProcessing
    {
        public string CV;
        public string Type;
        public int Begin = 0;
        public int End = 0;
        public MolecularProcessing Next = null;

        public void DebugPrint()
        {
            MolecularProcessing MPRunner = this.Next;
            while (MPRunner != null)
            {
                Console.WriteLine("MolecularProcessing\t" + MPRunner.CV + "\t" + MPRunner.Type + "\t" + MPRunner.Begin + "\t" + MPRunner.End);
                MPRunner = MPRunner.Next;
            }
        }
    }

    class ModifiedResidueUniProt
    {
        public PTMList Modification;
        public int Position = 0;
        public ModifiedResidueUniProt Next = null;

        // Absorbs every node from otherHead's list into this list, then orders the combined
        // list by sequence Position. Existing nodes are MOVED (their Next pointers re-linked),
        // never cloned, so LinkedModifications references into the disulfide list stay valid.
        // otherHead is emptied. (Implements the former "merge half cystines" TODO.)
        public void AbsorbAndSortFrom(ModifiedResidueUniProt otherHead)
        {
            if (ReferenceEquals(this, otherHead)) return;   // a list can't absorb itself; guard against data loss
            var nodes = new List<(ModifiedResidueUniProt Node, int SourceRank)>();
            for (var n = this.Next; n != null; )
            {
                var next = n.Next;
                n.Next = null;
                nodes.Add((n, 0));
                n = next;
            }
            for (var n = otherHead.Next; n != null; )
            {
                var next = n.Next;
                n.Next = null;
                nodes.Add((n, 1));
                n = next;
            }

            var sorted = nodes
                .OrderBy(x => x.Node.Position == 0)   // unknown (0) positions sort LAST
                .ThenBy(x => x.Node.Position)
                .ThenBy(x => x.SourceRank);           // modified residue before half cystine on ties
                                                      // (OrderBy/ThenBy are stable -> original order kept within a group)

            var runner = this;
            foreach (var item in sorted)
            {
                runner.Next = item.Node;
                runner = runner.Next;
            }
            runner.Next = null;
            otherHead.Next = null;                    // the disulfide list has been consumed
        }

        public void DebugPrint()
        {
            ModifiedResidueUniProt MRURunner = this.Next;
            while (MRURunner != null)
            {
                Console.WriteLine("ModifiedResidueUniProt\t" + MRURunner.Modification.Accession + "\t" + MRURunner.Modification.ID + "\t" + MRURunner.Position);
                MRURunner = MRURunner.Next;
            }
        }
    }
    class LinkedModifications
    {
        public ModifiedResidueUniProt Residue1;
        public ModifiedResidueUniProt Residue2;
        public LinkedModifications Next = null;
        public void DebugPrint()
        {
            LinkedModifications LMRunner = this.Next;
            while (LMRunner != null)
            {
                Console.WriteLine("LinkedModifications\t" + LMRunner.Residue1.Modification.ID + "\t" + LMRunner.Residue1.Position + "\t" + LMRunner.Residue2.Modification.ID + "\t" + LMRunner.Residue2.Position);
                LMRunner = LMRunner.Next;
            }
        }
    }
    
    class SequenceVariantSimple
    {
        public char NewAA;
        public int Position;
        public SequenceVariantSimple Next = null;
        public void DebugPrint()
        {
            SequenceVariantSimple SVSRunner = this.Next;
            while (SVSRunner != null)
            {
                Console.WriteLine("SequenceVariantSimple\t" + SVSRunner.NewAA + "\t" + SVSRunner.Position);
                SVSRunner = SVSRunner.Next;
            }
        }
    }

    class SequenceVariantComplex
    {
        public string NewSeq;
        public int Begin = 0;
        public int End = 0;
        public SequenceVariantComplex Next = null;

        public void DebugPrint()
        {
            SequenceVariantComplex SVCRunner = this.Next;
            while (SVCRunner != null)
            {
                Console.WriteLine("SequenceVariantComplex\t" + SVCRunner.NewSeq + "\t" + SVCRunner.Begin + "\t" + SVCRunner.End);
                SVCRunner = SVCRunner.Next;
            }
        }
    }

    class PTMList
    {
        public string  ID;
        public string  Accession;
        public string  PSIModAccession;
        public string  UnimodAccession;
        public PTMList Next;

        // A description-only "modification" with no CV accession, used for UniProt feature types
        // (glycosylation, cross-link, unmatched lipidation) that map to a generic \ModRes.
        public static PTMList Synthetic(string Description)
            => new PTMList { ID = Description, Accession = "", PSIModAccession = null, UnimodAccession = null };

        public PTMList Find(string Description)
        {
            PTMList PTMRunner = this.Next;
            while (PTMRunner != null)
            {
                if (PTMRunner.ID == Description)
                {
                    return PTMRunner;
                }
                PTMRunner = PTMRunner.Next;
            }
            return null;
        }
        public void DebugPrint()
        {
            PTMList PTMListRunner = this.Next;
            while (PTMListRunner != null)
            {
                Console.WriteLine("PTMList\t" + PTMListRunner.ID + "\t" + PTMListRunner.Accession + "\t" + PTMListRunner.PSIModAccession + "\t" + PTMListRunner.UnimodAccession);
                PTMListRunner = PTMListRunner.Next;
            }
        }
        public void ReadPTMListFromFile()
        {
            try
            {
                using (FileStream fs = File.Open("ptmlist.txt", FileMode.Open))
                using (BufferedStream bs = new BufferedStream(fs))
                using (StreamReader sr = new StreamReader(bs))
                {
                    PTMList PTMRunner = this;
                    // Add a PTM missing from ptmlist.txt:
                    PTMRunner.Next = new PTMList();
                    PTMRunner = PTMRunner.Next;
                    PTMRunner.ID = "Half cystine";
                    PTMRunner.Accession = "";
                    PTMRunner.PSIModAccession = "MOD:00798";
                    PTMRunner.UnimodAccession = "374";
                    string ThisLine = sr.ReadLine();
                    string RestOfLine;
                    while (ThisLine !=null)
                    {
                        if (ThisLine.Length > 2)
                        {
                            switch (ThisLine.Substring(0, 2))
                            {
                                case "ID":
                                    PTMRunner.Next = new PTMList();
                                    PTMRunner = PTMRunner.Next;
                                    PTMRunner.ID = ThisLine.Substring(5);
                                    break;
                                case "AC":
                                    PTMRunner.Accession = ThisLine.Substring(5);
                                    break;
                                case "DR":
                                    RestOfLine = ThisLine.Substring(5);
                                    if (RestOfLine.StartsWith("PSI-MOD;"))
                                    {
                                        PTMRunner.PSIModAccession = RestOfLine.Substring(9,RestOfLine.Length-10);
                                    }
                                    else if (RestOfLine.StartsWith("Unimod;"))
                                    {
                                        PTMRunner.UnimodAccession = RestOfLine.Substring(8, RestOfLine.Length-9);
                                    }
                                    break;
                            }

                        }
                        ThisLine = sr.ReadLine();
                    }
                }
            }
            catch (IOException e)
            {
                Console.WriteLine("Something went awry when I tried to read ptmlist.txt");
                System.Environment.Exit(1);
            }
        }
    }
    // Maps a CV accession (e.g. "MOD:00046" or "UNIMOD:34") to its OBO "name:" field, loaded from
    // a standard .obo file in the working directory. PEFF requires the OBO name for ModResPsi /
    // ModResUnimod (not the UniProt synonym). Absence-tolerant: a missing file leaves the map empty
    // and the writer falls back to the UniProt ptmlist ID.
    class OboNameMap
    {
        readonly Dictionary<string, string> Names = new Dictionary<string, string>();
        public int Count => Names.Count;
        public string Find(string Accession)
            => (Accession != null && Names.TryGetValue(Accession, out var Name)) ? Name : null;
        public void ReadFromFile(string Filename)
        {
            try
            {
                using var Reader = new StreamReader(Filename);
                string Id = null, Line;
                bool InTerm = false;   // only ingest [Term] stanzas, never [Typedef]/[Instance]
                while ((Line = Reader.ReadLine()) != null)
                {
                    if (Line.StartsWith("[")) { InTerm = Line.StartsWith("[Term]"); Id = null; }
                    else if (InTerm && Line.StartsWith("id: ")) Id = Line.Substring(4).Trim();
                    else if (InTerm && Line.StartsWith("name: ") && Id != null)
                    {
                        // Obsolete terms are kept on purpose: for an accession UniProt references,
                        // the OBO name: (obsolete or not) is the value PEFF requires; dropping it
                        // would force a UniProt-synonym fallback, the very thing we are fixing.
                        if (!Names.ContainsKey(Id)) Names[Id] = Line.Substring(6).Trim();
                        Id = null;
                    }
                }
            }
            catch (IOException) { /* absent file -> empty map -> writer falls back to UniProt IDs */ }
        }
    }
    class PeffHeader
    {
        public string DbName = "UniProtKB";
        public string Prefix = "sp";
        public string DbDescription = "";
        public string DbSource = "https://www.uniprot.org";
        public string DbVersion = "unknown";   // mandatory key; override with -dbversion
        public string DbDate = "";
        public int NumberOfEntries = 0;
        public bool HasAnnotationIdentifiers = false;
    }
    static class Program
    {
        static void Main(string[] args)
        {
            const string Version = "20260413 alpha";
            // Use periods to separate decimals
            CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
            Console.WriteLine("UniPEFF: Generate PEFF for supplied UniProtKB dtabase");
            Console.WriteLine("David L. Tabb, University Medical Center of Groningen");
            Console.WriteLine("Version " + Version);
            var RecordMolecularProcessing = true;
            var RecordAminoAcidModifications = true;
            var RecordSequenceVariations = true;
            var NextIsInput = false;
            var NextIsOutput = false;
            var NextIsPrefix = false;
            var NextIsDbVersion = false;
            var InputFile = "NA";
            string OutputFile = null;
            string UserPrefix = null;
            string UserDbVersion = null;
            var AnnotationIdentifiers = false;
            var PTMCV = new PTMList();
            var PsiModNames = new OboNameMap();
            var UnimodNames = new OboNameMap();
            foreach (var item in args)
            {
                if (NextIsInput)
                {
                    InputFile = item;
                    NextIsInput = false;
                }
                else if (NextIsOutput)
                {
                    OutputFile = item;
                    NextIsOutput = false;
                }
                else if (NextIsPrefix)
                {
                    UserPrefix = item;
                    NextIsPrefix = false;
                }
                else if (NextIsDbVersion)
                {
                    UserDbVersion = item;
                    NextIsDbVersion = false;
                }
                else switch (item)
                    {
                        case "-OmitMolecularProcessing":
                            RecordMolecularProcessing = false;
                            break;
                        case "-OmitAminoAcidModifications":
                            RecordAminoAcidModifications = false;
                            break;
                        case "-OmitSequenceVariations":
                            RecordSequenceVariations = false;
                            break;
                        case "-in":
                            NextIsInput = true;
                            break;
                        case "-out":
                            NextIsOutput = true;
                            break;
                        case "-prefix":
                            NextIsPrefix = true;
                            break;
                        case "-dbversion":
                            NextIsDbVersion = true;
                            break;
                        case "-AnnotationIdentifiers":
                            AnnotationIdentifiers = true;
                            break;
                        default:
                            Console.Error.WriteLine("\tError: I don't understand this argument: {0}.", item);
                            break;
                    }
            }
            if (NextIsInput || NextIsOutput || NextIsPrefix || NextIsDbVersion)
            {
                Console.Error.WriteLine("\tError: -in, -out, -prefix, and -dbversion each require a following value.");
                return;
            }
            if (RecordAminoAcidModifications)
            {
                Console.WriteLine("Reading ptmlist.txt from current directory.");
                PTMCV.ReadPTMListFromFile();
                PTMCV.DebugPrint();
                PsiModNames.ReadFromFile("psi-mod.obo");
                UnimodNames.ReadFromFile("unimod.obo");
                Console.WriteLine("OBO names loaded from current directory: psi-mod " + PsiModNames.Count + ", unimod " + UnimodNames.Count + " (absent files tolerated; names then fall back to UniProt IDs).");
            }
            Console.WriteLine("Input file is " + InputFile);
            PEFFmodel UniProtDB = null;
            XmlReader XReader = null;
            if (InputFile.EndsWith(".xml.gz"))
                {
                Console.WriteLine("Starting from gzipped XML");
                using FileStream gzFile = File.Open(InputFile, FileMode.Open);
                using var decompressor = new GZipStream(gzFile, CompressionMode.Decompress);
                XReader = XmlReader.Create(decompressor);
                UniProtDB = PEFFmodel.FromUniProtXML(XReader, RecordMolecularProcessing, RecordAminoAcidModifications, RecordSequenceVariations, PTMCV);
            }
            else if (InputFile.EndsWith(".xml"))
                {
                Console.WriteLine("Starting from uncompressed XML");
                using FileStream xmlFile = File.Open(InputFile, FileMode.Open);
                XReader = XmlReader.Create(xmlFile);
                UniProtDB = PEFFmodel.FromUniProtXML(XReader, RecordMolecularProcessing, RecordAminoAcidModifications, RecordSequenceVariations, PTMCV);
            }
            if (UniProtDB != null)
            {
                if (OutputFile == null)
                {
                    if (InputFile.EndsWith(".xml.gz")) OutputFile = InputFile.Substring(0, InputFile.Length - 7) + ".peff";
                    else if (InputFile.EndsWith(".xml")) OutputFile = InputFile.Substring(0, InputFile.Length - 4) + ".peff";
                    else OutputFile = InputFile + ".peff";
                }
                Console.WriteLine("Writing PEFF to " + OutputFile + (AnnotationIdentifiers ? " (with annotation identifiers)" : ""));
                using var PeffWriter = new StreamWriter(OutputFile);
                var Header = new PeffHeader();
                if (!string.IsNullOrEmpty(UserDbVersion)) Header.DbVersion = UserDbVersion;
                UniProtDB.WritePeff(PeffWriter, Header, UserPrefix, AnnotationIdentifiers, PsiModNames, UnimodNames);
            }
        }
    }
}