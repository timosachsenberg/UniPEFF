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
            if (Description.StartsWith("(Microbial infection)")) Description = Description.Substring(22);
            int SemiPosition = Description.IndexOf(';');
            if (SemiPosition > -1) Description = Description.Substring(0, SemiPosition);
            return Description;
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
                        PERunner.DataSet = InputStream["dataset"];
                    }
                    else if (InputStream.Name == "accession")
                    {
                        // Many accessions are possible for a UniProtKB entry; keep just the first
                        if (PERunner.Accession == null)
                        {
                            if (InputStream.Read())
                            {
                                PERunner.Accession = InputStream.Value;
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
                        if (InputStream.Read())
                        {
                            if (InputStream.Read())
                            {
                                if (InputStream["type"] == "primary")
                                {
                                    if (InputStream.Read())
                                    {
                                        PERunner.PrimaryGene = InputStream.Value;
                                    }
                                }
                            }
                        }
                    }
                    else if ((InputStream.Name == "sequence") && (InputStream["length"] != null))
                    {
                        if (InputStream.Read()) PERunner.BaseSequence = InputStream.Value;
                    }
                    if (InputStream.Name == "feature")
                    {
                        if (RecordMolecularProcessing)
                        {
                            switch (InputStream["type"])
                            {
                                case "chain":
                                    MPRunner.Next = new MolecularProcessing();
                                    MPRunner = MPRunner.Next;
                                    MPRunner.CV = "PEFF:0001020";
                                    MPRunner.Type = "mature protein";
                                    // Get begin and end positions from location element afterwards
                                    // Skip ahead until we're on a begin element (within a location).
                                    // As written, this code will die nastily if it hits end of file before it hits a begin element.
                                    while (InputStream.Name != "begin")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.Begin = Int32.Parse(InputStream["position"]);
                                    while (InputStream.Name != "end")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.End = Int32.Parse(InputStream["position"]);
                                    break;
                                case "initiator methionine":
                                    MPRunner.Next = new MolecularProcessing();
                                    MPRunner = MPRunner.Next;
                                    MPRunner.CV = "PEFF:0001035";
                                    MPRunner.Type = "initiator methionine";
                                    // As written, this code will die nastily if it hits end of file before it hits a begin element.
                                    while (InputStream.Name != "position")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                    {
                                        MPRunner.Begin = Int32.Parse(InputStream["position"]);
                                        MPRunner.End = MPRunner.Begin;
                                    }
                                    break;
                                case "propeptide":
                                    // See G5ECN9 for example of both propeptide and signal peptide
                                    MPRunner.Next = new MolecularProcessing();
                                    MPRunner = MPRunner.Next;
                                    MPRunner.CV = "PEFF:0001034";
                                    MPRunner.Type = "propeptide";
                                    // Get begin and end positions from location element afterwards
                                    // As written, this code will die nastily if it hits end of file before it hits a begin element.
                                    while (InputStream.Name != "begin")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.Begin = Int32.Parse(InputStream["position"]);
                                    while (InputStream.Name != "end")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.End = Int32.Parse(InputStream["position"]);
                                    break;
                                case "signal peptide":
                                    MPRunner.Next = new MolecularProcessing();
                                    MPRunner = MPRunner.Next;
                                    MPRunner.CV = "PEFF:0001021";
                                    MPRunner.Type = "signal peptide";
                                    // Get begin and end positions from location element afterwards
                                    // As written, this code will die nastily if it hits end of file before it hits a begin element.
                                    while (InputStream.Name != "begin")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.Begin = Int32.Parse(InputStream["position"]);
                                    while (InputStream.Name != "end")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.End = Int32.Parse(InputStream["position"]);
                                    break;
                                case "transit peptide":
                                    MPRunner.Next = new MolecularProcessing();
                                    MPRunner = MPRunner.Next;
                                    MPRunner.CV = "PEFF:0001022";
                                    MPRunner.Type = "transit peptide";
                                    // Get begin and end positions from location element afterwards
                                    // As written, this code will die nastily if it hits end of file before it hits a begin element.
                                    while (InputStream.Name != "begin")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.Begin = Int32.Parse(InputStream["position"]);
                                    while (InputStream.Name != "end")
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream["position"] != null)
                                        MPRunner.End = Int32.Parse(InputStream["position"]);
                                    break;
                            }
                        }
                        if (RecordAminoAcidModifications)
                        {
                            switch (InputStream["type"])
                            {
                                case "cross-link":
                                    // Isopeptide cross-link (e.g. ubiquitin/SUMO). A single <position>
                                    // is the common inter-chain case; <begin>/<end> is intra-chain. PEFF
                                    // has no generic cross-link-bond key, so we annotate the modified
                                    // residue(s) by position and do not pair them. Route to \ModResPsi/
                                    // \ModResUnimod when the description is a known ptmlist PTM, else \ModRes.
                                    {
                                        string Desc = CleanPtmDescription(InputStream["description"]);
                                        if (string.IsNullOrEmpty(Desc)) Desc = "cross-link";
                                        PTMList XLink = PTMCV.Find(Desc) ?? PTMList.Synthetic(Desc);
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        if (InputStream.Name == "position")
                                        {
                                            MRURunner.Next = new ModifiedResidueUniProt();
                                            MRURunner = MRURunner.Next;
                                            MRURunner.Position = Int32.Parse(InputStream["position"]);
                                            MRURunner.Modification = XLink;
                                        }
                                        else
                                        {
                                            MRURunner.Next = new ModifiedResidueUniProt();
                                            MRURunner = MRURunner.Next;
                                            try { MRURunner.Position = Int32.Parse(InputStream["position"]); }
                                            catch (ArgumentNullException) { MRURunner.Position = 0; }
                                            MRURunner.Modification = XLink;
                                            InputStream.Read();
                                            InputStream.Read();
                                            MRURunner.Next = new ModifiedResidueUniProt();
                                            MRURunner = MRURunner.Next;
                                            try { MRURunner.Position = Int32.Parse(InputStream["position"]); }
                                            catch (ArgumentNullException) { MRURunner.Position = 0; }
                                            MRURunner.Modification = XLink;
                                        }
                                    }
                                    break;
                                case "disulfide bond":
                                    // If a "begin" and "end" are supplied, the disulfide links two sites in the same chain.
                                    // If a "position" is supplied, the disulfide links to another chain (perhaps even another polypeptide?).
                                    // Note that we fake the entry for this mod since it isn't listed in ptmlist.txt
                                    PTMList HalfCystine = PTMCV.Find("Half cystine");
                                    InputStream.Read();
                                    InputStream.Read();
                                    InputStream.Read();
                                    InputStream.Read();
                                    // Console.WriteLine(InputStream.Value);
                                    if (InputStream.Name=="position")
                                    {
                                        // We only have one end of the bond here; don't create a LinkedModifications object.
                                        DSRunner.Next = new ModifiedResidueUniProt();
                                        DSRunner = DSRunner.Next;
                                        DSRunner.Position = Int32.Parse(InputStream["position"]);
                                        DSRunner.Modification = HalfCystine;
                                    }
                                    else
                                    {
                                        //InputStream.Name assumed to be "begin" indicating we have both ends of the bond
                                        DSRunner.Next = new ModifiedResidueUniProt();
                                        DSRunner = DSRunner.Next;
                                        try
                                        {
                                            DSRunner.Position = Int32.Parse(InputStream["position"]);
                                        }
                                        catch (ArgumentNullException failure)
                                        {
                                            // Sometimes UniProt throws us an "unknown"
                                            DSRunner.Position = 0;
                                        }
                                        DSRunner.Modification = HalfCystine;
                                        ModifiedResidueUniProt FirstPosition = DSRunner;
                                        InputStream.Read();
                                        InputStream.Read();
                                        DSRunner.Next = new ModifiedResidueUniProt();
                                        DSRunner = DSRunner.Next;
                                        try
                                        {
                                            DSRunner.Position = Int32.Parse(InputStream["position"]);
                                        }
                                        catch (ArgumentNullException failure)
                                        {
                                            DSRunner.Position = 0;
                                        }
                                        DSRunner.Modification = HalfCystine;
                                        ModifiedResidueUniProt SecondPosition = DSRunner; 
                                        DSLMRunner.Next = new LinkedModifications();
                                        DSLMRunner = DSLMRunner.Next;
                                        DSLMRunner.Residue1 = FirstPosition;
                                        DSLMRunner.Residue2 = SecondPosition;
                                    }
                                    break;
                                case "glycosylation site":
                                    // Glycan compositions are not in PSI-MOD/Unimod, so emit a generic
                                    // \ModRes carrying the UniProt description (e.g. "N-linked (GlcNAc...)").
                                    {
                                        string Desc = CleanPtmDescription(InputStream["description"]);
                                        if (string.IsNullOrEmpty(Desc)) Desc = "glycosylation site";
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        MRURunner.Next = new ModifiedResidueUniProt();
                                        MRURunner = MRURunner.Next;
                                        MRURunner.Position = Int32.Parse(InputStream["position"]);
                                        MRURunner.Modification = PTMList.Synthetic(Desc);
                                    }
                                    break;
                                case "lipid moiety-binding region":
                                    // Lipidation: route to \ModResPsi/\ModResUnimod if the description is a
                                    // known ptmlist PTM, else a generic \ModRes.
                                    {
                                        string Desc = CleanPtmDescription(InputStream["description"]);
                                        if (string.IsNullOrEmpty(Desc)) Desc = "lipid moiety-binding region";
                                        PTMList Hit = PTMCV.Find(Desc);
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        MRURunner.Next = new ModifiedResidueUniProt();
                                        MRURunner = MRURunner.Next;
                                        MRURunner.Position = Int32.Parse(InputStream["position"]);
                                        MRURunner.Modification = Hit ?? PTMList.Synthetic(Desc);
                                    }
                                    break;
                                case "modified residue":
                                    //  <feature type="modified residue" description="N-acetylserine" evidence="6">
                                    // <feature type="modified residue" description="(Microbial infection) O-acetylthreonine; by Yersinia YopJ; alternate" evidence="31">
                                    string Description = CleanPtmDescription(InputStream["description"]);
                                    PTMList CVHit = PTMCV.Find(Description);
                                    if (CVHit != null)
                                    {
                                        // Advance to the position in sequence where this occurs
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        MRURunner.Next = new ModifiedResidueUniProt();
                                        MRURunner = MRURunner.Next;
                                        MRURunner.Position = Int32.Parse(InputStream["position"]);
                                        MRURunner.Modification = CVHit;
                                    }
                                    else
                                    {
                                        Console.WriteLine("Failed to find this modified residue in ptmlist.txt:\t" + Description);
                                        Console.WriteLine("Please update your copy of ptmlist.txt from here:");
                                        Console.WriteLine("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt");
                                    }
                                    break;
                                case "non-standard amino acid":
                                    break;
                            }
                        }
                        if (RecordSequenceVariations)
                        {
                            switch (InputStream["type"])
                            {
                                case "sequence conflict":
                                    break;
                                case "sequence variant":
                                    // We don't yet know whether this variant is a simple or a complex
                                    // Skip to the contents of the Original or Location element.  In 5640, "VAR_083056" is an example of deletion.
                                    while ((InputStream.Name != "original") && (InputStream.Name != "location"))
                                    {
                                        InputStream.Read();
                                    }
                                    if (InputStream.Name == "location")
                                    {
                                        // This is a deletion
                                        SVCRunner.Next = new SequenceVariantComplex();
                                        SVCRunner = SVCRunner.Next;
                                        SVCRunner.NewSeq = "";
                                        while ((InputStream.Name != "begin") && (InputStream.Name != "position"))
                                        {
                                            InputStream.Read();
                                        }
                                        if (InputStream.Name == "position")
                                        {
                                            // This affects a single amino acid position
                                            SVCRunner.Begin = Int32.Parse(InputStream["position"]);
                                            SVCRunner.End = SVCRunner.Begin;
                                        }
                                        else
                                        {
                                            // This affects a range of amino acid positions
                                            SVCRunner.Begin = Int32.Parse(InputStream["position"]);
                                            while (InputStream.Name != "end")
                                            {
                                                InputStream.Read();
                                            }
                                            SVCRunner.End = Int32.Parse(InputStream["position"]);
                                        }
                                    }
                                    else
                                    {
                                        // This replaces one seq with another
                                        // We are currently at "original"
                                        string Original;
                                        string Variation;
                                        int    Begin=0;
                                        int    End=0;
                                        InputStream.Read();
                                        Original = InputStream.Value;
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        Variation = InputStream.Value;
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        InputStream.Read();
                                        if (InputStream.Name == "position")
                                        {
                                            // This is a single letter change
                                            SVSRunner.Next = new SequenceVariantSimple();
                                            SVSRunner = SVSRunner.Next;
                                            SVSRunner.NewAA = Variation[0];
                                            SVSRunner.Position = Int32.Parse(InputStream["position"]);
                                        }
                                        else
                                        {
                                            // This is a multiple letter change
                                            Begin = Int32.Parse(InputStream["position"]);
                                            InputStream.Read();
                                            InputStream.Read();
                                            End = Int32.Parse(InputStream["position"]);
                                            SVCRunner.Next = new SequenceVariantComplex();
                                            SVCRunner = SVCRunner.Next;
                                            SVCRunner.NewSeq = Variation;
                                            SVCRunner.Begin = Begin;
                                            SVCRunner.End = End;
                                        }
                                    }
                                    break;
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
            int EntryCount = 0;
            for (var ER = Entries.Next; ER != null; ER = ER.Next) EntryCount++;
            Header.NumberOfEntries = EntryCount;
            Header.HasAnnotationIdentifiers = AnnotationIdentifiers;
            // A user-supplied -prefix wins; otherwise derive one (sp/tr) from the first entry.
            // Either way a single prefix is used for the header and every entry line.
            if (!string.IsNullOrEmpty(PrefixOverride)) Header.Prefix = PrefixOverride;
            else if (Entries.Next != null) Header.Prefix = PEFFentry.PrefixForDataset(Entries.Next.DataSet);
            WritePeffHeader(Writer, Header);
            int Fallbacks = 0;
            for (var ER = Entries.Next; ER != null; ER = ER.Next)
            {
                Fallbacks += ER.WritePeffEntry(Writer, Header.Prefix, AnnotationIdentifiers, PsiModNames, UnimodNames);
            }
            if (Fallbacks > 0)
                Console.Error.WriteLine("\tNote: " + Fallbacks + " modification name(s) were not found in the supplied OBO files; emitted the UniProt ptmlist name instead (not strictly the PEFF-required OBO name).");
        }

        static void WritePeffHeader(TextWriter Writer, PeffHeader Header)
        {
            Writer.Write("# PEFF 1.0\n");
            Writer.Write("# GeneralComment=Generated by UniPEFF\n");
            Writer.Write("# //\n");
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
    class PEFFentry
    {
        public string DataSet;
        public string Accession;
        public string Name;
        public string FullName;
        public string PrimaryGene;
        public string BaseSequence;
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
            if (!string.IsNullOrEmpty(FullName))     Tag("PName", "(" + EscapePeff(FullName) + ")");
            if (!string.IsNullOrEmpty(PrimaryGene))  Tag("GName", EscapePeff(PrimaryGene));
            if (!string.IsNullOrEmpty(BaseSequence)) Tag("Length", BaseSequence.Length.ToString());
            // \DbUniqueId is intentionally NOT emitted: the PSI-MS CV (PEFF:0001001) states it
            // "shall not be used in the PEFF 1.0 serialization as it is redundant with the primary
            // identifier following the >".
            if (!string.IsNullOrEmpty(Name))         Tag("ID", EscapePeff(Name));

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
                        else if (OboNames.Count > 0) Fallbacks++;
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

            // Simple (single-residue) variants
            {
                var Builder = new StringBuilder();
                for (var V = SequenceVariantsSimple.Next; V != null; V = V.Next)
                {
                    if (V.Position == 0) continue;   // PEFF VariantSimple position MUST be > 0; "?" is ModRes-only
                    Builder.Append('(').Append(IdPrefix()).Append(V.Position).Append('|').Append(EscapePeff(V.NewAA.ToString())).Append(')');
                }
                if (Builder.Length > 0) Tag("VariantSimple", Builder.ToString());
            }
            // Complex (multi-residue / deletion) variants
            {
                var Builder = new StringBuilder();
                for (var V = SequenceVariantsComplex.Next; V != null; V = V.Next)
                {
                    if (V.Begin == 0 || V.End == 0) continue;   // positions count from 1; "?" is ModRes-only
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