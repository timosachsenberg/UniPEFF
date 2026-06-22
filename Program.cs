using System;
using System.IO;
using System.IO.Compression;
using System.Xml;
using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Text.Json;

namespace UniPEFF
{
    class PEFFmodel
    {
        string Source;
        PEFFentry Entries = new PEFFentry();
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
                                    // two positions!
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
                                    break;
                                case "lipid moiety-binding region":
                                    break;
                                case "modified residue":
                                    //  <feature type="modified residue" description="N-acetylserine" evidence="6">
                                    // <feature type="modified residue" description="(Microbial infection) O-acetylthreonine; by Yersinia YopJ; alternate" evidence="31">
                                    string Description = InputStream["description"];
                                    // Fix a weird mistake in PTM formatting from UniProtKB.
                                    if (Description.StartsWith("(Microbial infection)"))
                                    {
                                        Description = Description.Substring(22);
                                    }
                                    int     SemiPosition = Description.IndexOf(';');
                                    if (SemiPosition > -1)
                                    {
                                        Description = Description.Substring(0, SemiPosition);
                                    }
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
            PERunner.DebugPrint();
            return ThisPEFF;
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

        // TODO: merge sort the half cystines into this list after all the disulfide and modified residue lines are processed
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
            var InputFile = "NA";
            var PTMCV = new PTMList();
            foreach (var item in args)
            {
                if (NextIsInput)
                {
                    InputFile = item;
                    NextIsInput = false;
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
                        default:
                            Console.Error.WriteLine("\tError: I don't understand this argument: {0}.", item);
                            break;
                    }
            }
            if (RecordAminoAcidModifications)
            {
                Console.WriteLine("Reading ptmlist.txt from current directory.");
                PTMCV.ReadPTMListFromFile();
                PTMCV.DebugPrint();
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
        }
    }
}