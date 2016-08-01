#! /usr/bin/python
# -*- coding: utf-8 -*-
# Created by Roberto Preste

import requests, sys, json, csv, os, microdata, urllib
from datetime import date
from bs4 import BeautifulSoup
from random import randint

edamList = ["PCR experiment", "Secondary structure comparison", "GO concept ID (cellular component)", "Enzyme kinetics report format", "Protein post-translation modification site prediction", "Ensembl ID ('Rattus norvegicus')", "Low complexity sequences", "SCOP species", "Vector sequence detection", "Gene ID", "Genotyping", "Sequence tag profile (with gene assignment)", "Protein residue interactions", "Username", "Protein signal peptide detection (eukaryotes)", "CABRI accession", "Restriction map", "Trim vector", "Molecular modelling", "ConsensusPathDB entity ID", "GlobPlot domain image", "Microarray data processing", "EMBL/GenBank/DDBJ ID", "Number of output entities", "Multiple structure alignment construction", "Experiment report (genotyping)", "Sequence database search (by sequence using global alignment-based methods)", "Anonymisation", "MonosaccharideDB ID", "Sequence offset", "Protein repeat signature", "REDIdb ID", "Job type", "NACCESS log file", "DNA linear map rendering", "Drugs and target structures", "Gene ID (VectorBase)", "Sequence variation annotation format", "Microarray proximity map plotting", "Protein name (UniProt)", "Evolutionary biology", "Microbial collection", "SCF", "BCML", "markx3", "EMBOSS wordfinder log file", "Homology-based gene prediction", "Phylogenetic discrete states format", "ArrayExpress accession number", "ChIP-seq", "RNA sequence (raw)", "Sequence set (nucleic acid)", "Freshwater biology", "Sequence cluster processing", "Superfamily hidden Markov model number", "Restriction site recognition", "Position-specific scoring matrix", "Protein family signature", "Enzyme ID", "HMMER emission and transition", "ebwtl", "Taxonomy", "PDB residue number", "Score", "Differential protein expression analysis", "Phylogenomics", "Amino acid comparison matrix (floats)", "Mitochondria", "Chromosome annotation (aberration)", "Annotation processing", "Structure database search", "Mathematical modelling", "Gene symbol annotation", "Plasmid map", "ASTD ID (exon)", "rast", "EMBOSS megamerger log file", "Biological pathway or network report format", "Nucleic acid features (quadruplexes)", "SBRML", "GermOnline ID", "Gene ID (Xenbase)", "Nutritional science", "Sequence alignment (protein)", "Cardiology", "Alignment processing", "Gene structure", "HMMER NULL hidden Markov model", "Selenocysteine insertion sequence (SECIS) prediction", "Gene name (HUGO)", "Raw SCOP domain classification", "findkm", "mmCIF", "Matrix/scaffold attachment sites", "CATH domain sequences (COMBS)", "Cell line name", "Sequence version", "Transcriptomics", "Protein secondary structure processing", "Introns", "Secondary structure data", "consensus", "Human disease", "ebwt", "BLAST results", "Data retrieval (database metadata)", "Cell migration analysis", "Protein cysteine and disulfide bond assignment", "Genome assembly", "Structural data", "Sequence database hits evaluation data", "Nucleic acid features report (splice sites)", "SPSS", "Manchester OWL Syntax", "Ontology name", "Annotation track", "DNA methylation", "Protein domain classification format", "Sequence cluster ID (UniRef50)", "Tanimoto similarity score", "Immunogen design", "PSI MI TAB (MITAB)", "Nucleic acid sequence comparison", "Splitting", "CpG island and isochores", "match", "Recombination detection", "Light microscopy", "Nucleic acid sequence record", "Nucleic acid features report (promoters)", "Gene finding", "Quadruplex formation site detection", "Cell type name", "Variant pattern analysis", "Genome indexing", "Data retrieval (gene annotation)", "3D-1D scoring matrix format", "genePred", "Blattner number", "AB1", "Gene name (Arabidopsis)", "Image", "Structure alignment (protein C-alpha atoms)", "Multiple protein tertiary structure alignment (C-alpha atoms)", "Ontology comparison", "CIP strain data format", "Sequence profile data", "Molecule accession", "Ontology mapping", "Backbone torsion angle calculation", "Gene tree", "2D PAGE gel report", "Misspelling", "Disease pathways", "Sequence database search (by sequence using local alignment-based methods)", "Annotated URI", "File base name", "completely unambiguous pure nucleotide", "MEME background frequencies file", "Protein secondary structure prediction (disulfide bonds)", "Protein features report (super-secondary)", "Secondary structure alignment (RNA)", "Database name", "Sequence-to-3D-profile alignment", "Protein property comparison", "Gene expression", "Structure editing", "Job identifier", "Phylogenetic tree analysis", "Organism accession", "Phylogenetic continuous quantitative data", "Structure determination", "Virus identifier", "selex sequence format", "Small molecule report", "nbrf/pir", "Genome comparison", "Simulation experiment report", "Nucleic acid sequence", "DaliLite hit table", "Carbohydrate accession", "Chromatogram visualisation", "Frameshift detection", "Sequence cluster ID (UniRef)", "Taxonomic classification", "Locus ID (AGI)", "Protein crystallizability prediction", "Cell line report", "WormBase gene report format", "Database entry version information", "DragonDB gene report format", "Genome annotation", "Microarray experiment", "SAM", "Phylogenetic tree format", "Synonym", "Author ID", "DragonDB author identifier", "Primer and probe design", "MEME motif alphabet", "RNA secondary structure", "RNA", "EMBL format (XML)", "LSID", "Microarray scatter plot plotting", "Gene functional annotation", "RNAML", "Drug accession", "pdbseqres", "ig", "Multiple protein tertiary structure alignment (all atoms)", "Structure comparison", "Protein structural motif", "Protein-protein interaction prediction (from protein structure)", "RNAi report", "Gene component prediction", "GPML", "Pathway ID (ConsensusPathDB)", "DNA back-translation", "MaizeGDB gene report format", "de Novo sequencing", "Structural transformation matrix", "Tree-based sequence alignment", "Gene ID (FlyBase)", "Sequence database hits alignments list", "RNA structure prediction", "Amino acid index format", "Structure alignment (protein)", "RNA features report", "Sequencing", "Nucleic acid probability profile", "debug-feat", "Sequence accession", "UniGene entry format", "Locus ID (DictyBase)", "GO concept name (biological process)", "Structure alignment (protein)", "Protein microarrays", "SMILES", "Immunoprecipitation experiment", "pure rna", "Protein model validation", "Structure alignment", "Structure database search results", "Microarray tree or dendrogram rendering", "Core data", "Bibliography", "Promoter prediction", "Phylogenetic tree", "Statistical modelling", "Protein sequence alignment", "Protein residue surface calculation", "URL", "Genetic variation analysis", "Sequence-profile alignment (Domainatrix signature)", "Ontology term", "Sequence checksum generation", "Electron microscopy", "xpm", "customtrack", "Map accession", "Gene3D entry format", "Protein signal peptide detection", "spML", "Sequence assembly report", "Image metadata", "RNA secondary structure analysis", "DNA substitution modelling", "Ontologies, nomenclature and classification", "Standardization and normalization", "Codon usage table processing", "Amino acid pair-wise contact potentials", "Nucleic acid features (primers) format", "Immunoproteins, genes and antigens", "Hidden Markov model format", "HMM emission and transition counts", "Database cross-mapping", "Nucleic acid density plot", "Nucleic acid temperature profile plotting", "HIVDB entry format", "Structured RNA prediction and optimisation", "Structure identifier", "Protein structure raw data", "Biological model ID", "Data retrieval (ontology concept)", "Protein secondary structure assignment", "Chemical name (INN)", "KEGG DRUG entry format", "Nucleic acid comparison", "Microarray protocol annotation", "Regulatory element prediction", "EMAP concept ID", "MIRIAM data type synonymous name", "Biological model name", "Gene identifier (NCBI RefSeq)", "Structure visualisation", "Structure processing (protein)", "Data retrieval (feature table)", "Feature table processing", "Sequence-MEME profile alignment", "Protein non-covalent interactions report", "Molecular replacement", "Chemical registry number", "Classification", "Protein-ligand interaction report", "Data resource definition name", "Literature data resources", "Stockholm format", "Ensembl gene ID", "Tool metadata", "Amino acid index (White-Wimley data)", "genpept", "Medical imaging", "Sequence database search", "Protein expression", "Genome assembly", "Window size", "Drug ID (TTD)", "Sequence accession (nucleic acid)", "Sequence assembly ID", "Sequence feature key", "Exactly 2", "Differential expression analysis", "Protein sites and features", "Gene name (KEGG GENES)", "Biodiversity", "Neurobiology", "Motif database search", "Sequence length range", "Gene Expression Atlas Experiment ID", "Structure retrieval", "dbSNP ID", "Nucleic acid classification", "Protein-nucleic acid binding site analysis", "Phylogenetic tree image", "EMBOSS sequence type", "CATH identifier", "IPI protein ID", "Protein-protein interaction prediction (from protein sequence)", "Data retrieval (genotype and phenotype annotation)", "Chemical structure image", "PDF", "Sequencing quality control", "Microarray experiment report", "Microarray experiment data format", "EMBL accession", "qualsolid", "Consensus-based sequence alignment", "Phylogenetic tree generation (parsimony methods)", "dhf", "Ensembl ID ('Macaca mulatta')", "Ensembl ID ('Ornithorhynchus anatinus')", "CATH node ID (family)", "Tag-based peptide identification", "SQLite", "Molecular interaction ID", "Nucleic acid structure comparison", "Protein databases", "Ensembl ID ('Oryctolagus cuniculus')", "GCDML", "Gene expression data format", "UniProt ID", "NeuroMorpho ID", "Diffraction data reduction", "Protein domain recognition", "Concentration", "Sequence record format", "Monosaccharide identifier", "Protein dipole moment", "Reproductive health", "Protein-metal contact calculation", "CATH topology", "Enzyme ID (BioCyc)", "Sequence features (comparative)", "SGD ID", "Gene identifier", "LocARNA PP", "Toxins and targets", "Protein modelling (side chains)", "Mass spectrometry data format", "Gene ID (HGNC)", "FlyBase primary identifier", "PS", "Zinc finger prediction", "Protein ID (CORUM)", "Sequence record format (XML)", "Regression analysis", "Protein structural quality report", "Codon usage table generation", "PeptideAtlas ID", "Mobile genetic element ID", "ChEBI ID", "unambiguous pure dna", "Sequence feature name", "Pathway ID (Panther)", "Ear, nose and throat medicine", "Amino acid name (single letter)", "Ramachandran plot", "GeneIlluminator gene report format", "Gene expression matrix", "Data integration and warehousing", "Protein hydrophobic moment plotting", "Hydrogen bond calculation", "Nucleotide identifier", "Michaelis Menten plot", "Bisulfite mapping", "Protein extinction coefficient", "Validation", "UniProt-like (text)", "CATH architecture", "CATH functional category", "Reaction ID (BioCyc)", "Protein titration curve", "Molecular property identifier", "Protein feature identifier", "XML", "Drug ID (KEGG)", "ACLAME ID", "DNA translation", "Sequence features", "TCDB ID", "Sequence database search (by property)", "Codon name", "Neurology", "Protein domains", "Document clustering", "Community profiling", "Codon usage table ID", "Simulated gene expression data generation", "EC number", "ppm", "Phylip tree distance format", "dbEST accession", "Protein-ligand interactions", "Protein subcellular localization prediction", "Article ID", "Sequence alignment visualisation", "TREMBL accession", "Transposon prediction", "SCOP protein", "Format detection", "Disease report", "SBS experimental data", "SGD gene report format", "Sequence variation ID", "Protein distance matrix calculation", "Sequence processing (nucleic acid)", "Nucleic acid melting temperature", "Ensembl ID ('Pan troglodytes')", "Probes and primers", "DNA variation", "Genome map", "Nucleic acid features report (primers)", "Sequence editing (nucleic acid)", "Sequence cutting", "SMART protein schematic", "Spot ID (HSC-2DPAGE)", "Gramene primary identifier", "Ensembl ID ('Danio rerio')", "Protein modification ID", "Nucleic acid features (restriction sites) format", "Base pairing probability matrix dotplot", "DNA structure prediction", "PCR experiment report", "Protein property calculation (from structure)", "MHC peptide immunogenicity report", "STRIDE log file", "Ecology", "Structural profile alignment generation (multiple)", "Physical map", "Ensembl ID ('Gallus gallus')", "GWAS report", "IDAT", "Locus annotation", "qcML", "Sequence identifier (nucleic acid)", "Sequence position", "Workflow format", "GPCR prediction", "mase format", "Article analysis", "phylipnon sequence format", "Protein SNP mapping", "Pharmacogenomic test report", "nexus alignment format", "PCR primer design", "Protein sequence feature detection", "Ensembl ID ('Myotis lucifugus')", "Clone ID (IMAGE)", "Phylogenetic tree distances calculation", "Isotopic distributions calculation", "Bibliography generation", "Protein surface calculation", "Threading", "Beta diversity data", "EMBL-like (XML)", "Genetic linkage report", "Protein secondary structure assignment (from coordinate data)", "Probabilistic sequence generation", "Codon adaptation index", "Amino acid comparison matrix (integers)", "Workflows", "Protein surface calculation (accessible molecular)", "Mobile genetic elements", "HET group detection", "Sequence assembly component", "Genome indexing (Burrows-Wheeler)", "uniprotkb-like format", "Membrane and lipoproteins", "Multiple sequence alignment", "ASTD ID", "Ordination plot", "Aggregation", "MEME background Markov model", "NeXML", "Gene transcription", "EMBOSS supermatcher error file", "2D PAGE spot report", "FIG ID", "Identifier with metadata", "Gene map", "GEO accession number", "mzXML", "HPA antibody id", "Sequence-profile alignment format", "Gene regulatory network analysis", "EMBL format", "AraC-XylS ID", "srspair", "Protein folding simulation", "CATH domain sequences (ATOM)", "Sequence alignment validation", "Tumor annotation", "Molecular model refinement", "Medicine", "Nucleotide comparison matrix (floats)", "SBOL", "RNA sequence", "Data retrieval (phylogenetic tree)", "protXML", "Gene3D ID", "PHYLIP format", "RNAi experiment", "Document, record and content management", "Protein aliphatic index calculation", "Virus annotation", "Rigid body refinement", "Feature table query", "Chemical formula format", "Microbial ecology", "Gap separation penalty", "Nucleic acid features (microRNA)", "Nucleic acid signature", "RNA-seq time series data analysis", "Nucleosome exclusion sequences", "Sequence database search (by amino acid composition)", "cdsxml", "Sequence profile database search", "Sequence cluster ID (UniRef90)", "Protein secondary structure prediction (integrated)", "ZTR", "im", "Nucleic acid structure comparison", "Pharmacokinetics and pharmacodynamics", "Medical toxicology", "Protein sequence record", "Secondary structure processing", "Structure prediction", "SMART accession number", "Genome build identifier", "Atomic y coordinate", "PATIKA entry format", "Biological model accession", "Base-calling", "COMBINE OMEX", "Protein features report (membrane regions)", "Web portal", "GelML", "Workflow", "PDB residue name", "ArrayExpress entry format", "Ontology concept comment", "Gene expression QTL analysis", "CATH domain report", "Residue packing validation", "Repeat sequence organisation analysis", "Marine biology", "Gene regulatory network prediction", "Protein structural motifs and surfaces", "Mass spectrometry", "Structure retrieval (by keyword)", "Molecular surface analysis", "Nucleic acid features report (restriction sites)", "Nucleic acid thermodynamics", "Microarray raw data analysis", "Chemical class enrichment", "Data architecture, analysis and design", "BioModel mathematical model format", "STRING entry format (HTML)", "Locus ID (UTR)", "CATH chain report format", "Dentistry", "Protein-nucleic acid binding sites", "Pfam accession number", "Hybrid sequence alignment construction", "Sample comparison", "Protein region signature", "Microarray dendrograph plotting", "Pfam clan ID", "EST accession", "Ab-initio gene prediction", "Position frequency matrix", "qual", "Atomic x coordinate", "Sequence features (repeats) format", "Sequence variations", "Job ID", "Thermo RAW", "bed12", "Molecular property (general)", "Databank", "Phylip tree format", "Sequence editing", "Ontology concept definition", "Atom name", "KRSS2 Syntax", "TRANSFAC accession number", "BEL", "Position weight matrix", "mega", "Gene ID (MfunGD)", "Protein sequence alignment analysis", "Phylogenetic tree report (invariants) format", "Phylogenetic tree analysis (natural selection)", "Gene transcriptional features report", "Sequencing error detection", "GO (biological process)", "Residue contact calculation (residue-negative ion)", "Mapping", "MHC Class II epitopes report", "MSDchem ligand dictionary entry format", "Ontology and terminology", "mira", "Haplotype mapping", "Protein folding, stability and design", "SNP calling", "PED", "Protein modelling (loops)", "completely unambiguous", "Linucs ID", "Gene regulatory networks", "EST and cDNA sequence analysis", "Editing", "Pairwise protein tertiary structure alignment (C-alpha atoms)", "Alignment format (XML)", "Coding region prediction", "Rare diseases", "Protein complex", "Nucleic acid features report (PolyA signal or site)", "Proteomics experiment report", "Representative sequence identification", "nii", "Data management", "Ontology concept name", "ambiguous", "MGF", "GPCR coupling selectivity prediction", "FASTA-aln", "prosite-profile", "Pathway ID (NCI-Nature)", "Structure alignment (RNA)", "Gene ID (EcoGene)", "IntEnz enzyme report format", "Sequence composition calculation (protein)", "Raw CATH domain classification format", "GCG", "clustal sequence format", "Nucleic acid melting curve", "nexusnon", "Structure alignment (multiple)", "DNA features", "DNA mapping", "Methylation level analysis (gene-specific)", "Metabolites", "Data security", "Mapping assembly", "Phylogenetic tree generation (minimum distance methods)", "Misnomer", "Visualisation", "Notation3", "Sequence record", "Secondary structure alignment metadata (protein)", "Gene and protein families", "Taxon", "Tertiary structure record", "Nucleic acid melting profile plotting", "Sequence report", "Protein modelling (backbone)", "TIGRFam entry format", "Protein fold recognition", "Protein pH-dependent property calculation", "BIND entry format", "Phylogenetic character data", "Nucleic acid folding family identification", "File name", "rgb", "Population genetics", "Database search", "Protein modifications", "Gene ID (GeneDB Plasmodium falciparum)", "Sequence-profile alignment (fingerprint)", "Sequence alignment report", "DNA base structural data", "Sequence alignment report (site correlation)", "Ensembl variation file format", "EMBOSS database resource definition", "lhf", "Tabix index file format", "Data retrieval (RNA family annotation)", "Sequence cluster visualisation", "Sequence length", "Data reference", "Amino acid index", "Genome accession", "GTF", "Version information", "Translation phase specification", "Amino acid index (chemical classes)", "Microarray image", "Structure alignment (protein all atoms)", "Vertebrates", "WormBase name", "Data retrieval (identifier)", "Metabolic pathways", "Gene ID (GeneFarm)", "Diffraction data analysis", "Protein folding pathway prediction", "Data index data", "Homology modelling", "Data resource definition ID", "RNA-Seq alignment", "Worms", "Sequence composition plot", "Scents", "Experimental data", "SNP annotation", "Microarray metadata", "Hydrogen bond calculation (inter-residue)", "Cross-assembly", "Mass spectrometry data", "Data retrieval (protein interaction annotation)", "MIRIAM data type primary name", "Protein features report (secondary structure)", "MAF", "Listfile processing", "restover format", "XSLT stylesheet", "Medical biotechnology", "CiteXplore-all", "Chemical formula", "Amino acid property", "Genome feature comparison", "PlasMapper TextMap", "Compound name", "iTRAQ", "Retention times prediction", "Database metadata", "DIALIGN format", "Sequence tagged site (STS) mapping", "Amino acid index (hydropathy)", "ATC code", "Splice transcript prediction", "Lipids", "Word size", "Drug ID (PharmGKB)", "Gap separation penalty (integer)", "Protein ID (EMBL/GenBank/DDBJ)", "Sequence accession (protein)", "Protein quaternary structure prediction", "Protein signal peptides", "pdbatom", "Protein tertiary structure", "Format (typed)", "Protein-nucleic acid binding prediction", "CAF", "Reference sample report", "Structure report", "Sequence similarity", "dasdna", "Sequence feature analysis", "Transcriptome assembly (de novo)", "Electron microscopy model ID", "Sequence database search (by sequence using profile-based methods)", "Sequence alignment refinement", "Sequence feature identifier", "CopasiML", "XML Schema", "Protein features report (binding sites)", "Protein topology", "Data index format", "Structure databases", "Loading", "UniProt accession", "axt", "Target-Decoy", "Protein sequence composition", "KEGG DISEASE entry format", "Repeat sequences", "Pain medicine", "Codon usage table formatting", "JSON", "dat", "Molecular charge", "Data index report", "Search parameter", "User metadata", "Nucleic acid sequence visualisation", "Gap extension penalty (float)", "Dotplot", "Terminal gap extension penalty", "Data retrieval (ontology annotation)", "Microarray data standardization and normalization", "Taxonomy", "Protein feature prediction (from structure)", "Protein binding sites", "Read pre-processing", "mzML", "Protein secondary structure assignment (from CD data)", "Sequence identity", "Nucleic acid report", "Sequence mutation and randomization", "Sequence database hits scores list", "Sequence comparison", "Peptide immunogenicity prediction", "Amino acid name (full name)", "Workflow ID", "RNA secondary structure prediction (shape-based)", "Ab initio structure prediction", "Protein ionization curve", "tRNA", "Cell and tissue culture", "GenBank feature", "Protein identifier", "Map data", "TIGR gene report format", "Gene expression analysis", "Biodiversity report", "VNTR", "Linked data format", "ChemSpider ID", "LIPID MAPS ID", "Sequence search", "Functional genomics", "System metadata", "Parameter", "FlyBase gene report format", "Protein surface report", "Phylogenetic tree generation (consensus)", "Pathway ID (aMAZE)", "Pathway ID (CPDB)", "Enzyme kinetics data", "Protein modelling (mutation)", "Codon usage", "Pfam entry format", "Amino acid identifier", "CATH class", "EMBL-like (text)", "CATH structurally similar group", "Residue validation", "Protein hydrogen exchange rate", "Nucleic acid feature identifier", "Processed microarray data", "Phenomenon identifier", "Nucleic acid sequence alignment analysis", "TMT-tag", "Protein geometry validation", "ID retrieval", "Incident curve plotting", "EMBOSS Uniform Feature Object", "KEGG Glycan ID", "RNA secondary structure format", "Compound ID (BioCyc)", "Infectious disease", "prosite-pattern", "Nucleic acid features (siRNA)", "Protein flexibility and motion analysis", "FlyBase secondary identifier", "Ontology concept ID", "Enzyme identifier", "Clinical trial report", "Alignment data", "Ontology concept identifier", "GFF3-seq", "Protein secondary structure prediction (helices)", "WormBase class", "siRNA duplex prediction", "DNA replication and recombination", "pir", "treecon sequence format", "Comparison matrix (amino acid)", "Genomics", "Restriction map drawing", "Chemical name (IUPAC)", "Structure file processing", "mzIdentML", "SRA format", "Chromatographic alignment", "Spot serial number", "Systems medicine", "Ensembl ID ('Ciona savignyi')", "Ensembl ID ('Erinaceus europaeus')", "Sequence alignment comparison", "EMBOSS report", "Gene expression report format", "Signaling Gateway protein ID", "InterPro architecture image", "Vienna RNA calculated energy", "RNA splicing", "VCF", "Genetic organisation", "SVG", "Protein-protein interaction report", "Nucleic acid feature detection", "Ensembl ID ('Felis catus')", "Chemical identifier", "Compound ID (HMDB)", "HMMER profile alignment (sequences versus HMMs)", "Proteome", "Linkage disequilibrium calculation", "bedgraph", "PubMed citation", "Phylip character frequencies format", "Sequence identifier (protein)", "Prediction and recognition (nucleic acid)", "Sequence feature ID", "Protein function prediction", "SNP", "KEGG organism code", "Design", "Humans", "CATH node", "Nucleic acid features report (STS)", "PIR identifier", "UTRdb taxon", "K-mer countgraph", "KEGG ENZYME enzyme report format", "SCOP family", "MPSS experimental data", "Gramene secondary identifier", "Structural genomics target selection", "Alpha diversity data", "Protein report (transcription factor)", "Gene ID (MIPS)", "Genetic code", "Radiation Hybrid (RH) scores", "Data", "Raw sequence format", "Molecular biology", "STRING ID", "Sequence map", "Sequence assembly", "Ensembl ID ('Loxodonta africana')", "Database ID", "ASTD ID (polya)", "Protein hydropathy calculation (from structure)", "DaliLite log file", "Secondary structure alignment", "Toggle", "Analytical chemistry", "Gender medicine", "ASTD ID (tss)", "EMBOSS sites log file", "Gene name (Genolist)", "HAMAP ID", "DNA transduction map", "HIT ID", "R file format", "Sequence feature ID (SwissRegulon)", "Alignment format", "Fingerprint", "ABS ID", "Protein secondary structure", "HIX ID", "FASTQ-sanger", "Sequence database search", "Ontology format", "Sequence assembly", "Gene ID (ECK)", "Codon number", "CATH representative domain sequences (COMBS)", "Molecule interaction report", "E-value", "Tool version information", "Sequence cluster ID (SYSTERS)", "Nucleotide comparison matrix (integers)", "Toxin structure", "Gene features report (intron)", "Data retrieval", "Gene family report", "Nucleic acid sequence analysis", "InterPro secondary accession", "PubChem entry format", "Validation of peptide-spectrum matches", "GDE", "Protein structure assignment", "Nucleic acid property", "GCG MSF", "Genetic map", "Resource type", "refseqp", "Nucleic acid features report (signal or transit peptide)", "Protein extinction coefficient calculation", "Transcription factor binding site prediction", "Ensembl ID ('Homo sapiens')", "MMDB ID", "URI", "Methylation calling", "Nucleic acid features (immunoglobulin gene structure)", "File format name", "Pathway or network accession", "Structural variation discovery", "Tool log", "Biomarkers", "Protein atom surface calculation", "Allergy, clinical immunology and immunotherapeutics.", "ChIP-on-chip", "Results sort order", "Term ID list", "Cytogenetic map", "TAIR gene report format", "Genetic information processing pathways", "Protein families", "Protein globularity", "Reaction ID (KEGG)", "EMBOSS domainatrix log file", "Quantification", "unpure", "Secondary structure alignment metadata", "Database version information", "Pathway or network processing", "Sequence motif recognition", "Gene cluster", "Conversion", "Sequence conversion", "Sequence distance matrix generation", "Mass spectrometry experiment", "Carbohydrate conformational map", "Protein key folding sites", "completely unambiguous pure", "Ensembl ID ('Echinops telfairi')", "Sequence cluster ID (UniGene)", "EMBOSS simple format", "PharmGKB ID", "Amino acid name", "UniProtKB format", "Dotplot plotting", "Ecological data", "Imaging", "Codon usage table", "Urology and nephrology", "Sequence assembly", "Phylogenetic tree editing", "Report", "Genetic mapping and linkage", "DNA base pair twist angle data", "Gene ID (Genolist)", "Data identity and mapping", "Phylogeny", "Chemical structure sketch", "Protein comparison", "Gene classification", "Gel ID", "Localized reassembly", "Sequence feature detection", "Database field name", "Locus ID (PseudoCAP)", "Protein signature type", "Protein hydrogen exchange rate calculation", "Sequence record lite format", "Protein-nucleic acid interactions", "PRINTS code", "Database hits (sequence) format", "Protein modelling", "TIGRFam ID", "Promoter ID", "Sequence cluster format (protein)", "Pfam domain name", "Sequence signature report", "pbm", "Sequence accession (hybrid)", "WIG", "Job status", "ENCODE peak format", "Sequencing-based expression profile data analysis", "rcc", "Protein residue", "quicktandem", "Protein secondary database search", "RNA-seq read count analysis", "GIF", "Family name", "est2genome format", "tiff", "ProDom accession number", "Sequence visualisation", "Comparison matrix type", "myGrid", "Ontology relation type", "ASTD ID (intron)", "Structure analysis", "GeneCards gene report format", "Protein sequence visualisation", "Mutation annotation (basic)", "Opthalmology", "Constrained sequence alignment", "Ribosomes", "Expression signals", "Dimethyl", "STRING entry format", "Gene and transcript structure (report)", "ORF identifier", "Gene ID (SGN)", "Structural similarity search", "Gene name (TGD)", "Multiple nucleotide sequence alignment", "Alignment score or penalty", "Mathematics", "Secondary structure image", "Tropical medicine", "Data types and objects", "Ontology identifier", "CATH domain report format", "Physicochemical property data processing", "Documentation and help", "BIND accession number", "Structural profile processing", "Gene ID (PlasmoDB)", "Protein sequence analysis", "Pathway ID (KEGG)", "Carbohydrate property", "giFASTA format", "Secondary structure report", "DNA sequence", "BLAST XML results format", "Protein interaction raw data analysis", "mega variant", "Protein features (domains) format", "Sequence database search (by motif or pattern)", "PIRSF entry format", "Epigenomics", "Protein topological domains", "Article report", "Sequence distance matrix", "Assembly", "Sequence merging", "Gap extension penalty", "Taxonomic classification", "Pairwise structure alignment", "Gene ID (GeneDB Leishmania major)", "ID mapping", "DTD", "EMBL feature location", "InterPro protein view report format", "Structure alignment (pair)", "Nucleic acid features report (polymorphism)", "Panther Pathways entry format", "Protein molecular weight calculation", "Cheminformatics", "Nucleic acid features report (microsatellite)", "SNP", "EMBL-HTML", "Acronym", "TreeCon format", "Bit score", "KEGG PLANT entry format", "Gene name (CGSC)", "Nucleic acid curvature calculation", "Oligonucleotide probe sets annotation", "Multiple protein tertiary structure alignment", "ClustalW dendrogram", "MHC Class I epitopes report", "bgzip", "UniProt keywords", "Protein contact map calculation", "Literature and reference", "Operation", "EBI Application Result XML", "Fish", "OWL/XML", "Species frequency estimation", "Workflow data", "Protein fold name", "Epitope mapping", "Peptides and amino acids", "Lipid report", "Gap opening penalty", "Gene ID (GeneDB Glossina morsitans)", "Schema", "OME-TIFF", "Conserved transcription regulatory sequence identification", "Abstract", "Plot", "FASTA-like (text)", "Genome annotation", "Nucleic acid features report (VNTR)", "Microscope image visualisation", "PED/MAP", "Gene name (dictyBase)", "Document format", "Sequence assembly format", "OBO", "Sequence alignment", "Structure alignment (protein pair)", "tRNA structure", "Smith-Waterman format", "C-alpha trace", "ID list", "Virtual PCR", "Text processing", "Protein interaction prediction", "Primer removal", "xlsx", "Protein and peptide identification", "Map type", "TSV", "Genus name", "MAGE-TAB", "Ensembl protein ID", "Atomic data format", "Safety sciences", "ORF name", "Residue interaction calculation", "Sequence processing (protein)", "Sequence sites, features and motifs", "Surgery", "MaizeDB ID", "Microarrays", "FASTA", "DNA transcription", "Protein features report (active sites)", "Alignment format (pair only)", "Sequence analysis", "Protein structure analysis", "exp", "Molecular mass", "Cell migration track image", "Comparison matrix", "SBtab", "markx2", "Gap extension penalty (integer)", "Terminal gap opening penalty", "Phylogenetic continuous quantitative character format", "Phylogenetic footprinting / shadowing", "Structural (3D) profiles", "Sequence alignment (words)", "scores format", "Sequence search results", "NGS experiment", "Literature analysis", "Topic", "Data retrieval (protein annotation)", "nrrd", "Residue distance calculation", "Alignment data", "Sequence trimming", "Spectral counting", "Hopp and Woods plot", "Protein geometry report", "Genotype and phenotype annotation ID", "Peptide identification", "DrugBank ID", "NCBI taxon", "Nucleic acid property calculation", "Sequence-profile alignment", "Metabolic labeling", "Geotemporal metadata", "Chemical redundancy removal", "Omics", "Organism report", "Genes and proteins resources", "Protein interaction ID", "Splicing model analysis", "Operation (typed)", "Article data", "Drug discovery", "PDB", "Phylogenetic character contrasts", "FSSP entry format", "Sequence similarity plot", "Sequence mask parameter", "Entity collection identifier", "Molecular property", "EMBASSY domain classification", "Raw microarray data", "Sequence alignment metadata", "Pathway ID (reactome)", "Experiment report", "Sequence alignment report (site conservation)", "Genetic mapping", "Chimeric sequence detection", "Sequence feature label", "Statistical estimate score", "BioCyc ID", "Genotyping experiment", "Molecular surface calculation", "MEME Dirichlet prior", "pmc", "Sequence motif or profile", "Sequence trace image", "Sequence assembly visualisation", "Phylogenetic tree construction (from polymorphism data)", "Compound identifier", "Mutation identifier", "Sequence set (stream)", "Protein hydrophobic region calculation", "Comparison matrix identifier", "DeprecatedClass", "jackknifer", "Annotation", "Ligand identifier", "RGD gene report format", "Fickett testcode plot", "Modelling and simulation", "Residue contact calculation (residue-nucleic acid)", "Raw sequence", "Date", "Analysis", "Restriction sites", "CATH homologous superfamily", "Nucleic acid sequences", "ProQ report format", "Pathway ID (PATIKA)", "Sequence tag profile", "Laboratory animal science", "Protein report (enzyme) format", "RGD ID", "Protein family report format", "MIRIAM URI", "Protein pKa value", "Protein secondary structure prediction (coils)", "RDF format", "Sequence clusters and classification", "Protein fold recognition report", "Relax-NG schema", "HGVbase ID", "Peptide immunogenicity data", ".nib", "Read depth analysis", "pure protein", "Text mining", "Structural data processing", "DIP ID", "Molecule type", "mf", "TAIR accession (protein)", "Sequence-to-profile alignment", "Protein aliphatic index", "Whole gene prediction", "EMDB entry format", "protein", "Biophysics", "Disease ID (PharmGKB)", "FASTA search results format", "Protein name", "CAS number", "Repeat sequence detection", "Phylogeny reconstruction", "Neutron diffraction", "Atomic property", "Type", "Sequence processing", "Phylogeny visualisation", "Genome visualisation", "Enzyme ID (CAZy)", "Structure processing (RNA)", "Transcription regulatory sequence analysis", "User ID", "Sequence position specification", "Correlation", "UniGene taxon", "OWL Functional Syntax", "Protein structural motifs and surfaces", "Classification", "GenBank accession", "Sequence alignment generation (multiple profile)", "Experiment annotation ID", "Metabolomics", "jackknifernon", "SMART domain assignment report format", "Cell cycle", "Lipoproteins", "SCOP superfamily", "Genotype/phenotype report", "Protein surface calculation (accessible)", "Comparison matrix (nucleotide)", "MIRIAM data type name", "Protein isoelectric point", "Primer or probe design", "Molecule name", "Functional profiling", "Protein subcellular localization", "INOH entry format", "rna", "Yeast", "MINT ID", "Nucleic acid design", "ConsensusPathDB identifier", "SED-ML", "Protein non-canonical interactions", "Disease ID", "Polypeptide chain ID", "Map drawing", "Protein domain classification node", "Sequence annotation", "pair", "GPCR classification", "Electron microscopy model format", "Peptide hydrophobic moment", "Locus ID (CMR)", "Bioengineering", "Restriction enzyme name", "Molecular docking", "Gene expression data", "Structural biology", "Nucleic acid classification", "Molecule report", "EMBLXML", "medline", "Ancestral reconstruction", "Transcription factor identifier", "Chemical name (brand)", "MSF", "Comparison matrix (floats)", "Annotation retrieval (sequence)", "Protein sequence properties plot", "Database comparison", "Entity feature identifier", "Ordered locus name", "Pathway ID (INOH)", "Functional mapping", "Musculoskeletal medicine", "Amino acid name (three letter)", "Pathway or network identifier", "Sequence annotation", "Sequence attribute", "Medicines research and development", "CDD ID", "Pathway ID (DQCS)", "2D PAGE report", "Nucleic acid enthalpy", "Logical operator", "Structure-based sequence alignment", "Drug identifier", "Rate of association", "KEGG GENES gene report format", "Peptide database search", "18O labeling", "Identifier (typed)", "Proteomics experiment", "Structure prediction", "Ensembl ID (Homo sapiens)", "Sequence range", "Protein geometry calculation", "CABRI catalogue name", "Sequence cluster ID (CluSTr)", "Protein crystallizability", "Chromosome report", "BAM", "Gene expression profile analysis", "Protein features report (repeats)", "Protein folding analysis", "Protein residue surface calculation (accessible molecular)", "Data retrieval (pathway or network)", "Ensembl ID ('Dasypus novemcinctus')", "HET group name", "Format validation", "Hepatic and biliary medicine", "Nucleic acid folding report", "G protein-coupled receptors (GPCR)", "Run number", "Ontology visualisation", "Disease identifier", "Sequence trace", "Protein sequence record (lite)", "Peptide identifier", "Nucleotide base annotation", "DNA substitution model", "Sequence signature map", "SS", "Protein folds and structural domains", "Protein family report", "Gene transcript report", "BioXSD", "Ensembl ID ('Gasterosteus aculeatus')", "NCI-Nature pathway entry format", "Image annotation", "Transcriptome assembly (mapping)", "pdbseqresnuc", "Transmembrane protein prediction", "DAS sequence feature annotation", "Pathway ID (Unipathway)", "Ensembl ID ('Ciona intestinalis')", "Tool name (FASTA)", "InterPro entry abstract format", "xbm", "Gap penalty", "Vienna local RNA secondary structure format", "EMBOSS graph", "RFAM accession", "Sequencing metadata name", "Sequence motifs", "Phylogenetic tree bootstrapping", "P-value", "Server metadata", "Protein feature prediction (from sequence)", "Domain-domain interaction (indirect)", "Vienna RNA concentration data", "SILAC", "PIRSF ID", "Gene expression profile processing", "Alignment format (text)", "Pure mathematics", "Inhibitor annotation", "Map format", "Protein-drug interactions", "COGEME unisequence ID", "TCID", "Map", "Phasing", "Ontology metadata", "JASPAR profile ID", "Sequence database search (by isoelectric point)", "Sequence database cross-references", "Rotamer likelihood prediction", "Nucleic acid density plotting", "Pharmacogenomics", "Sequence tagged sites", "Protein residue surface calculation (accessible)", "Gamma diversity data", "Protein secondary structure alignment generation", "Protein signature", "Sequence classification", "Pathology", "UTRSite ID", "unambiguous pure rna sequence", "Named entity recognition", "Disease pathway or network report", "De-novo assembly", "Protein family ID (GeneFarm)", "TreeFam accession number", "Sequence retrieval", "Transcriptome assembly", "URI format", "Protein atom", "Hierarchy", "Transcriptional features (report)", "ipi", "DOI", "ENZYME enzyme report format", "Lipid identifier", "siRNA binding specificity prediction", "DNA base pair stacking energies data", "Ensembl gene report format", "Gene report", "Split read mapping", "Data acquisition", "sif", "Arabidopsis", "Metabolic network modelling", "Sequence clustering", "Microarray cluster textual view generation", "Biological model format", "ORF ID", "SMILES string", "FMA", "QSAR descriptor (constitutional)", "Nucleic acid design", "Salt bridge calculation", "Gene ID (Virginia microbial)", "Protein structure alignment", "Phylogenetic species tree construction", "Data handling", "Data index analysis", "CGD gene report format", "Public health and epidemiology", "Generation", "Database name (SwissRegulon)", "Prosite protein pattern", "AAindex ID", "Document similarity calculation", "Radiation Hybrid Mapping", "GO (cellular component)", "HGVbase entry format", "Mutation annotation (functional)", "RNA-Seq analysis", "Gene expression and microarray", "Hidden Markov model", "Workflow metadata", "NMR spectrum", "HumanCyc entry format", "Gene expression profile comparison", "completely unambiguous pure protein", "Annotation retrieval", "InterPro accession", "Protein interaction network prediction", "Protein binding site prediction (from structure)", "Regular expression", "Integrated gene prediction", "UMLS", "UniSTS accession", "Statistical inference", "Nucleic acid stitch profile plotting", "MeSH", "Gene name (NCBI)", "Scaffold gap completion", "Protein ID (ConoServer)", "Staden format", "Protein ID (TopDB)", "Gene resources", "Genetic code name", "Linkage analysis", "SISYPHUS ID", "Protein hydropathy calculation", "Stock number (TAIR)", "Strain name", "Prosite nucleotide pattern", "Family name (virus)", "REBASE withrefm enzyme report format", "DAS format", "Base frequencies table", "Local sequence alignment", "PDB ID", "CPDB entry format", "Protein family name", "Oligonucleotide probe annotation", "Sequence motif discovery", "Mutation ID", "GO (molecular function)", "DASGFF", "Ontology", "Mutation annotation (prognostic)", "Molecular dynamics", "Data resource definition", "Restriction site creation", "Information content matrix", "TreeCon-seq", "InterPro entry name", "Directory metadata", "GlycoMap ID", "Ontology concept reference", "Protein features report (disordered structure)", "Reaction kinetics ID (SABIO-RK)", "PDB atom record format", "PSL", "Phylogenetic tree report (tree distances) format", "Sequence alignment file processing", "GO", "Phylogenetic tree comparison", "Plant ontology term", "Personalized medicine", "Software engineering", "Protein secondary structure comparison", "Domain-domain interactions", "CiteXplore-core", "Phylip distance matrix", "Clone ID (RefSeq)", "EMBOSS repeat", "HMMER-aln", "Protein features", "Labeled quantification", "COGEME EST ID", "mFLJ/mKIAA number", "Ion counting", "Gene ID (GeneDB Trypanosoma brucei)", "Phylogenetic consensus tree", "Phylogenetic tree generation (from continuous quantitative characters)", "Epigenetics", "Aligned sequence order", "Protein peeling", "Preclinical and clinical studies", "hssp", "Animal study", "PCR primer design (for gene transcription profiling)", "Sequence masking", "Virus ID", "ChEBI entry format", "Functional, regulatory and non-coding RNA", "Protein binding site signature", "Northern blot report", "3D profile-to-3D profile alignment (pairwise)", "Radiation hybrid map", "pcx", "Sequence retrieval (by keyword)", "Prediction and recognition (protein)", "Amino acid identifier format", "Structure alignment (nucleic acid)", "Gene name (CGD)", "Protein structure prediction", "DNA structural variation", "Sequence cluster (nucleic acid)", "Gene order", "Protein features (PEST sites)", "Sequence generation", "KEGG LIGAND entry format", "Data retrieval (restriction enzyme annotation)", "MGED", "Alignment analysis", "Gene features report (exon)", "Cell line name (assonant)", "Endocrinology and metabolism", "Ramachandran plot calculation", "Protein solubility prediction", "Protein ID (LGICdb)", "acedb", "Protein design", "Protein chemical modifications", "OWL format", "Mass spectra calibration", "geneseq", "HTML", "Gene synonym", "QSAR descriptor (molecular)", "Sequence alignment analysis", "Coding RNA", "AGP", "Gene expression profile pathway mapping", "imzML", "TreeFam format", "Tau angle calculation", "Sequence coordinate conversion", "Protein residue surface calculation (vacuum molecular)", "phylip sequence format", "Gene annotation format", "Gene name (SGD)", "Clone library", "Hybrid sequence alignment (pair)", "GCT/Res format", "msf alignment format", "Sequence signature data", "Protein structure (all atoms)", "DRCAT resource", "Nucleic acid sites, features and motifs", "Protein structure prediction", "Resource metadata", "Gene name (Bacillus subtilis)", "Matrix", "CRAM", "Pathway or network visualisation", "NMR", "Pathway or network name", "unambiguous pure", "Organism identifier", "Raw image", "Protein secondary structure report", "Protein domain (C-alpha atoms)", "Nucleic acid temperature profile", "Blind peptide database search", "Exonic splicing enhancer prediction", "RNA secondary structure prediction", "Proteolytic digest", "Sequence generation (nucleic acid)", "REBASE enzyme number", "Genbank common name", "Atomic occupancy", "Phylogenetic tree generation (method centric)", "InChIKey", "Primer3 mispriming library file", "Identifier", "Profile-profile alignment", "Gene regulatory network report", "qualillumina", "Protein features (post-translation modifications)", "MaxQuant APL peaklist format", "Transcription", "Nucleic acid feature detection", "Gene expression correlation analysis", "Biological pathway or network format", "dbid", "Format", "plain text format (unformatted)", "Vienna RNA structural data", "Nucleic acid structure data", "Sequence alignment metadata (quality report)", "Sequence similarity score", "srs format", "Protein features report (key folding sites)", "UMLS vocabulary", "Match reward score", "Protein structure comparison", "Mathematical model", "Microbiology", "BIOM format", "CT", "Immunogenicity prediction", "Phylogenetic character weights", "Accession", "Mismatch penalty score", "Nucleic acid structure analysis", "Phylogenetic tree analysis (gene family prediction)", "Phylogenetic tree format (XML)", "Structure formatting", "FASTA-like", "Enzymes", "Proline mutation value calculation", "Protein features report (nucleic acid binding sites)", "QSAR descriptor", "Functional enrichment", "Repeat sequence analysis", "UNII", "Gap separation penalty (float)", "Enrichment", "KEGG REACTION enzyme report format", "pkl", "Amino acid index ID", "Protein-drug interaction report", "Pairwise structure alignment generation (global)", "Transcription factor binding sites", "Structure retrieval (by code)", "Protein isoelectric point calculation", "Sequence alignment parameter", "Free cysteine detection", "Terminal gap penalty", "Sequence mask character", "Comparison", "CHP", "KEGG GLYCAN entry format", "Image analysis", "Pathway or network analysis", "Profile-to-profile alignment", "Promoters", "Network simulation", "Protein secondary structure visualisation", "DPVweb ID", "Protein classification", "Quality affairs", "completely unambiguous pure rna sequence", "Structural (3D) profile ID", "Mammals", "NCBI genetic code ID", "RNA structure", "Sequence record format (text)", "Map feature", "Protein hydropathy cluster calculation", "Sequence motif data", "Protein features (epitopes)", "Sequence signature model", "Nucleic acid sequence composition (report)", "Synthetic chemistry", "Gynaecology and obstetrics", "Sequence feature annotation format", "Genome alignment", "FASTQ-illumina", "Graph format", "X!Tandem XML", "Sequence quality report format (text)", "Tool topic", "Database entry", "Color", "Gap opening penalty (float)", "Citation", "bed6", "Small molecule data processing", "Protein optical density", "BDML", "phyloXML", "Sequence composition, complexity and repeats", "Sequence feature source", "Protein secondary structure analysis", "Locus ID (CGD)", "Protein secondary structure prediction", "Isotope-coded protein label", "Flies", "Metagenomics", "Protein interaction networks", "Molecular medicine", "Toxicology", "Nucleic acid thermodynamic data", "Enzyme name", "Molecular data", "Metadata retrieval", "Codon usage table name", "Matrix format", "Gene expression report ID", "Protein analysis", "pure", "Sequence composition calculation (nucleic acid)", "File name extension", "Protein property", "Protein binding site prediction", "Epitope mapping (MHC Class I)", "Comparison matrix (integers)", "Helical net", "SAGE data processing", "Sequence profile name", "Methylation level analysis (global)", "Molecule identifier", "Protein atom surface calculation (accessible)", "Protein tertiary structure prediction", "Pathway ID (BioCyc)", "Nucleic acid features report (mutation)", "MIRIAM identifier", "Drug name", "Sequence set ID", "Protein cleavage sites", "Nucleic acid thermodynamic property calculation", "netCDF", "DNA", "Peptide molecular weights", "Sequence data", "Phylogenetic tree report (tree stratigraphic)", "Text mining report format", "Molecular dynamics simulation", "fitch program", "Protein sites, features and motifs", "SCOP node", "Data resource identifier", "Ensembl ID", "debug", "KEGG object identifier", "NCBI format", "RFAM name", "Ramachandran plot validation", "ICD identifier", "Trim ends", "Cell biology", "Protein solubility", "Phylogenetic tree generation (quartet methods)", "Phylogenetic tree generation (from gene frequencies)", "Variant prioritization", "bigBed", "affymetrix", "Bioinformatics", "HMMER2", "Catalogue ID", "RFLP", "MAP", "Tool name (EMBASSY package)", "Oscar3", "Gene ID (KOME)", "Residue symmetry contact calculation", "Chemical registry number (Gmelin)", "ABI", "Protein cleavage sites and proteolysis", "Phylogenetic tree ID", "GO concept name (cellular component)", "Km", "Phylogenetic report", "Viruses", "Text", "Chromosome name", "Disease (specific)", "Nucleic acid features (stem loop)", "Structural profile", "DICOM format", "Amino acid annotation", "Sequence editing", "unambiguous sequence", "Sequence feature detection", "Nucleic acid data", "Z-value", "RNA family identifier", "Protein-ligand interactions", "GFF", "Unicellular eukaryotes", "Sequence alignment formatting", "Window step size", "Statistics and probability", "BSML", "Residue contact calculation (residue-residue)", "Vienna RNA structure constraints", "Nucleic acid identifier", "Physiology", "TAIR accession (gene)", "SCOP concise classification string (sccs)", "Transcription factor name", "Sequence generation (protein)", "SCOP sunid", "Machine learning", "Swiss-Prot to PDB mapping", "N-Triples", "Cysteine bridge detection", "Chromosomes", "Sequence classification", "Cellular process pathways report", "Database entry metadata", "Nc statistic", "BlotBase blot ID", "Phenotype name", "Alignment", "Non-coding RNA", "Superfamily entry format", "Nucleic acid structure prediction", "Tomography", "Sequence range format", "3D-1D scoring matrix generation", "ModelDB ID", "PubMed ID", "Tool name (signature)", "BRENDA organism ID", "Tool", "Gene regulatory network processing", "Codon usage fraction difference", "Respiratory medicine", "Vienna RNA parameters", "snpeffdb", "Nucleic acid stitch profile", "SCOP fold", "PubChem ID", "Hidden Markov model", "Directory name", "Molecular interactions, pathways and networks", "Ensembl ID ('Cavia porcellus')", "Tool name (BLAST)", "Protein function prediction (from sequence)", "Reaction ID (Rhea)", "Amino acid word frequencies table", "Sequence width", "Nucleic acid features report (binding)", "microRNA detection", "Atom ID", "Nucleic acid folding analysis", "CATH domain ID", "KEGG COMPOUND entry format", "Nomenclature", "NeuronDB ID", "DDBJ accession", "Phylip cliques format", "Nucleotide code", "Gene homology (report)", "CATH version information", "Biodiversity data format", "BioPAX", "Genome index", "Experiment annotation format", "Training material", "RDF/XML", "Protein variants", "Splice sites", "Trim to reference", "MRM/SRM", "Blot ID", "Genotype and phenotype data", "Operon prediction", "PDBML", "MGI accession", "Tool name", "Cancer type", "Sequence features (compositionally-biased regions)", "Vmax", "Error", "Protein sequence cleavage", "Codon usage bias plot", "Protein signal peptide detection (bacteria)", "Nucleic acid melting temperature", "Data quality management", "DNA base trimer roll angles data", "Codon usage data", "SCOP class", "iRefIndex ID", "myGrid concept ID", "Emission matrix", "InterPro hits format", "Ensembl ID ('Canis familiaris')", "PCR primer design (for methylation PCRs)", "Discrete entity identifier", "Gene ID (ZFIN)", "Protein features report (domains)", "Plasmid map drawing", "COG sequence cluster format", "Phylogenetic invariants", "Taverna workflow format", "GenBank-like format (text)", "BioNumbers ID", "Protein interaction data processing", "RNA family report", "Systems biology", "Evidence", "newick", "FASTA-HTML", "Nucleic acid melting profile", "Genus name (virus)", "Protein structure assignment (from NMR data)", "Mutation type", "MHT", "Cell culture collection", "Signal transduction pathway report", "Physiology parameter", "Gene expression profile clustering", "Species name", "Compound ID (3DMET)", "Raw SCOP domain classification format", "Dinucleotide property", "Sequence composition table", "Database", "Comparison matrix name", "IntAct accession number", "DBD ID", "Protein dipole moment calculation", "Plasmid identifier", "PSI-PAR", "Linkage disequilibrium (report)", "Data retrieval (protein family annotation)", "Mutation annotation (prevalence)", "RNA-Seq", "MIRIAM datatype", "EMBL feature", "Primers", "Peak detection", "Protein interaction format", "DiProDB ID", "Protein interaction experiment", "Nucleic acid structure", "k-mer counting", "GenomeReviews ID", "Nucleic acid structure report", "pgSnp", "Optimisation and refinement", "KEGG PATHWAY entry format", "Sequence motif identifier", "ACE", "Expressed gene list", "Anaesthesiology", "Domain-nucleic acid interaction report", "Lane identifier", "Sequence record full", "LGICdb identifier", "Annotation", "Nucleic acid sequence record (lite)", "Protein ID (EcID)", "Metal-bound cysteine detection", "Nucleic acid features", "Protein interaction report", "HGNC concept ID", "3D profile-to-3D profile alignment", "Sequence profiles and HMMs", "Protein feature detection", "Molecular weights standard fingerprint", "2D PAGE data", "Sequence complexity report", "LipidBank ID", "Isolation report", "Protein residue 3D cluster", "Transmembrane protein visualisation", "Genetic code identifier", "MPSS data processing", "Chemical biology", "Functional clustering", "RNAVirusDB ID", "Literature search", "Sequence data processing", "Ensembl transcript ID", "rRNA", "NCBI taxonomy vocabulary", "GFF2-seq", "Pharmacology", "Protein interaction network processing", "Sequence alignment (nucleic acid)", "Affymetrix probe sets information library file", "Psychiatry", "nexus-seq", "Protein circular dichroism (CD) spectroscopic data", "Metabolic pathway report", "HMMER hidden Markov model ID", "Microsatellites", "Sequence motif comparison", "Sequence profile", "Ensembl ID ('Tupaia belangeri')", "TreeBASE study accession number", "ZFIN gene report format", "Fungi annotation (anamorph)", "ChEBI", "EcoCyc gene report format", "Strain identifier", "Panther Families and HMMs entry format", "Sequence alignment image", "Peptide mass fingerprint", "Helical wheel drawing", "Spectrum", "Cellular process pathways", "Protein sequence comparison", "QSAR descriptor (topological)", "Plants", "Stock number", "Sequence motif rendering", "Gene ID (WormBase)", "Sequence database search (by sequence using word-based methods)", "completely unambiguous pure dna", "Score end gaps control", "HMM emission and transition counts format", "Sequence alignment conversion", "Gene ID (Gramene)", "Gene name (HGNC)", "Nucleic acid features report (RFLP)", "Sequence alignment (protein pair)", "Nucleic acid restriction digest", "Transmembrane protein prediction", "Sequence databases", "Sequence set (polymorphic)", "Gene name (MaizeGDB)", "SNP detection", "Gene name (MGD)", "Cell line name (no punctuation)", "unambiguous pure nucleotide", "PCR primer design (based on gene structure)", "smarts", "Protein domain", "Protein domain (all atoms)", "GFF3", "Veterinary medicine", "Structural genomics", "DNA binding sites", "Article metadata", "HMMER synthetic sequences set", "TraML", "markx1", "Comparative genomics", "vectorstrip cloning vector definition file", "Microarray data rendering", "iHOP organism ID", "Protein structure", "Sequence motif analysis", "QTL map", "Sequence motif format", "TAIR accession (At gene)", "PDB insertion code", "Primer3 internal oligo mishybridizing library", "PNG", "Protein features report (cleavage sites)", "Sample collections", "restrict format", "Phylogenetic tree generation (AI methods)", "Sequence trace format", "Molecular genetics", "Medical informatics", "Nucleic acid restriction", "Sequence profile data", "ENCODE broad peak format ", "Residue contact calculation (residue-ligand)", "Cytoscape input file format", "afg", "Phylogenetic property values format", "Variant calling", "Phylogenetic gene frequencies data", "Methylation analysis", "Structure processing", "Gene name (EcoGene primary)", "Nucleic acid sequence alignment", "Kinase name", "Cell line name (truncated)", "Virulence prediction", "Phylogenetic property values", "Number of iterations", "Protein chain", "Protein chain (all atoms)", "Biological pathway map", "Full torsion angle calculation", "Psiblast checkpoint file", "csfasta", "Protein chain (C-alpha atoms)", "Structural distance matrix", "PCR primers", "CellML", "Nucleic acid folding", "Biological system modelling", "Heat map generation", "Protein sequences", "Protein function analysis", "Protein residue surface calculation (vacuum accessible)", "Local structure alignment", "Pathway or network prediction", "Pcons report format", "TAIR accession", "Protein report format", "BioPax term", "X-ray diffraction", "Codon usage fraction calculation", "Protein features (sites)", "BCF", "unambiguous pure protein", "Computational biology", "Cell type ontology ID", "Reference identification", "Polymorphism report format", "NCBI genome accession", "Genetics", "Protein threading", "InChI", "Gene ID (GeneDB Schizosaccharomyces pombe)", "nucleotide", "Sequence alignment editing", "UniParc accession", "EMBOSS vectorstrip log file", "NCBI gene report format", "PTM identification", "Protein ID (PeroxiBase)", "Matrix/scaffold attachment site prediction", "Multiple protein sequence alignment", "Experimental data (proteomics)", "Sequence alignment (nucleic acid pair)", "Sequence ambiguity calculation", "Data retrieval (codon usage table)", "Sequence set", "Domainatrix 3D-1D scoring matrix format", "Protein report (membrane protein)", "Protein active site signature", "Protein post-translational modifications", "Informatics", "Diffraction data integration", "GI number (protein)", "Protein alignment", "selex", "strider format", "PDB atom name", "Protein atom surface calculation (accessible molecular)", "Gene name (ASPGD)", "Nucleic acid structure", "Gene ID (VBASE2)", "Complementary medicine", "Codon usage bias", "Sequence cluster (protein)", "Biological imaging", "Genetic codes and codon usage", "Information retrieval", "Nucleic acid sites and features", "DNA structure", "Deletion map", "raw", "Map set data", "Residue non-canonical interaction detection", "NCBI version", "tRNA gene prediction", "Sequence signature matches", "Gene set testing", "HGNC vocabulary", "Laboratory techniques", "Gene features (exonic splicing enhancer)", "Protein targeting and localization", "HMMER3", "Sequence set (protein)", "GeneDB gene report format", "Genome identifier", "Toxin annotation", "Protein super-secondary structure", "Structure classification", "Gene features (SECIS element)", "Phylogenetic tree analysis (shape)", "Sequence alignment", "Indexing", "Biomaterials", "Database search results", "Peptide property", "Root-mean-square deviation", "Gap opening penalty (integer)", "EPS", "GTrack", "Sequence metadata", "hdf5", "Protein-protein interactions", "Locus ID (MGG)", "Global structure alignment", "Secondary structure", "Word composition", "Paediatrics", "Experimental design and studies", "Data mining", "Pathway or network report", "Cytoband format", "Microarray Box-Whisker plot plotting", "Protein family identifier", "dbGaP format", "Gene expression profiling", "Sequence editing (protein)", "Database search (by sequence)", "Compound libraries and screening", "Enumerated file name", "Helical wheel", "High-throughput sequencing", "Scaffolding", "Sequence record full format", "Sequence cluster format", "Residue interaction prediction", "Sequence identifier", "Monosaccharide accession", "Data submission, annotation and curation", "GenBank-like format", "Pileup", "Protein-protein interactions", "Sequence coordinates", "qual454", "Global sequence alignment", "GeminiSQLite", "Clustering", "Protein accession", "Enzyme kinetics calculation", "Entity identifier", "Codon usage bias plotting", "OMIM ID", "Imputation", "Compound ID (KEGG)", "Phylogenetic tree annotation", "Protein charge plot", "RNA secondary structure alignment", "Gene name (DragonDB)", "Proteomics", "PHD", "markx0 variant", "Strain accession", "Reactome entry format", "Microarray wave graph plotting", "Protein sequence", "Sequence alignment", "PO", "Chemical registry number (Beilstein)", "Transcription regulatory element prediction (trans)", "Phenomics", "Mouse clinic", "Nucleic acids", "Phylogenetic tree report (cliques) format", "Phylogenetic distance matrix identifier", "Phylogenetic tree report (tree evaluation)", "Sequence name", "Online course", "Theoretical biology", "Protein architecture analysis", "Protein sequence-structure scoring matrix", "Ensembl gene tree ID", "Structure similarity score", "SBS data processing", "Protein data", "Computer science", "Structure alignment ID", "HMMER format", "Protein structure comparison", "Sequence profile processing", "Sequence-structure alignment", "Phylogenetic tree distances", "Gene features (coding region) format", "Chemical name (synonymous)", "Compound accession", "Protein sequence hydropathy plot", "Nucleic acid probability profile plotting", "Pathogens", "Gene ID (MIPS Medicago)", "SBML", "Transcription factor accession", "Protein structure image", "PolyA signal or sites", "Nucleic acid sequence (raw)", "JPG", "Chemical registry number (CAS)", "Nucleosome formation potential prediction", "PRIDE XML", "Protein feature detection", "Promoters", "Differential binding analysis", "mega-seq", "Structure ID", "Protein architecture comparison", "Sequence alignment ID", "Tertiary structure format", "Locus ID", "Phylogenetic tree report (tree shape)", "Structure alignment report", "Bibliographic reference format", "Sequence-3D profile alignment", "DictyBase gene report format", "Phylogenetic data", "Text mining report", "Phylogenetic tree generation (from molecular sequences)", "Classification report", "Protein family accession", "PCR primer design (for genotyping polymorphisms)", "Nucleic acid sequence record (full) ", "MetaCyc entry format", "Prediction and recognition", "HIVDB identifier", "UniProt accession (extended)", "Sequence feature type", "Protein property calculation (from sequence)", "RSF", "mspcrunch", "Data retrieval (sequence profile)", "GPCR analysis", "docx", "Pathway or network", "Copy number estimation", "Disease name", "Microarray probe design", "Sequence checksum", "MAGE-ML", "Hierarchy identifier", "Hopp and Woods plotting", "Small molecule data", "Protein fold recognition", "Staden experiment format", "Pharmacovigilence", "Northern blot image", "Sample ID", "PseudoCAP gene report format", "Protein features report (topological domains)", "Brite hierarchy ID", "Map annotation", "DNA vaccine design", "Protein features report (chemical modifications)", "Locus ID (SGD)", "Phylogenetic tree generation (data centric)", "Parasitology", "PMML", "Transition matrix", "MGED concept ID", "Dirichlet distribution", "Invertebrates", "Protein report", "Nucleic acid features report (expression signal)", "Feature table", "Ensembl ID ('Bos taurus')", "Phylogenetic character cliques", "Nucleic acid property processing", "Locus ID (MMP)", "Molecular interaction data processing", "DNA mutation", "Sequence database search (by sequence for primer sequences)", "Protein structure report (quality evaluation) format", "Structure alignment", "MRI", "Tool identifier", "UniProt keyword", "Protein domain classification", "ColiCard report format", "Genome version information", "Strain data format", "Cysteine torsion angle calculation", "Enzyme report", "Gramene identifier", "Sequence feature table format", "Gene ID (NCBI)", "Topology diagram drawing", "Animals", "Base position variability plotting", "Base position variability plot", "Operating system name", "Gene ID (CGD)", "Drug formulation and delivery", "Variant classification", "Protein interaction networks", "Trauma medicine", "Methylated DNA immunoprecipitation", "pcd", "dssp", "WormBase wormpep ID", "1 or more", "Gene ID (DictyBase)", "Nucleosome formation or exclusion sequence prediction", "Antimicrobial resistance prediction", "Format identifier", "Genetic marker identification", "Medline UI", "affymetrix-exp", "Tool name (EMBOSS)", "Nucleic acid features report (replication and recombination)", "Structure analysis", "Translational medicine", "HGMD ID", "DSSP secondary structure assignment", "ArachnoServer ID", "Quantitative trait locus", "ISA-TAB", "Compound ID (ChEMBL)", "Sequence composition calculation", "SCOP domain identifier", "Protein interaction raw data", "Northern blot experiment", "markx10", "UMLS concept ID", "Pairwise sequence alignment generation (local)", "GO concept ID (molecular function)", "Flow cell identifier", "Enzyme ID (MEROPS)", "phylipnon", "Nucleic acid features (d-loop)", "Pathway ID (SMPDB)", "MeSH concept ID", "Text mining", "Sequence similarity search", "GlycomeDB ID", "FMA concept ID", "Protein secondary structure prediction", "Phylogenetic character data format", "Translation frame specification", "Protein contact map", "Sequence profile type", "Biosafety report", "Identifier (hybrid)", "Phylogenetic tree generation", "Isolation source", "mhd", "GVF", "Environmental information processing pathway report", "Kingdom name", "Sequence property (nucleic acid)", "ISBN", "Structure alignment processing", "Embryo report", "TIGR identifier", "HGNC", "Affymetrix probe sets library file", "Gastroenterology", "Sequence feature comparison", "Pearson format", "igstrict", "OMIM entry format", "2D PAGE image", "GeneSNP ID", "Protein NMR data", "Prosite accession number", "Sequence assembly validation", "Relationship inference", "Epitope mapping (MHC Class II)", "BLAST sequence alignment type", "Sequence motif", "Fungi annotation", "PCR primer design (for large scale sequencing)", "Phylogenetic tree type", "Ensembl ID ('Takifugu rubripes')", "MGD gene report format", "Biobank", "EMAP", "Operon", "Structure image", "Mass spectrometry spectra", "debug-seq", "Plotting", "ProDom entry format", "QSAR descriptor (geometrical)", "GWAS study", "Profile-to-profile alignment (pairwise)", "Maximum occurence analysis", "PubChem bioassay ID", "Database management", "Amino acid frequencies table", "EGA accession", "Person identifier", "Indel detection", "dbProbe ID", "meme-motif", "Protein family ID (PANTHER)", "Physical mapping", "Sequence property", "Phylogenetic tree format (text)", "FASTQ-like format", "phylip property values", "Electron microscopy model", "Environmental information processing pathways", "Sequence alignment type", "Chromosome name (BioCyc)", "Protein interaction network rendering", "CleanEx dataset code", "Biology", "ELM ID", "dta", "Ensembl ID ('Spermophilus tridecemlineatus')", "Genome report", "JASPAR format", "BRENDA enzyme report format", "Gene name (AceView)", "Single particle analysis", "Nucleic acid analysis", "codata", "Biosafety classification", "GFF2", "Tag mapping", "GenBank format", "Nucleic acid entropy", "Surface rendering", "QSAR descriptor (electronic)", "Multiple sample visualisation", "cPath ID", "Transmembrane protein database search", "Base word frequencies table", "Cultivation parameter", "Sequence annotation track format", "Map identifier", "Structure processing (nucleic acid)", "HET group dictionary entry format", "Ontology data", "Protein interactions", "Oligonucleotide ID", "Gene ID (miRBase)", "Individual genetic data format", "Structure comparison", "Protein ID (DisProt)", "Cell type identifier", "Heat map", "Gene ID (GeneDB)", "Protein ID (CuticleDB)", "Nucleic acid folding prediction (alignment-based)", "Transmembrane protein analysis", "Sequence profile format", "Person name", "Phylogenetic tree generation (maximum likelihood and Bayesian methods)", "Pathway or network annotation", "HMMER Dirichlet prior", "Residue bump detection", "Protein secondary structure image", "Sequin format", "Exons", "Mapping", "Embryology", "Protein domain signature", "Nucleic acid sequence reverse and complement", "Orpha number", "Sequence features (repeats)", "EPD ID", "Protein-ligand docking", "UNITE accession", "Microarray spots image", "Ensembl ID ('Mus musculus')", "Demonstration", "NCBI Genome Project ID", "dna", "Sequence cluster ID (COG)", "Eukaryotes", "InterPro detailed match image", "Ensembl ID ('Otolemur garnettii')", "Formatting", "Whole microarray graph plotting", "Transcript ID", "UTR accession", "EMBOSS whichdb log file", "Email address", "Raw CATH domain classification", "Heterogen annotation", "Geriatric medicine", "nexusnon alignment format", "Reference map name (SWISS-2DPAGE)", "IntAct entry format", "Sequence alignment (hybrid)", "Quantitative genetics", "Gene ID (JCVI)", "Protein architecture recognition", "Protein interaction network analysis", "Applied mathematics", "EMDB ID", "Protein structural motif recognition", "Reaction ID (SABIO-RK)", "Sequence alignment (multiple)", "primersearch primer pairs sequence record", "daf", "Carbohydrates", "Infectious tropical disease", "Whole genome sequencing", "RNA structure", "Sequence database search (by molecular weight)", "Protein hydropathy calculation (from sequence)", "Small molecule structure", "Isotropic B factor", "Structure alignment (RNA)", "Protein sequence repeats", "pdbatomnuc", "GI number", "Nucleic acid folding energy calculation", "Statistical calculation", "Haplotype map", "Biotechnology", "Sequence profile generation", "Gene expression data analysis", "PDB database entry format", "Microarray principal component plotting", "Protein-motif interaction", "dbSNP polymorphism report format", "CCAP strain number", "2D PAGE experiment", "Database name (Osteogenesis)", "Structure retrieval (water)", "Phylogenetic tree visualisation", "Sequence set (bootstrapped)", "Read mapping", "Ab initio structure prediction", "Coding region", "DNA-Seq", "Sequence redundancy removal", "Drug report", "2bit", "RNA structure covariance model generation", "Protein super-secondary structure prediction", "Restriction digest", "Protein binding site prediction (from sequence)", "Protein conserved site signature", "Sequence alignment (pair)", "Locus ID (ASPGD)", "Binary format", "PCR primer design (for conserved primers)", "BacMap gene card format", "Data handling", "Neurite measurement", "Sequence cluster ID (UniRef100)", "UniProt format", "Pathway or network data", "BioCyc enzyme report format", "Haematology", "Protein hydropathy data", "Rice", "Phylogenetics", "xgmml", "Protein disordered structure", "QSAR descriptor name", "psd", "Dot-bracket format", "Sequence word comparison", "Clone ID", "InterPro entry format", "SFF", "Translation initiation site prediction", "Sequence comparison", "Splice site prediction", "Label-free quantification", "Sequence feature table format (XML)", "Data governance", "Gene prediction", "Protein classification", "Specific protein resources", "REBASE restriction sites", "Protein features report (signal peptides)", "Medicinal chemistry", "Toxin accession", "Nucleic acid structure analysis", "Protein hydropathy", "EST assembly", "Sequence length specification", "Drop off score", "Reaction data", "Molecular similarity score", "pure dna", "FASTQ-like format (text)", "Protein folding site prediction", "Structural clustering", "doc2loc document information", "Cell line name (exact)", "Spectral analysis", "Pairwise sequence alignment generation (global)", "DNA sequence (raw)", "Sequence cluster", "Index format", "Metabolic disease", "Oncology", "Filtering", "Calculation", "Protein secondary structure format", "Distance matrix", "EMBOSS sequence pattern", "Data retrieval (database cross-reference)", "Atomic coordinate", "STRING entry format (XML)", "Sequence alignment analysis (site correlation)", "Protein-nucleic acid interactions", "GCG format variant", "Sequence mask type", "Sequence file editing", "Genotype and phenotype", "Toxin name", "DNA sense specification", "Pathway ID (BioSystems)", "JSON-LD", "Protein structure assignment (from X-ray crystallographic data)", "Turtle", "EMAGE ID", "Gramene gene report format", "Protein-protein interaction prediction", "Sequencing-based expression profile data processing", "Nucleic acid features report (CpG island and isochore)", "Domainatrix signature", "SMART entry format", "SRF", "Protein membrane regions", "MHC peptide immunogenicity prediction", "Gene ID (MIPS Maize)", "Protein secondary structure prediction (turns)", "iHOP symbol", "MRI image", "Sequence motif (nucleic acid)", "Sequence image", "Pairwise protein tertiary structure alignment (all atoms)", "Multiple nucleic acid tertiary structure alignment", "Protein sequence (raw)", "HMMER profile alignment (HMM versus sequences)", "Sequence motif (protein)", "Sequence retrieval (by code)", "Phylip discrete states format", "Protein image", "Drug development", "Physics", "Medline Display Format", "bedstrict", "Gene name (GeneFarm)", "Sequence signature identifier", "xls", "Rendering parameter", "Secondary structure alignment generation", "Alignment", "Nucleic acid sequence feature detection", "Operon drawing", "Genetic variation", "Protein structure analysis", "Structural (3D) profile alignment", "Codon usage table comparison", "Map data", "Nucleic acid features report (matrix/scaffold attachment sites)", "pepXML", "MeSH vocabulary", "aaindex", "Protein function prediction", "Phosphorylation sites", "Phylogenetic sub/super tree detection", "REBASE proto enzyme report format", "PDB file sequence retrieval", "swiss feature", "Nucleic acid repeats", "URN", "ASN.1 sequence format", "Transcription regulatory element prediction (DNA-cis)", "mzQuantML", "3D-1D scoring matrix", "Model organisms", "Sequence contamination filtering", "VDB", "Protein sequence record (full)", "Sequence-profile alignment (HMM) format", "Single particle alignment and classification", "Data retrieval (sequence alignment)", "aMAZE entry format", "T-Coffee format", "Regnerative medicine", "Structure alignment (nucleic acid pair)", "Cell type accession", "Genotype and phenotype annotation format", "Pairwise sequence alignment", "Scent annotation", "Molecular interaction analysis", "Data resource definition accession", "Matrix identifier", "FASTQ-solexa", "Species tree", "OBO format", "insdxml", "PolyA signal detection", "Sequence analysis", "Sequence feature table format (text)", "NCBI locus tag", "Textual format", "MAT", "Sequence database name", "Protein interaction network comparison", "MEME motifs directive file", "Deisotoping", "Peptide annotation", "Study topic", "Phylogenetic tree reconstruction", "Secondary structure alignment (protein)", "ChemSpider entry format", "Secondary structure alignment metadata (RNA)", "Organism name", "Phylogenetic discrete data", "YAML", "BLC", "Database entry identifier", "Nucleic acid structure prediction", "Mice or rats", "DNA packaging", "Name", "A2M", "Sequence features metadata", "Sequence assembly format (text)", "Amino acid index (van der Waals radii)", "Data retrieval (tool metadata)", "Genetic code prediction", "Job metadata", "Random sequence generation", "Structure", "Protein function comparison", "Locus ID (Tropgene)", "Dermatology", "Drug structure relationship map", "Molecular surface comparison", "Sequence alignment analysis (conservation)", "WormBase identifier", "3D profile generation", "Exactly 1", "Amino acid index (molecular weight)", "qualsolexa", "EMBOSS Uniform Sequence Address", "COSMIC ID", "Nucleic acid repeats (report)", "Computational chemistry", "Locus ID (MaizeGDB)", "Sequence alignment processing", "Biochemistry", "Mascot .dat file", "Sequence", "chrominfo", "Accessible surface calculation", "mzTab", "SBGN-ML", "TreeBASE format", "Protein domain ID", "Sequence motif matches (nucleic acid)", "Genome indexing (suffix arrays)", "Gene transcripts", "Plant Ontology concept ID", "GO concept ID (biological process)", "meganon", "ABCD format", "Compound ID (ChemIDplus)", "GO concept ID", "Small molecule report format", "Protein features (mutation)", "Article comparison", "Restriction enzyme report", "Two-dimensional gel electrophoresis", "Codon usage bias calculation", "Book ID", "Pairwise structure alignment generation (local)", "Sequence composition report", "Mass spectrum visualisation", "Carbohydrate identifier", "CATH node ID", "Keyword", "Peptide ID", "CATH PDB report format", "Peptide molecular weight hits", "Pathway ID (PharmGKB)", "Geographic location", "Gene identifier (Entrez)", "Query and retrieval", "2 or more", "Molecules", "bigWig", "Simulation experiment", "Sequence distance matrix format", "Cardinality", "Protein sequence analysis", "Data index", "Sequence cluster format (nucleic acid)", "pgm", "Toxin identifier", "Sequence type", "Developmental biology", "Sequence submission", "BindingDB Monomer ID", "Locus ID (EntrezGene)", "GO concept name (molecular function)", "ENCODE narrow peak format", "Exome assembly", "arff", "Sequence record lite", "Sanger inverted repeats", "RNA secondary structure visualisation", "Phylip continuous quantitative characters", "Protein architecture report", "bmp", "Sequence motif matches (protein)", "Structure alignment", "Differential gene expression analysis", "Experimental measurement", "BioPax concept ID", "GO concept name", "Probabilistic data generation", "Exome sequencing", "Sample annotation", "Protein folding report", "Hanes Woolf plot", "CpG island and isochore detection", "pure nucleotide", "Protein globularity prediction", "Sequence-profile alignment (HMM)", "Sequence ambiguity report", "T3DB ID", "Ontology concept data", "Article", "Sequence complexity calculation", "markx0", "Lipid accession", "Mobile genetic elements", "Alignment", "Biomedical science", "Protein properties", "Structure database search", "Function analysis", "FlyBase ID", "Protein titration curve plotting", "Nucleic acid melting curve plotting", "Regulatory RNA", "Chemistry", "Electron microscopy volume map", "Human genetics", "Sequence information report", "Sequence profile ID", "Microarray tree-map rendering", "NCBI taxonomy ID", "Protein hydrogen bonds", "Pathogenicity report", "Virtual PCR", "Protein secondary structure", "hennig86", "FASTQ", "Ontology comparison", "Structure database search (by sequence)", "Phylip format variant", "Disease genes and proteins", "Microarray hybridisation data", "MatrixDB interaction ID", "ChEBI concept ID", "Nucleic acid features (difference and change)", "RESID ID", "Molecular docking", "Taverna workflow ID", "Gene identifier (NCBI UniGene)", "BioModel ID", "Protein cleavage site prediction", "Morphology parameter", "Nucleic acid features (codon)", "Prokaryotes and Archaea", "Protein structure report", "Reaction ID (MACie)", "OBO file format name", "Protein flexibility or motion report", "Protein-nucleic acid interactions report", "Phylogenetic tree processing", "Organelles", "Gene ID (SGD)", "Codon usage data processing", "Target ID (TTD)", "PubChem CID", "Proteins", "Regulatory affairs", "RefSeq accession", "Sequence feature qualifier", "Protein solvent accessibility report", "OBO-XML", "InterPro match table format", "Peptide identification", "Ensembl ID ('Xenopus tropicalis')", "CDD PSSM-ID", "Map data processing", "Immunology", "Sequence motif matches", "Nexus format", "RNA secondary structure image", "Transcription regulatory element prediction (RNA-cis)", "EMBOSS listfile", "Fate map", "Gene cluster format", "Database category name", "RefSeq accession (protein)", "Gene expression profile", "Oligonucleotide probe data", "tabix", "WIFF format", "Anatomy", "Phylip tree raw", "SAGE experimental data", "Ensembl ID ('Monodelphis domestica')", "Polymorphism detection", "Sequence cluster ID", "Cytoband position", "PSI MI XML (MIF)", "BioGRID interaction ID", "InterPro compact match image", "Ensembl ID ('Oryzias latipes')", "Password", "EMBL-like format", "GenBank-HTML", "Data index processing", "Spot ID", "PRINTS entry format", "Kinetic model", "Chemical name (ChEBI)", "Protein distance matrix", "Deposition", "Dirichlet distribution format", "IMGT/HLA ID", "Signaling pathways", "Codon usage analysis", "Protein fragment weight comparison", "Reaction ID", "Sequence property (protein)", "CleanEx entry name", "Gene symbol", "DNA polymorphism", "Genetic information processing pathway report", "Gene name", "Obsolete concept (EDAM)", "PDB model number", "Peak calling", "Critical care medicine", "Protein-ligand complex", "Fungi", "Protein active sites", "Protein X-ray crystallographic data", "Virus annotation (taxonomy)", "Lipid structure", "Carbohydrate report", "Image format", "Carbohydrate structure", "CATH representative domain sequences (ATOM)", "Small molecules", "Gene features report (operon)", "Protein interaction analysis", "Protein property calculation", "Drug structure", "ConsensusPathDB entity name", "Variant filtering", "Laboratory information management", "meganon sequence format", "Database name (CMD)", "Article format", "ClustalW format", "Genotype experiment ID", "Hit sort order", "Gene regulation", "Threshold", "Sequence clustering", "BED", "Ontology concept", "Protein site signature", "Signal or transit peptide", "Transcription factors and regulatory sites", "Data visualisation", "BAI", "Protein post-translational modification signature", "SMART domain name", "LAV", "Sequence formatting", "Atomic z coordinate", "PRIDE experiment accession number", "iHOP text mining abstract format", "Nucleic acid sequence analysis", "Residue cluster calculation", "cel", "Protein report (function)", "Primer3 primer"]

termsDict = {"event_bioschemas": "Event (Bioschemas)",
             "event_schema": "Event (schema.org)",
             "organization_bioschemas": "Organization (Bioschemas)",
             "organization_schema": "Organization (schema.org)",
             "person_bioschemas": "Person (Bioschemas)",
             "person_schema": "Person (schema.org)",
             "training_bioschemas": "Training (Bioschemas)",
             "training_schema": "Training (schema.org)"}

class TagsReference:
    """Extract all the Bioschemas and schema.org properties from the Bioschemas website."""

    def __init__(self):
        self.refTags = {}

        self.getTags()
        self.writeCSV()
        self.writeJSON()

    def getTags(self):
        """Get the latest properties from the Bioschemas website."""

        print "Updating Bioschemas properties..."

        kinds = ["event", "organization", "person", "training"]

        for kind in kinds:
            base_url = "http://bioschemas.org/groups"
            end_url = "%ss/full_%s.html" % (kind, kind)
            url = "%s/%s" % (base_url, end_url)

            response = requests.get(url)
            html = response.content
            soup = BeautifulSoup(html, "lxml")

            bioschemas_table = soup.find("table", attrs={"class": "bioschemas"})
            schema_table = soup.find("table", attrs={"class": "schema"})
            thing_table = soup.find("table", attrs={"class": "thing"})

            td_list = []
            type_list = []
            for row in bioschemas_table.findAll("tr"):
                row_list = []
                for cell in row.findAll("td"):
                    row_list.append(cell.text.encode("utf-8"))
                td_list.append(row_list)
            for i in xrange(1, len(td_list)):
                prop_dict = {}
                prop_dict[td_list[i].pop(0)] = td_list[i]
                type_list.append(prop_dict)
            self.refTags["%s_bioschemas" % kind] = type_list

            td_list = []
            type_list = []
            for row in schema_table.findAll("tr"):
                row_list = []
                for cell in row.findAll("td"):
                    row_list.append(cell.text.encode("utf-8"))
                td_list.append(row_list)
            for i in xrange(1, len(td_list)):
                prop_dict = {}
                prop_dict[td_list[i].pop(0)] = td_list[i]
                type_list.append(prop_dict)
            self.refTags["%s_schema" % kind] = type_list

            td_list = []
            type_list = []
            for row in thing_table.findAll("tr"):
                row_list = []
                for cell in row.findAll("td"):
                    row_list.append(cell.text.encode("utf-8"))
                td_list.append(row_list)
            for i in xrange(1, len(td_list)):
                prop_dict = {}
                prop_dict[td_list[i].pop(0)] = td_list[i]
                type_list.append(prop_dict)
            self.refTags["%s_schema" % kind] += type_list

    def writeCSV(self):
        """Save all the Bioschemas properties in a CSV file."""

        print "Saving Bioschemas properties in bioschemasTags.csv..."

        # For each type, save their properties in a specific CSV file

        sorted_props = sorted(self.refTags.keys())
        for el in sorted_props:
            try:
                f = open("bioschemas/%s.csv" % el, "wb")
            except IOError:
                os.system("mkdir bioschemas")
                f = open("bioschemas/%s.csv" % el, "wb")
            writer = csv.writer(f)
            writer.writerow(["Property", "Type", "Description", "Cardinality", "Guideline", "Vocabulary"])
            for prop in self.refTags[el]:
                row_to_write = []
                for key in prop:
                    row_to_write.append(key.encode("utf-8"))
                    for value in prop[key]:
                        row_to_write.append(value.encode("utf-8"))
                writer.writerow(row_to_write)
            f.close()

        # Save the number and details of the properties for each Bioschemas type, along with the number of minimum, recommended and optional properties

        f = open("bioschemas/bioschemasTags.csv", "wb")
        writer = csv.writer(f)
        writer.writerow(["Type", "Properties", "Details", "Minimum", "Recommended", "Optional"])
        for el in sorted_props:
            countMin = 0
            countRec = 0
            countOpt = 0
            details_list = []
            for prop in self.refTags[el]:
                details_list.append(prop.keys()[0])
                if prop.values()[0][3] == "Minimum":
                    countMin += 1
                elif prop.values()[0][3] == "Recommended":
                    countRec += 1
                elif prop.values()[0][3] == "Optional":
                    countOpt += 1
            writer.writerow([el, len(details_list), details_list, countMin, countRec, countOpt])
        f.close()

    def writeJSON(self):
        """Save all the Bioschemas properties in a JSON file."""

        endfile = open("bioschemas/bioschemasTags.json", "wb")
        json.dump(self.refTags, endfile, indent=4, sort_keys=True)
        endfile.close()

class WebsiteTags:
    """Extract all the properties found in the scraped website."""

    def __init__(self, url, nome, sibs=1):
        self.url = url
        self.nome = nome
        self.sibs = sibs
        self.websiteTypes = set()
        self.foundTags = {
            "event_bioschemas": [],
            "event_schema": [],
            "organization_bioschemas": [],
            "organization_schema": [],
            "person_bioschemas": [],
            "person_schema": [],
            "training_bioschemas": [],
            "training_schema": []
        }
        self.tagsGuide = {
            "event_bioschemas": [],
            "event_schema": [],
            "organization_bioschemas": [],
            "organization_schema": [],
            "person_bioschemas": [],
            "person_schema": [],
            "training_bioschemas": [],
            "training_schema": []
        }
        self.tagsData = {
            "event_bioschemas": [],
            "event_schema": [],
            "organization_bioschemas": [],
            "organization_schema": [],
            "person_bioschemas": [],
            "person_schema": [],
            "training_bioschemas": [],
            "training_schema": []
        }
        self.allProps = []          # All the props found in the website
        self.validProps = []        # Valid props found in the website

        self.checkFiles()
        self.checkWebsite()

    def checkFiles(self):
        """Check if the requested files exist, otherwise they will be created."""

        try:
            sourcefile = open("scrapedWebsites.csv", "rb")
        except IOError:
            sourcefile = open("scrapedWebsites.csv", "wb")
            siteWriter = csv.writer(sourcefile)
            siteWriter.writerow(["Website", "Name"])
        finally:
            sourcefile.close()

    def checkWebsite(self):
        """Check whether the website has already been scraped before or not."""

        try:
            scrapedSites = open("scrapedWebsites.csv", "rb")
            siteReader = csv.reader(scrapedSites)

            for row in siteReader:
                if row[0] == self.url:
                    self.updateWebsiteTags()
                    break
                else:
                    continue
            else:
                self.saveWebsite()
            scrapedSites.close()
        except IOError:
            self.saveWebsite()

    def saveWebsite(self):
        """Save the current website URL and name in the scrapedWebsites.csv file."""

        scrapedSites = open("scrapedWebsites.csv", "ab")
        siteWriter = csv.writer(scrapedSites)
        siteWriter.writerow(["%s" % self.url, "%s" % self.nome])
        scrapedSites.close()

        self.newWebsiteTags()

    def newWebsiteTags(self):
        """Create a new directory named after the domain of the scraped website."""

        if self.nome == self.url:
            websiteComps = self.url.split("/")
            websiteName = "_".join(websiteComps[2:])
        else:
            websiteName = self.nome

        os.system("mkdir %s" % websiteName)

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        if soup.findAll(itemprop=True):
            self.scrapeMicrodata(websiteName)
            self.validateWebsiteMicrodata(websiteName)
            if self.sibs == 0:
                self.validateSiblings()
        else:
            self.scrapeRDFa(websiteName)
            self.validateWebsiteRDFa(websiteName)
            if self.sibs == 0:
                self.validateSiblings()

    def updateWebsiteTags(self):
        """Update the data from previously scraped websites."""

        if self.nome == self.url:
            websiteComps = self.url.split("/")
            websiteName = "_".join(websiteComps[2:])
        else:
            websiteName = self.nome

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        if soup.findAll(itemprop=True):
            self.scrapeMicrodata(websiteName)
            self.validateWebsiteMicrodata(websiteName)
            if self.sibs == 0:
                self.validateSiblings()
        else:
            self.scrapeRDFa(websiteName)
            self.validateWebsiteRDFa(websiteName)
            if self.sibs == 0:
                self.validateSiblings()

    def scrapeMicrodata(self, sitename):
        """Get the microdata from the website and save them into a JSON file."""

        reference = TagsReference()
        print "Connecting to %s..." % self.url

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        print "Getting the Types found in %s..." % self.url

        for el in soup.findAll(itemtype=True):
            webtype = el.get("itemtype").split("/")[3]
            self.websiteTypes.add(webtype)

        endfile = open("%s/typesFound.txt" % sitename, "wb")
        for el in self.websiteTypes:
            endfile.write(el + "\n")
        endfile.close()

        print "Scraping microdata from %s..." % self.url

        for typeProp in reference.refTags:
            prop_dict = {}
            guide_dict = {}
            for property in reference.refTags[typeProp]:
                for el in soup.findAll(itemprop=True):

                    prop = el.get("itemprop")

                    if prop == "image":
                        cont = el.get("src")
                    elif prop == "logo":
                        cont = el.get("content")
                    elif prop == "sameAs" or prop == "url":
                        cont = el.get("href")
                    else:
                        cont = el.text.strip().replace("\n", "").replace("\r", "").replace("\t", "").encode("utf-8")

                    for key in property:
                        if prop == key:
                            try:
                                prop_dict[prop] += [cont]
                            except KeyError:
                                prop_dict[prop] = [cont]

                for key in property:
                    if key in prop_dict:
                        cont = prop_dict[key]
                        if cont == None:
                            pass
                        else:
                            cont = str(cont[0])

                            # Type values
                            if property[key][0] == "Integer":
                                if cont.isdigit():
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            elif property[key][0] == "Number":
                                if cont.replace(".", "", 1).isdigit() or cont.replace(",", "", 1).isdigit():
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            elif property[key][0] == "Boolean":
                                if cont == "True" or cont == "False" or cont == "true" or cont == "false" or cont == "T" or cont == "F" or cont == "t" or cont == "f":
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            elif property[key][0] == "URL":
                                if "/" in cont or "." in cont:
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            else:
                                val0 = "OK"

                            # Cardinality
                            if property[key][2] == "One":
                                if type(cont) == str:
                                    val1 = "OK"
                                else:
                                    val1 = "Wrong"
                            else:
                                if type(cont) != str:
                                    val1 = "OK"
                                else:
                                    val1 = "Wrong"

                            # Vocabulary
                            if property[key][4] == "Yes":
                                if cont.capitalize() in edamList:
                                    val2 = "OK"
                                else:
                                    val2 = "Wrong"
                            else:
                                val2 = "OK"


                            # Guide dict
                            guide_dict[key] = [property[key][0], property[key][1], property[key][2], property[key][3], property[key][4], val0, val1, val2]

                    else:
                        guide_dict[key] = [property[key][0], property[key][1], property[key][2], property[key][3], property[key][4], "Not found", "Not found", "Not found"]


                self.foundTags[typeProp] = [prop_dict]
                self.tagsGuide[typeProp] = [guide_dict]

        endfile = open("%s/tagsFound.json" % sitename, "wb")
        json.dump(self.foundTags, endfile, indent=4, sort_keys=True)
        endfile.close()
        endfile2 = open("%s/tagsGuide.json" % sitename, "wb")
        json.dump(self.tagsGuide, endfile2, indent=4, sort_keys=True)
        endfile2.close()

    def validateWebsiteMicrodata(self, sitename):
        """Validate the entries found in the website for Bioschemas markup using Microdata."""

        found_types = open("%s/typesFound.txt" % sitename, "rb")    # All the types found in the website
        f = open("%s/tagsFound.json" % sitename, "rb")
        r = open("bioschemas/bioschemasTags.json", "rb")
        found = json.load(f)
        ref = json.load(r)

        report = open("%s/report.csv" % sitename, "wb")
        w_rep = csv.writer(report)
        w_rep.writerow(["Report"])

        typesToAdd = []
        propsToAdd = []
        details = []

        for line in found_types:

            miss = open("%s/missingMin.csv" % sitename, "ab")  # Required props missing in the website
            w_miss = csv.writer(miss)

            prov = open("%s/providedRec.csv" % sitename, "ab")  # Not required props provided in the website
            w_prov = csv.writer(prov)

            line = line.strip()
            if line == "Event" or line == "Organization" or line == "Person" or line == "Training":
                subj_type = line.lower()
                w_rep.writerow(["Website %s belonging to the %s type." % (self.url, line)])
                typeProps = set()
                this_detail = []
                foundMin = 0
                foundRec = 0
                foundOpt = 0
                # Found Bioschemas type

                for ref_key in ref:
                    foundDict = {}
                    if ref_key.startswith(subj_type):
                        for i in ref[ref_key]:
                            for ref_prop in i:
                                expType = i[ref_prop][0]
                                cardinality = i[ref_prop][2]
                                guideline = i[ref_prop][3]
                                vocabulary = i[ref_prop][4]

                                # Validate Minimum props

                                for k in found[ref_key]:
                                    if guideline == "Minimum" and not ref_prop in k:
                                        w_miss.writerow([ref_prop, ref_key])
                                    elif guideline == "Minimum" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundMin += 1
                                    elif guideline == "Recommended" and ref_prop in k:
                                        w_prov.writerow([ref_prop, ref_key])
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundRec += 1
                                    elif guideline == "Optional" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundOpt += 1

                                    if ref_prop in k:
                                        self.allProps.append(ref_prop)
                                        self.tagsData[ref_key].append(ref_prop)

                typesToAdd.append(line)
                propsToAdd.append(list(typeProps))
                this_detail.append([foundMin, foundRec, foundOpt])
                details.append(this_detail)

                miss.close()
                prov.close()

                w_rep.writerow(["Found %d valid properties out of %d properties provided." % (len(self.validProps), len(self.allProps))])

                miss = open("%s/missingMin.csv" % sitename, "rb")
                miss_r = csv.reader(miss)

                prov = open("%s/providedRec.csv" % sitename, "rb")
                prov_r = csv.reader(prov)

                for line in prov_r:
                    w_rep.writerow(["Missing %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                for line in miss_r:
                    w_rep.writerow(["Provided %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                miss.close()
                prov.close()

        f.close()
        r.close()
        found_types.close()
        report.close()

        UpdateRegistry().updateRegistryFile(sitename, typesToAdd, propsToAdd, details)
        UpdateRegistry().createChartFile(sitename, typesToAdd)

    def scrapeRDFa(self, sitename):
        """Get the RDFa data from the website and save them into a JSON file."""

        reference = TagsReference()
        print "Connecting to %s..." % self.url

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        print "Getting the Types found in %s..." % self.url

        for el in soup.findAll(typeof=True):
            webtype = el.get("typeof")
            self.websiteTypes.add(webtype)

        endfile = open("%s/typesFound.txt" % sitename, "wb")
        for el in self.websiteTypes:
            endfile.write(el + "\n")
        endfile.close()

        print "Scraping microdata from %s..." % self.url

        for typeProp in reference.refTags:
            prop_dict = {}
            guide_dict = {}
            for property in reference.refTags[typeProp]:
                for el in soup.findAll(property=True):

                    prop = el.get("property")

                    if len(prop.split(":")) > 1:
                        prop = prop.split(":")[-1]
                    else:
                        prop = prop.strip()

                    if prop == "image":
                        cont = el.get("src")
                    elif prop == "logo":
                        cont = el.get("content")
                    elif prop == "sameAs" or prop == "url":
                        cont = el.get("href")
                    else:
                        cont = el.text.strip().replace("\n", " ").replace("\r", " ").replace("\t", " ").encode("utf-8")

                    for key in property:
                        if prop == key:
                            try:
                                prop_dict[prop] += [cont]
                            except KeyError:
                                prop_dict[prop] = [cont]

                for key in property:
                    if key in prop_dict:
                        cont = prop_dict[key]
                        if cont == None:
                            pass
                        else:
                            cont = str(cont[0])

                            # Type values
                            if property[key][0] == "Integer":
                                if cont.isdigit():
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            elif property[key][0] == "Number":
                                if cont.replace(".", "", 1).isdigit():
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            elif property[key][0] == "Boolean":
                                if cont == "True" or cont == "False" or cont == "true" or cont == "false" or cont == "T" or cont == "F" or cont == "t" or cont == "f":
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            elif property[key][0] == "URL":
                                if "/" in cont or "." in cont:
                                    val0 = "OK"
                                else:
                                    val0 = "Wrong"
                            else:
                                val0 = "OK"

                            # Cardinality
                            if property[key][2] == "One":
                                if type(cont) == str:
                                    val1 = "OK"
                                else:
                                    val1 = "Wrong"
                            else:
                                if type(cont) != str:
                                    val1 = "OK"
                                else:
                                    val1 = "Wrong"

                            # Vocabulary
                            if property[key][4] == "Yes":
                                if cont.capitalize() in edamList:
                                    val2 = "OK"
                                else:
                                    val2 = "Wrong"
                            else:
                                val2 = "OK"

                            # Guide dict
                            guide_dict[key] = [property[key][0], property[key][1], property[key][2], property[key][3], property[key][4], val0, val1, val2]

                    else:
                        guide_dict[key] = [property[key][0], property[key][1], property[key][2], property[key][3], property[key][4], "Not found", "Not found", "Not found"]


                self.foundTags[typeProp] = [prop_dict]
                self.tagsGuide[typeProp] = [guide_dict]

        endfile = open("%s/tagsFound.json" % sitename, "wb")
        json.dump(self.foundTags, endfile, indent=4, sort_keys=True)
        endfile.close()
        endfile2 = open("%s/tagsGuide.json" % sitename, "wb")
        json.dump(self.tagsGuide, endfile2, indent=4, sort_keys=True)
        endfile2.close()

    def validateWebsiteRDFa(self, sitename):
        """Validate the entries found in the website for Bioschemas markup using RDFa."""

        found_types = open("%s/typesFound.txt" % sitename, "rb")    # All the types found in the website
        f = open("%s/tagsFound.json" % sitename, "rb")
        r = open("bioschemas/bioschemasTags.json", "rb")
        found = json.load(f)
        ref = json.load(r)

        report = open("%s/report.csv" % sitename, "wb")
        w_rep = csv.writer(report)
        w_rep.writerow(["Report"])

        typesToAdd = []
        propsToAdd = []
        details = []

        for line in found_types:

            miss = open("%s/missingMin.csv" % sitename, "ab")   # Required props missing in the website
            w_miss = csv.writer(miss)

            prov = open("%s/providedRec.csv" % sitename, "ab")  # Not required props provided in the website
            w_prov = csv.writer(prov)

            line = line.strip()
            if line == "Event" or line == "Organization" or line == "Person" or line == "Training":
                subj_type = line.lower()
                w_rep.writerow(["Website %s belonging to the %s type." % (self.url, line)])
                typeProps = set()
                this_detail = []
                foundMin = 0
                foundRec = 0
                foundOpt = 0
                # Found Bioschemas type

                for ref_key in ref:
                    foundDict = {}
                    if ref_key.startswith(subj_type):
                        for i in ref[ref_key]:
                            for ref_prop in i:
                                expType = i[ref_prop][0]
                                cardinality = i[ref_prop][2]
                                guideline = i[ref_prop][3]
                                vocabulary = i[ref_prop][4]

                                # Validate Minimum props

                                for k in found[ref_key]:
                                    if guideline == "Minimum" and not ref_prop in k:
                                        w_miss.writerow([ref_prop, ref_key])
                                    elif guideline == "Minimum" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundMin += 1
                                    elif guideline == "Recommended" and ref_prop in k:
                                        w_prov.writerow([ref_prop, ref_key])
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundRec += 1
                                    elif guideline == "Optional" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundOpt += 1

                                    if ref_prop in k:
                                        self.allProps.append(ref_prop)
                                        self.tagsData[ref_key].append(ref_prop)

                typesToAdd.append(line)
                propsToAdd.append(list(typeProps))
                this_detail.append([foundMin, foundRec, foundOpt])
                details.append(this_detail)

                miss.close()
                prov.close()

                w_rep.writerow(["Found %d valid properties out of %d properties provided." % (len(self.validProps), len(self.allProps))])

                miss = open("%s/missingMin.csv" % sitename, "rb")
                miss_r = csv.reader(miss)

                prov = open("%s/providedRec.csv" % sitename, "rb")
                prov_r = csv.reader(prov)

                for line in prov_r:
                    w_rep.writerow(["Missing %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                for line in miss_r:
                    w_rep.writerow(["Provided %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                miss.close()
                prov.close()

        f.close()
        r.close()
        found_types.close()
        report.close()

        UpdateRegistry().updateRegistryFile(sitename, typesToAdd, propsToAdd, details)
        UpdateRegistry().createChartFile(sitename, typesToAdd)

    def validateSiblings(self):
        """Collect any sibling pages and validate them."""

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")
        targetList = []

        for el in soup.findAll(href=True):
            if el.get("href").startswith(self.url):
                print el.get("href")
                targetList.append(el.get("href"))
        for i in len(targetList):
            WebsiteTags(targetList[i], "%s_%d" % (self.nome, i), 1)


class UpdateRegistry:
    """Update the registry file every time a new website is added to the scraped websites file."""

    def __init__(self):
        self.properties = []
        self.websites = []
        self.types = []

    def updateRegistryFile(self, website, type_bs, props, details):
        """Create and update the registry file."""

        finDict = {}
        elDict = {}
        for i in range(len(type_bs)):
            elDict[type_bs[i]] = details[i]
            elDict[type_bs[i]] += props[i]
        finDict[website] = [elDict]

        with open("registry.json", "rb") as f:
            data = json.load(f)

        data.update(finDict)

        with open("registry.json", "wb") as f:
            json.dump(data, f, indent=4)

    def createChartFile(self, website, type_bs):
        """Create the csv file needed to build the chart in the website page."""

        okColor = "#6CC4A4"
        noColor = "#E1514B"
        warnColor = "#FEC574"

        for typ in type_bs:
            graphList = []
            maxList = []
            okCount = 0
            noCount = 0
            warnCount = 0
            ord = 0

            with open("%s/tagsGuide.json" % website, "rb") as f, open("%s/tagsFound.json" % website, "rb") as g:
                data = json.load(f)
                data2 = json.load(g)

                if typ == "Event":
                    for i in data["event_bioschemas"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["event_bioschemas"][0][i][3]
                        actType = data["event_bioschemas"][0][i][5]
                        actCard = data["event_bioschemas"][0][i][6]
                        actVocab = data["event_bioschemas"][0][i][7]

                        try:
                            scoreGraph = len(data2["event_bioschemas"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                    for i in data["event_schema"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["event_schema"][0][i][3]
                        actType = data["event_schema"][0][i][5]
                        actCard = data["event_schema"][0][i][6]
                        actVocab = data["event_schema"][0][i][7]

                        try:
                            scoreGraph = len(data2["event_schema"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                elif typ == "Organization":
                    for i in data["organization_bioschemas"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["organization_bioschemas"][0][i][3]
                        actType = data["organization_bioschemas"][0][i][5]
                        actCard = data["organization_bioschemas"][0][i][6]
                        actVocab = data["organization_bioschemas"][0][i][7]

                        try:
                            scoreGraph = len(data2["organization_bioschemas"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                    for i in data["organization_schema"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["organization_schema"][0][i][3]
                        actType = data["organization_schema"][0][i][5]
                        actCard = data["organization_schema"][0][i][6]
                        actVocab = data["organization_schema"][0][i][7]

                        try:
                            scoreGraph = len(data2["organization_schema"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                elif typ == "Person":
                    for i in data["person_bioschemas"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["person_bioschemas"][0][i][3]
                        actType = data["person_bioschemas"][0][i][5]
                        actCard = data["person_bioschemas"][0][i][6]
                        actVocab = data["person_bioschemas"][0][i][7]

                        try:
                            scoreGraph = len(data2["person_bioschemas"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                    for i in data["person_schema"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["person_schema"][0][i][3]
                        actType = data["person_schema"][0][i][5]
                        actCard = data["person_schema"][0][i][6]
                        actVocab = data["person_schema"][0][i][7]

                        try:
                            scoreGraph = len(data2["person_schema"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                elif typ == "Training":
                    for i in data["training_bioschemas"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["training_bioschemas"][0][i][3]
                        actType = data["training_bioschemas"][0][i][5]
                        actCard = data["training_bioschemas"][0][i][6]
                        actVocab = data["training_bioschemas"][0][i][7]

                        try:
                            scoreGraph = len(data2["training_bioschemas"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                    for i in data["training_schema"][0]:
                        label = i
                        scoreGraph = 0
                        weightGraph = 0.0
                        colorGraph = ""

                        stdGuide = data["training_schema"][0][i][3]
                        actType = data["training_schema"][0][i][5]
                        actCard = data["training_schema"][0][i][6]
                        actVocab = data["training_schema"][0][i][7]

                        try:
                            scoreGraph = len(data2["training_schema"][0][i])
                        except KeyError:
                            scoreGraph = 0
                        if stdGuide == "Minimum":
                            weightGraph = 0.5
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        elif stdGuide == "Recommended":
                            weightGraph = 0.3
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        else:
                            weightGraph = 0.2
                            if actType == "OK" and actCard == "OK" and actVocab == "OK":
                                colorGraph = okColor
                                ord += 1
                                
                            elif actType == "Not found" and actCard == "Not found" and actVocab == "Not found":
                                colorGraph = noColor
                                ord += 1
                                
                            else:
                                colorGraph = warnColor
                                ord += 1
                                
                        graphList.append([ord, label, scoreGraph, weightGraph, colorGraph])
                        maxList.append(scoreGraph)

                meanVal = max(maxList)

                endfile = open("%s/chartData_%s.csv" % (website, typ), "wb")
                w = csv.writer(endfile)
                w.writerow(["order", "label", "score", "weight", "color", "mean"])

                maxFile = open("%s/maxData_%s.csv" % (website, typ), "wb")
                s = csv.writer(maxFile)
                s.writerow(["a", "b"])

                for el in graphList:
                    el.append(meanVal)
                    w.writerow(el)
                    # append the mean!
                    s.writerow([meanVal, "val"])
                    #print el
                endfile.close()
                maxFile.close()

                #with open("%s/chartData_%s.json" % (website, typ), "wb") as endfile:
                #    json.dump(graphList, endfile, indent=4)







WebsiteTags(sys.argv[1], sys.argv[2], sys.argv[3])
#UpdateRegistry().createFiles()