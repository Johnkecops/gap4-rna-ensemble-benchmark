"""Reference metadata, transcribed from PubMed records retrieved 2026-11-04 (see reference_validation.tsv).
Fields: authors (list, 'Initials Surname'), title, journal (abbrev.), year, volume, pages, doi, pmid."""
REFS = {
 "mccaskill1990": (["J.S. McCaskill"], "The equilibrium partition function and base pair binding probabilities for RNA secondary structure", "Biopolymers", 1990, "29", "1105–1119", "10.1002/bip.360290621", "1695107"),
 "ding2003": (["Y. Ding", "C.E. Lawrence"], "A statistical sampling algorithm for RNA secondary structure prediction", "Nucleic Acids Res.", 2003, "31", "7280–7301", "10.1093/nar/gkg938", "14654704"),
 "waymentsteele2022": (["H.K. Wayment-Steele", "W. Kladwang", "A.I. Strom", "J. Lee", "A. Treuille", "A. Becka", "et al."], "RNA secondary structure packages evaluated and improved by high-throughput experiments", "Nat. Methods", 2022, "19", "1234–1242", "10.1038/s41592-022-01605-0", "36192461"),
 "he2024": (["S. He", "R. Huang", "J. Townley", "R.C. Kretsch", "T.G. Karagianes", "D.B.T. Cox", "et al."], "Ribonanza: deep learning of RNA structure through dual crowdsourcing", "bioRxiv [Preprint]", 2024, "", "2024.02.24.581671", "10.1101/2024.02.24.581671", "38464325"),
 "cordero2012": (["P. Cordero", "J.B. Lucks", "R. Das"], "An RNA Mapping DataBase for curating RNA structure mapping experiments", "Bioinformatics", 2012, "28", "3006–3008", "10.1093/bioinformatics/bts554", "22976082"),
 "lorenz2011": (["R. Lorenz", "S.H. Bernhart", "C. Höner zu Siederdissen", "H. Tafer", "C. Flamm", "P.F. Stadler", "I.L. Hofacker"], "ViennaRNA Package 2.0", "Algorithms Mol. Biol.", 2011, "6", "26", "10.1186/1748-7188-6-26", "22115189"),
 "zhang2023": (["H. Zhang", "S. Li", "L. Zhang", "D.H. Mathews", "L. Huang"], "LazySampling and LinearSampling: fast stochastic sampling of RNA secondary structure with applications to SARS-CoV-2", "Nucleic Acids Res.", 2023, "51", "e7", "10.1093/nar/gkac1029", "36401871"),
 "mei2025": (["X. Mei", "H. Liu", "Y. Zhu", "E. Zhao", "L. Li", "H. Zhang"], "TVAE-RNA: ensemble-based RNA secondary structure prediction via transformer variational autoencoders", "Bioinformatics", 2025, "41", "btaf527", "10.1093/bioinformatics/btaf527", "40981507"),
 "sacco2026": (["G. Sacco", "J. Li", "R.P. Smyth", "G. Sanguinetti", "G. Bussi"], "MERGE-RNA: a physics-based model to predict RNA secondary structure ensembles with chemical probing", "Nucleic Acids Res.", 2026, "54", "gkag917", "10.1093/nar/gkag917", "42770328"),
 "mathews2004rna": (["D.H. Mathews"], "Using an RNA secondary structure partition function to determine confidence in base pairs predicted by free energy minimization", "RNA", 2004, "10", "1178–1190", "10.1261/rna.7650904", "15272118"),
 "lyngso2000": (["R.B. Lyngsø", "C.N.S. Pedersen"], "RNA pseudoknot prediction in energy-based models", "J. Comput. Biol.", 2000, "7", "409–427", "10.1089/106652700750050862", "11108471"),
 "turner2010": (["D.H. Turner", "D.H. Mathews"], "NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure", "Nucleic Acids Res.", 2010, "38", "D280–D282", "10.1093/nar/gkp892", "19880381"),
 "mathews2004pnas": (["D.H. Mathews", "M.D. Disney", "J.L. Childs", "S.J. Schroeder", "M. Zuker", "D.H. Turner"], "Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure", "Proc. Natl. Acad. Sci. U.S.A.", 2004, "101", "7287–7292", "10.1073/pnas.0401799101", "15123812"),
 "ditommaso2017": (["P. Di Tommaso", "M. Chatzou", "E.W. Floden", "P.P. Barja", "E. Palumbo", "C. Notredame"], "Nextflow enables reproducible computational workflows", "Nat. Biotechnol.", 2017, "35", "316–319", "10.1038/nbt.3820", "28398311"),
 "li2006": (["W. Li", "A. Godzik"], "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", "Bioinformatics", 2006, "22", "1658–1659", "10.1093/bioinformatics/btl158", "16731699"),
 "fu2012": (["L. Fu", "B. Niu", "Z. Zhu", "S. Wu", "W. Li"], "CD-HIT: accelerated for clustering the next-generation sequencing data", "Bioinformatics", 2012, "28", "3150–3152", "10.1093/bioinformatics/bts565", "23060610"),
 "huang2019": (["L. Huang", "H. Zhang", "D. Deng", "K. Zhao", "K. Liu", "D.A. Hendrix", "D.H. Mathews"], "LinearFold: linear-time approximate RNA folding by 5'-to-3' dynamic programming and beam search", "Bioinformatics", 2019, "35", "i295–i304", "10.1093/bioinformatics/btz375", "31510672"),
 "zhang2020": (["H. Zhang", "L. Zhang", "D.H. Mathews", "L. Huang"], "LinearPartition: linear-time approximation of RNA folding partition function and base-pairing probabilities", "Bioinformatics", 2020, "36", "i258–i267", "10.1093/bioinformatics/btaa460", "32657379"),
 "incarnato2018": (["D. Incarnato", "E. Morandi", "L.M. Simon", "S. Oliviero"], "RNA Framework: an all-in-one toolkit for the analysis of RNA structures and post-transcriptional modifications", "Nucleic Acids Res.", 2018, "46", "e97", "10.1093/nar/gky486", "29893890"),
 "zubradt2017": (["M. Zubradt", "P. Gupta", "S. Persad", "A.M. Lambowitz", "J.S. Weissman", "S. Rouskin"], "DMS-MaPseq for genome-wide or targeted RNA structure probing in vivo", "Nat. Methods", 2017, "14", "75–82", "10.1038/nmeth.4057", "27819661"),
 "siegfried2014": (["N.A. Siegfried", "S. Busan", "G.M. Rice", "J.A.E. Nelson", "K.M. Weeks"], "RNA motif discovery by SHAPE and mutational profiling (SHAPE-MaP)", "Nat. Methods", 2014, "11", "959–965", "10.1038/nmeth.3029", "25028896"),
 "deigan2009": (["K.E. Deigan", "T.W. Li", "D.H. Mathews", "K.M. Weeks"], "Accurate SHAPE-directed RNA structure determination", "Proc. Natl. Acad. Sci. U.S.A.", 2009, "106", "97–102", "10.1073/pnas.0806929106", "19109441"),
 "lowweeks2010": (["J.T. Low", "K.M. Weeks"], "SHAPE-directed RNA secondary structure prediction", "Methods", 2010, "52", "150–158", "10.1016/j.ymeth.2010.06.007", "20554050"),
 "matthews1975": (["B.W. Matthews"], "Comparison of the predicted and observed secondary structure of T4 phage lysozyme", "Biochim. Biophys. Acta", 1975, "405", "442–451", "10.1016/0005-2795(75)90109-9", "1180967"),
 "hanley1982": (["J.A. Hanley", "B.J. McNeil"], "The meaning and use of the area under a receiver operating characteristic (ROC) curve", "Radiology", 1982, "143", "29–36", "10.1148/radiology.143.1.7063747", "7063747"),
 "schober2018": (["P. Schober", "C. Boer", "L.A. Schwarte"], "Correlation coefficients: appropriate use and interpretation", "Anesth. Analg.", 2018, "126", "1763–1768", "10.1213/ANE.0000000000002864", "29481436"),
 "kraemer2006": (["H.C. Kraemer"], "Correlation coefficients in medical research: from product moment correlation to the odds ratio", "Stat. Methods Med. Res.", 2006, "15", "525–545", "10.1177/0962280206070650", "17260922"),
 "rosner2007": (["B. Rosner", "R.J. Glynn"], "Interval estimation for rank correlation coefficients based on the probit transformation with extension to measurement error correction of correlated ranked data", "Stat. Med.", 2007, "26", "633–646", "10.1002/sim.2547", "16596580"),
 "henderson2005": (["A.R. Henderson"], "The bootstrap: a technique for data-driven statistics. Using computer-intensive analyses to explore experimental data", "Clin. Chim. Acta", 2005, "359", "1–26", "10.1016/j.cccn.2005.04.002", "15936746"),
 "capra2007": (["J.A. Capra", "M. Singh"], "Predicting functionally important residues from sequence conservation", "Bioinformatics", 2007, "23", "1875–1882", "10.1093/bioinformatics/btm270", "17519246"),
 "rivas1999": (["E. Rivas", "S.R. Eddy"], "A dynamic programming algorithm for RNA structure prediction including pseudoknots", "J. Mol. Biol.", 1999, "285", "2053–2068", "10.1006/jmbi.1998.2436", "9925784"),
 "lu2009": (["Z.J. Lu", "J.W. Gloor", "D.H. Mathews"], "Improved RNA secondary structure prediction by maximizing expected pair accuracy", "RNA", 2009, "15", "1805–1813", "10.1261/rna.1643609", "19703939"),
 "wilkinson2016": (["M.D. Wilkinson", "M. Dumontier", "I.J. Aalbersberg", "G. Appleton", "M. Axton", "A. Baak", "et al."], "The FAIR Guiding Principles for scientific data management and stewardship", "Sci. Data", 2016, "3", "160018", "10.1038/sdata.2016.18", "26978244"),
}


def fmt(key):
    a, t, j, y, v, p, doi, pmid = REFS[key]
    names = ", ".join(x for x in a if x != "et al.") + (", et al." if "et al." in a else "")
    vol = f" {v}" if v else ""
    return f"{names}, {t}, {j}{vol} ({y}) {p}. https://doi.org/{doi}. PMID: {pmid}."
