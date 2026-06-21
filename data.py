"""
Cascade - mini biological knowledge base.

This is a small, self-contained demo dataset so the app works fully offline.
Swap this out later for real data (Reactome, KEGG, STRING, DGIdb) via their
APIs once you're ready to go from "demo" to "real" — the algorithms in
algorithms.py don't need to change, only the data feeding them.
"""

# gene -> list of pathways it participates in
GENE_PATHWAYS = {
    "TP53":   ["Apoptosis", "Cell Cycle", "DNA Damage Response"],
    "MDM2":   ["Apoptosis", "Cell Cycle"],
    "CDKN1A": ["Cell Cycle", "DNA Damage Response"],
    "RB1":    ["Cell Cycle"],
    "CCND1":  ["Cell Cycle", "PI3K-AKT Signaling"],
    "CDK4":   ["Cell Cycle"],
    "CDK6":   ["Cell Cycle"],
    "E2F1":   ["Cell Cycle", "Apoptosis"],
    "BAX":    ["Apoptosis"],
    "BCL2":   ["Apoptosis"],
    "CASP3":  ["Apoptosis"],
    "CASP9":  ["Apoptosis"],
    "ATM":    ["DNA Damage Response", "Cell Cycle"],
    "ATR":    ["DNA Damage Response"],
    "BRCA1":  ["DNA Damage Response", "Apoptosis"],
    "BRCA2":  ["DNA Damage Response"],
    "PIK3CA": ["PI3K-AKT Signaling"],
    "AKT1":   ["PI3K-AKT Signaling", "Apoptosis"],
    "PTEN":   ["PI3K-AKT Signaling"],
    "MTOR":   ["PI3K-AKT Signaling"],
    "EGFR":   ["PI3K-AKT Signaling", "MAPK Signaling"],
    "KRAS":   ["MAPK Signaling", "PI3K-AKT Signaling"],
    "BRAF":   ["MAPK Signaling"],
    "MAP2K1": ["MAPK Signaling"],
    "MAPK1":  ["MAPK Signaling"],
    "MYC":    ["Cell Cycle", "MAPK Signaling"],
    "VEGFA":  ["Angiogenesis"],
    "FLT1":   ["Angiogenesis"],
    "KDR":    ["Angiogenesis"],
    "HIF1A":  ["Angiogenesis", "Metabolism"],
    "PDGFRB": ["Angiogenesis", "MAPK Signaling"],
    "IL6":    ["Inflammation"],
    "TNF":    ["Inflammation", "Apoptosis"],
    "NFKB1":  ["Inflammation", "Apoptosis"],
    "STAT3":  ["Inflammation", "MAPK Signaling"],
    "IL1B":   ["Inflammation"],
    "HK2":    ["Metabolism"],
    "PKM":    ["Metabolism"],
    "LDHA":   ["Metabolism"],
    "SLC2A1": ["Metabolism"],
}

# protein-protein interaction edges (undirected), used to build the network
PPI_EDGES = [
    ("TP53", "MDM2"), ("TP53", "CDKN1A"), ("TP53", "BAX"), ("TP53", "ATM"),
    ("TP53", "BRCA1"), ("MDM2", "RB1"), ("CDKN1A", "CDK4"), ("CDKN1A", "CDK6"),
    ("RB1", "E2F1"), ("CDK4", "CCND1"), ("CDK6", "CCND1"), ("CCND1", "MYC"),
    ("E2F1", "MYC"), ("BAX", "BCL2"), ("BAX", "CASP9"), ("CASP9", "CASP3"),
    ("BCL2", "AKT1"), ("ATM", "ATR"), ("ATM", "BRCA1"), ("ATR", "BRCA2"),
    ("BRCA1", "BRCA2"), ("PIK3CA", "AKT1"), ("AKT1", "MTOR"), ("AKT1", "PTEN"),
    ("PTEN", "PIK3CA"), ("MTOR", "HIF1A"), ("EGFR", "PIK3CA"), ("EGFR", "KRAS"),
    ("KRAS", "BRAF"), ("BRAF", "MAP2K1"), ("MAP2K1", "MAPK1"), ("MAPK1", "MYC"),
    ("KRAS", "PIK3CA"), ("VEGFA", "FLT1"), ("VEGFA", "KDR"), ("FLT1", "KDR"),
    ("KDR", "PDGFRB"), ("HIF1A", "VEGFA"), ("PDGFRB", "MAPK1"), ("IL6", "STAT3"),
    ("TNF", "NFKB1"), ("NFKB1", "IL1B"), ("STAT3", "MAPK1"), ("NFKB1", "BAX"),
    ("TNF", "CASP3"), ("HIF1A", "HK2"), ("HK2", "PKM"), ("PKM", "LDHA"),
    ("HIF1A", "SLC2A1"), ("SLC2A1", "HK2"),
]

# pathway -> drug targets known/explored for it (demo mapping)
PATHWAY_DRUG_TARGETS = {
    "Apoptosis":              ["BCL2 (Venetoclax)", "MDM2 (Idasanutlin)"],
    "Cell Cycle":             ["CDK4/6 (Palbociclib)", "CDK4/6 (Ribociclib)"],
    "DNA Damage Response":    ["PARP (Olaparib)", "ATR (Berzosertib)"],
    "PI3K-AKT Signaling":     ["PI3K (Alpelisib)", "mTOR (Everolimus)", "AKT (Capivasertib)"],
    "MAPK Signaling":         ["BRAF (Vemurafenib)", "MEK (Trametinib)"],
    "Angiogenesis":           ["VEGFA (Bevacizumab)", "VEGFR2 (Ramucirumab)"],
    "Inflammation":           ["TNF (Adalimumab)", "IL-6R (Tocilizumab)"],
    "Metabolism":             ["LDHA (GSK2837808A, investigational)"],
}

# Total number of genes considered "in the genome" for enrichment background.
# In a real tool this would be ~20,000; we use a smaller number so the demo
# dataset produces meaningful p-values.
BACKGROUND_GENE_COUNT = 500
