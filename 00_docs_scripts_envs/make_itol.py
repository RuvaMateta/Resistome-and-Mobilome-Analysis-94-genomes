import pandas as pd

palette = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf',
           '#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#1b9e77','#d95f02']

meta = pd.read_csv("merged/metadata_final.csv")
mlst = pd.read_csv("mlst_forced.tsv", sep="\t", header=None, usecols=[0, 2], names=["file", "ST"], dtype=str)
mlst["sample"] = mlst["file"].str.replace("input/", "", regex=False).str.replace(".fasta", "", regex=False)
meta = meta.merge(mlst[["sample", "ST"]], on="sample", how="left")
meta["year"] = meta["year_final"].round().astype("Int64").astype(str)

pc = pd.read_csv("merged/plasmid_cluster_amr_integron.tsv", sep="\t")
amr = pd.read_csv("merged/amr_all.tsv", sep="\t")
acq = amr[(amr["type"] == "AMR") & (amr["subtype"] == "AMR")]
els = pd.read_csv("merged/integron_elements_all.tsv", sep="\t")
complete_int = set(els[els["type"] == "complete"]["sample"])

PLASMIDS = ["plasmid_AC082", "plasmid_AE928", "plasmid_AB042", "plasmid_AB595", "plasmid_AA277", "plasmid_AA402"]
GENES = ["blaCTX-M-15", "blaCTX-M-27", "qnrB19", "sul1", "sul2", "tet(A)", "dfrA7", "dfrA14"]

def strip(path, label, series):
    series = series.dropna()
    cats = sorted(series.unique())
    cols = {c: palette[i % len(palette)] for i, c in enumerate(cats)}
    with open(path, "w") as f:
        f.write(f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t{label}\nCOLOR\t#000000\nCOLOR_BRANCHES\t0\n")
        f.write(f"LEGEND_TITLE\t{label}\nLEGEND_SHAPES\t" + "\t".join(["1"] * len(cats)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(cols[c] for c in cats) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(cats) + "\nDATA\n")
        for s, v in series.items():
            f.write(f"{s}\t{cols[v]}\t{v}\n")

def binary(path, label, table):
    table = table.loc[:, table.any()]
    cols = [palette[i % len(palette)] for i in range(table.shape[1])]
    n = table.shape[1]
    with open(path, "w") as f:
        f.write(f"DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\t{label}\nCOLOR\t#000000\n")
        f.write("FIELD_SHAPES\t" + "\t".join(["1"] * n) + "\n")
        f.write("FIELD_LABELS\t" + "\t".join(str(c) for c in table.columns) + "\n")
        f.write("FIELD_COLORS\t" + "\t".join(cols) + "\nDATA\n")
        for s, row in table.iterrows():
            f.write(s + "\t" + "\t".join("1" if v else "-1" for v in row) + "\n")

trees = {"E_coli": ("tree/E_coli.treefile", "K12_MG1655_ref"),
         "K_pneumoniae": ("tree/K_pneumoniae_no119.treefile", "K_pneumoniae_ERR13744124"),
         "Salmonella": ("tree/Salmonella.treefile", "LT2_ref")}

for sp, (tf, newname) in trees.items():
    txt = open(tf).read()
    assert "Reference:" in txt, "no Reference tip in " + tf
    open(f"itol/{sp}.nwk", "w").write(txt.replace("Reference:", newname + ":"))
    m = meta[meta["species"] == sp].set_index("sample")
    strip(f"itol/{sp}_country.txt", "Country", m["country"])
    strip(f"itol/{sp}_source.txt", "Source", m["source_category"])
    strip(f"itol/{sp}_year.txt", "Year", m["year"])
    cnt = m["ST"].value_counts()
    strip(f"itol/{sp}_ST.txt", "ST (singletons grey = other)",
          m["ST"].where(m["ST"].map(cnt) >= 2, "other").fillna("none"))
    samples = list(m.index)
    pt = pd.DataFrame({c: [s in set(pc[pc["cluster"] == c]["sample"]) for s in samples] for c in PLASMIDS}, index=samples)
    pt["complete_integron"] = [s in complete_int for s in samples]
    binary(f"itol/{sp}_plasmids.txt", "Plasmid clusters and integron", pt)
    gt = pd.DataFrame({g: [s in set(acq[acq["gene"] == g]["sample"]) for s in samples] for g in GENES}, index=samples)
    binary(f"itol/{sp}_genes.txt", "AMR genes (genome-level)", gt)
    print(sp, "genomes:", len(samples), "-> itol/")
