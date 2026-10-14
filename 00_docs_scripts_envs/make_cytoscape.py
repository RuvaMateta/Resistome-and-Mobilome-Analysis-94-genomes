import pandas as pd

meta = pd.read_csv("merged/metadata_final.csv")
mlst = pd.read_csv("mlst_forced.tsv", sep="\t", header=None, usecols=[0, 2], names=["file", "ST"], dtype=str)
mlst["sample"] = mlst["file"].str.replace("input/", "", regex=False).str.replace(".fasta", "", regex=False)
nodes = meta[["sample", "species", "country", "source_category", "year_final"]].merge(
    mlst[["sample", "ST"]], on="sample", how="left")
nodes["year_final"] = nodes["year_final"].round().astype("Int64")

# isolate - isolate edges: near-identical shared plasmid clusters (the pairs that passed the cutoffs)
p = pd.read_csv("merged/plasmid_pairs.tsv", sep="\t")
p = p[p["pass"]].copy()
m = meta.set_index("sample")
for side in ("q", "s"):
    p[side + "_country"] = p[side + "sample"].map(m["country"])
    p[side + "_year"] = p[side + "sample"].map(m["year_final"])

def tier(r):
    if r.q_country != r.s_country: return "broad regional dispersal"
    return "local co-occurrence" if r.q_year == r.s_year else "persistent reservoir"

p["tier"] = p.apply(tier, axis=1)
p["min_plasmid_size"] = p[["qsize", "ssize"]].min(axis=1)
p["plasmid_pair"] = p["qcl"] + "~" + p["scl"]
edges = p.rename(columns={"qsample": "source", "ssample": "target"})[
    ["source", "target", "plasmid_pair", "ident", "cov_smaller", "min_plasmid_size", "cross_species", "tier"]]
edges.to_csv("cytoscape/edges_plasmids.csv", index=False)

# isolate - AMR gene edges, for resistance genes found on plasmid contigs
amr = pd.read_csv("merged/amr_all.tsv", sep="\t")
acq = amr[(amr["type"] == "AMR") & (amr["subtype"] == "AMR")]
rows = []
for l in open("plsdb/plasmid_queries.fasta"):
    if l.startswith(">"):
        s, cl, c = l[1:].split()[0].split("|", 2)
        rows.append((s, cl, c))
cm = pd.DataFrame(rows, columns=["sample", "cluster", "contig"])
g = cm.merge(acq[["sample", "contig", "gene", "class"]], on=["sample", "contig"]).drop_duplicates(["sample", "gene"])
ge = g.rename(columns={"sample": "source", "gene": "target"})[["source", "target", "class"]]
ge.to_csv("cytoscape/edges_isolate_gene.csv", index=False)

gn = pd.DataFrame({"sample": sorted(ge["target"].unique())})
gn["species"] = "AMR gene"; gn["country"] = ""; gn["source_category"] = ""; gn["year_final"] = pd.NA; gn["ST"] = ""
nodes = pd.concat([nodes, gn], ignore_index=True).rename(columns={"sample": "id"})
nodes.to_csv("cytoscape/nodes.csv", index=False)

print(len(edges), "plasmid-sharing edges;", int(edges["cross_species"].sum()), "cross-species")
print(edges["tier"].value_counts().to_string())
print(len(ge), "isolate-gene edges;", ge["target"].nunique(), "genes")
