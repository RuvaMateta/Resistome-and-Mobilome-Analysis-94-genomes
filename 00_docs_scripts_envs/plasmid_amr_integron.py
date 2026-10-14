import pandas as pd

rows = []
for line in open("plsdb/plasmid_queries.fasta"):
    if line.startswith(">"):
        s, cl, contig = line[1:].strip().split()[0].split("|", 2)
        rows.append((s, cl, contig))
cmap = pd.DataFrame(rows, columns=["sample", "cluster", "contig"])

amr = pd.read_csv("merged/amr_all.tsv", sep="\t")
acq = amr[(amr["type"] == "AMR") & (amr["subtype"] == "AMR")]
els = pd.read_csv("merged/integron_elements_all.tsv", sep="\t")
ig = els[["sample", "contig", "type"]].drop_duplicates().rename(columns={"type": "integron_type"})

contigs = cmap.groupby(["sample", "cluster"]).size().rename("contigs")
genes = (cmap.merge(acq, on=["sample", "contig"]).groupby(["sample", "cluster"])["gene"]
         .agg(lambda s: ",".join(sorted(set(s)))).rename("amr_genes"))
integ = (cmap.merge(ig, on=["sample", "contig"]).groupby(["sample", "cluster"])["integron_type"]
         .agg(lambda s: ",".join(sorted(set(s)))).rename("integron"))
pc = pd.concat([contigs, genes, integ], axis=1).reset_index().fillna("-")
pc.to_csv("merged/plasmid_cluster_amr_integron.tsv", sep="\t", index=False)

print(len(pc), "plasmid clusters;", (pc["amr_genes"] != "-").sum(), "carry acquired AMR genes;",
      (pc["integron"] != "-").sum(), "carry an integron\n")

x = pd.read_csv("merged/cross_species_tiers.tsv", sep="\t")
mem = []
for _, r in x.iterrows():
    t = f"{r.qcl}~{r.scl}"
    mem += [(t, r.qsample, r.qcl), (t, r.ssample, r.scl)]
mem = pd.DataFrame(mem, columns=["type", "sample", "cluster"]).drop_duplicates()
out = mem.merge(pc, on=["sample", "cluster"], how="left").sort_values(["type", "sample"])
print(out[["type", "sample", "cluster", "contigs", "amr_genes", "integron"]].to_string(index=False))
