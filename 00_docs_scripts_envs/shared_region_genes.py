import pandas as pd

cols = ["qseqid","sseqid","pident","length","qlen","slen","qstart","qend","sstart","send","evalue","bitscore"]
h = pd.read_csv("plsdb/plasmid_allvsall.tsv", sep="\t", header=None, names=cols)
for side, c in (("q", "qseqid"), ("s", "sseqid")):
    p = h[c].str.split("|", n=2, expand=True)
    h[side+"sample"], h[side+"cl"], h[side+"contig"] = p[0], p[1], p[2]
pairs = pd.read_csv("merged/cross_species_tiers.tsv", sep="\t")[["qsample","qcl","ssample","scl"]].drop_duplicates()
h = h.merge(pairs, on=["qsample","qcl","ssample","scl"])

amr = pd.read_csv("merged/amr_all.tsv", sep="\t")
acq = amr[(amr["type"] == "AMR") & (amr["subtype"] == "AMR")][["sample","contig","start","stop","gene","in_integron"]]

res = []
for side in ("q", "s"):
    A = acq.rename(columns={"sample": side+"sample", "contig": side+"contig"})
    m = h.merge(A, on=[side+"sample", side+"contig"])
    lo = m[[side+"start", side+"end"]].min(axis=1); hi = m[[side+"start", side+"end"]].max(axis=1)
    g1 = m[["start","stop"]].min(axis=1); g2 = m[["start","stop"]].max(axis=1)
    m = m[(g1 >= lo) & (g2 <= hi)].copy()
    m["side"], m["cluster"], m["isolate"] = side, m[side+"cl"], m[side+"sample"]
    res.append(m[["qcl","scl","side","cluster","isolate","gene","in_integron"]])
r = pd.concat(res)
out = r.groupby(["qcl","scl","side","cluster"]).agg(
    isolates=("isolate","nunique"),
    genes_in_shared_region=("gene", lambda s: ",".join(sorted(set(s)))),
    any_in_integron=("in_integron","any")).reset_index()
out.to_csv("merged/shared_region_genes.tsv", sep="\t", index=False)
pd.set_option("display.max_colwidth", 120)
print(out.to_string(index=False))
