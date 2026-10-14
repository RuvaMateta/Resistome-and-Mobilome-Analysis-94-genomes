import pandas as pd

MIN_IDENT = 99.0     # percent identity cutoff (your step-10 threshold)
MIN_COV = 0.80       # fraction of the plasmid cluster covered; set before looking at results

cols = ["qseqid","sseqid","pident","length","qlen","slen","qcovhsp","evalue","bitscore","stitle"]
df = pd.read_csv("plsdb/plsdb_hits.tsv", sep="\t", header=None, names=cols)
p = df["qseqid"].str.split("|", n=2, expand=True)
df["sample"], df["cluster"], df["contig"] = p[0], p[1], p[2]

# 1. best hit per contig (highest bitscore)
best_contig = (df.sort_values("bitscore", ascending=False)
                 .groupby("qseqid", as_index=False).head(1))
best_contig.to_csv("plsdb/plsdb_best_per_contig.tsv", sep="\t", index=False)

# 2. per cluster: combine hits to the same PLSDB plasmid across the cluster's contigs
sizes = (df.drop_duplicates("qseqid")
           .groupby(["sample","cluster"])["qlen"].sum().rename("cluster_size"))
df["id_x_len"] = df["pident"] * df["length"]
g = (df.groupby(["sample","cluster","sseqid"])
       .agg(aligned=("length","sum"), id_x_len=("id_x_len","sum"),
            stitle=("stitle","first")).reset_index())
g["ident"] = g["id_x_len"] / g["aligned"]
g = g.merge(sizes.reset_index(), on=["sample","cluster"])
g["cov"] = (g["aligned"] / g["cluster_size"]).clip(upper=1.0)
g = g.sort_values(["cov","ident"], ascending=False)
best_cl = g.groupby(["sample","cluster"], as_index=False).head(1).copy()
best_cl["pass"] = (best_cl["ident"] >= MIN_IDENT) & (best_cl["cov"] >= MIN_COV)
best_cl = best_cl.drop(columns="id_x_len")
best_cl.to_csv("plsdb/plsdb_best_per_cluster.tsv", sep="\t", index=False)

best_cl["species"] = best_cl["sample"].str.rsplit("_", n=1).str[0]
print(len(best_cl), "plasmid clusters with a PLSDB hit")
print(best_cl.groupby("species")["pass"].agg(clusters="size", passing="sum"))
