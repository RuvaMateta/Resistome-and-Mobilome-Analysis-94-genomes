import pandas as pd

MIN_COV = 0.80   # fraction of the smaller cluster covered; fix before looking at results

cols = ["qseqid","sseqid","pident","length","qlen","slen","qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_csv("plsdb/plasmid_allvsall.tsv", sep="\t", header=None, names=cols)
for side, col in (("q","qseqid"), ("s","sseqid")):
    p = df[col].str.split("|", n=2, expand=True)
    df[side+"sample"], df[side+"cl"] = p[0], p[1]

# cluster sizes (sum of contig lengths), taken before any filtering
sizes = (df.drop_duplicates("qseqid").groupby(["qsample","qcl"])["qlen"].sum())
size = sizes.to_dict()

df = df[df["qsample"] != df["ssample"]].copy()
df["lo"] = df[["qstart","qend"]].min(axis=1)
df["hi"] = df[["qstart","qend"]].max(axis=1)

def union_len(g):
    tot, cur_lo, cur_hi = 0, None, None
    for lo, hi in sorted(zip(g["lo"], g["hi"])):
        if cur_hi is None or lo > cur_hi + 1:
            if cur_hi is not None: tot += cur_hi - cur_lo + 1
            cur_lo, cur_hi = lo, hi
        else:
            cur_hi = max(cur_hi, hi)
    return tot + (cur_hi - cur_lo + 1)

per_contig = (df.groupby(["qsample","qcl","ssample","scl","qseqid"])
                .apply(union_len).rename("covered").reset_index())
pairs = per_contig.groupby(["qsample","qcl","ssample","scl"])["covered"].sum().reset_index()

df["idlen"] = df["pident"] * df["length"]
idn = df.groupby(["qsample","qcl","ssample","scl"]).agg(idlen=("idlen","sum"), alen=("length","sum")).reset_index()
pairs = pairs.merge(idn, on=["qsample","qcl","ssample","scl"])
pairs["ident"] = pairs["idlen"] / pairs["alen"]
pairs["qsize"] = [size[(a,b)] for a,b in zip(pairs["qsample"], pairs["qcl"])]
pairs["ssize"] = [size[(a,b)] for a,b in zip(pairs["ssample"], pairs["scl"])]
pairs["cov_smaller"] = (pairs["covered"] / pairs[["qsize","ssize"]].min(axis=1)).clip(upper=1.0)

sp = lambda s: s.str.rsplit("_", n=1).str[0]
pairs["qspecies"], pairs["sspecies"] = sp(pairs["qsample"]), sp(pairs["ssample"])
pairs["cross_species"] = pairs["qspecies"] != pairs["sspecies"]

# one row per unordered pair of clusters (keep the direction with most coverage)
a = pairs["qsample"] + "|" + pairs["qcl"]; b = pairs["ssample"] + "|" + pairs["scl"]
pairs["pair"] = [" ~ ".join(sorted(x)) for x in zip(a, b)]
pairs = pairs.sort_values("covered", ascending=False).drop_duplicates("pair")
pairs["pass"] = (pairs["ident"] >= 99) & (pairs["cov_smaller"] >= MIN_COV)
pairs.to_csv("merged/plasmid_pairs.tsv", sep="\t", index=False)

ok = pairs[pairs["pass"]]
print(len(pairs), "cluster pairs from different samples;", len(ok), "pass")
print(ok.groupby(["qspecies","sspecies"]).size())
