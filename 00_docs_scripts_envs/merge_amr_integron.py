import glob, os
import pandas as pd

def read(path, **kw):
    try:
        d = pd.read_csv(path, sep="\t", dtype=str, **kw)
    except pd.errors.EmptyDataError:
        return None
    return None if d.empty else d

# AMRFinderPlus
amr = []
for f in sorted(glob.glob("amr_results/*_amr.tsv")):
    d = read(f)
    if d is not None:
        d.insert(0, "sample", os.path.basename(f)[:-len("_amr.tsv")]); amr.append(d)
amr = pd.concat(amr, ignore_index=True).rename(columns={
    "Contig id": "contig", "Start": "start", "Stop": "stop", "Element symbol": "gene",
    "Type": "type", "Subtype": "subtype", "Class": "class", "Subclass": "subclass",
    "% Identity to reference": "pct_identity", "% Coverage of reference": "pct_coverage"})
amr = amr[["sample","contig","start","stop","gene","type","subtype","class","subclass",
           "pct_identity","pct_coverage"]]
amr["start"] = amr["start"].astype(int); amr["stop"] = amr["stop"].astype(int)

# IntegronFinder annotations (one row per element)
els = []
for f in sorted(glob.glob("integron_results/annotations/*/*.tsv")):
    d = read(f, comment="#")
    if d is not None:
        d.insert(0, "sample", os.path.basename(f)[:-4]); els.append(d)
els = pd.concat(els, ignore_index=True).rename(columns={"ID_replicon": "contig"})
els["pos_beg"] = els["pos_beg"].astype(int); els["pos_end"] = els["pos_end"].astype(int)
els.to_csv("amr_integron/integron_elements_all.tsv", sep="\t", index=False)

# IntegronFinder summary (keep only contigs that carry an integron)
sm = []
for f in sorted(glob.glob("integron_results/summary/*/*.tsv")):
    d = read(f, comment="#")
    if d is not None:
        d.insert(0, "sample", os.path.basename(f)[:-4]); sm.append(d)
sm = pd.concat(sm, ignore_index=True).rename(columns={"ID_replicon": "contig"})
for c in ("CALIN", "complete", "In0"):
    sm[c] = pd.to_numeric(sm[c])
sm = sm[(sm["CALIN"] + sm["complete"] + sm["In0"]) > 0]
sm.to_csv("amr_integron/integron_contigs_all.tsv", sep="\t", index=False)

# one span per integron, then flag AMR genes inside it
span = (els.groupby(["sample","contig","ID_integron","type"])
           .agg(lo=("pos_beg","min"), hi=("pos_end","max")).reset_index()
           .rename(columns={"type": "integron_type"}))
m = amr.merge(span, on=["sample","contig"], how="left")
gmin = m[["start","stop"]].min(axis=1); gmax = m[["start","stop"]].max(axis=1)
hits = (m[(gmin <= m["hi"]) & (gmax >= m["lo"])]
        [["sample","contig","start","stop","gene","ID_integron","integron_type"]]
        .drop_duplicates(subset=["sample","contig","start","stop","gene"]))
amr = amr.merge(hits, on=["sample","contig","start","stop","gene"], how="left")
amr["in_integron"] = amr["ID_integron"].notna()
amr.to_csv("amr_integron/amr_all.tsv", sep="\t", index=False)

print(amr["sample"].nunique(), "samples with AMRFinderPlus rows;",
      els["sample"].nunique(), "samples with an integron")
print(span.groupby("integron_type")["sample"].agg(integrons="size", samples="nunique"), "\n")
acq = amr[(amr["type"] == "AMR") & (amr["subtype"] == "AMR")]
print(len(acq), "acquired AMR gene hits;", int(acq["in_integron"].sum()), "inside an integron")
print(acq[acq["in_integron"]].groupby(["integron_type","gene"]).size()
         .sort_values(ascending=False).head(15))
