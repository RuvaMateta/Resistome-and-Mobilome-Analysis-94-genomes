import pandas as pd

x = pd.read_csv("merged/cross_species_pairs.tsv", sep="\t")
meta = pd.read_csv("merged/metadata_final.csv")[["sample","country","year_final","source_category"]]
for side in ("q", "s"):
    mm = meta.rename(columns={"sample": side+"sample", "country": side+"_country",
                              "year_final": side+"_year", "source_category": side+"_source"})
    x = x.merge(mm, on=side+"sample", how="left")

def tier(r):
    if r.q_country != r.s_country: return "broad regional dispersal"
    if pd.isna(r.q_year) or pd.isna(r.s_year): return "year unknown"
    return "local co-occurrence" if r.q_year == r.s_year else "persistent reservoir"

x["tier"] = x.apply(tier, axis=1)
x.to_csv("merged/cross_species_tiers.tsv", sep="\t", index=False)
print(x["tier"].value_counts(), "\n")

rows = []
for (a, b), g in x.groupby(["qcl", "scl"]):
    rows.append({"q_plasmid": a, "s_plasmid": b, "pairs": len(g),
        "tiers": ", ".join(sorted(set(g["tier"]))),
        "countries": ", ".join(sorted(set(g["q_country"]) | set(g["s_country"]))),
        "years": ", ".join(str(int(v)) for v in sorted(set(g["q_year"].dropna()) | set(g["s_year"].dropna()))),
        "sources": ", ".join(sorted(set(g["q_source"].dropna()) | set(g["s_source"].dropna())))})
print(pd.DataFrame(rows).to_string(index=False))
