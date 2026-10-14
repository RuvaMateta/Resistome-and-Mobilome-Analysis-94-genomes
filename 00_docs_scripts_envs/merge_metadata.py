import glob, re
import pandas as pd

MISSING = {"", "missing", "not collected", "not applicable", "na", "n/a",
           "unknown", "not available", "not provided", "none", "null"}
SPECIES = {"Ecoli": "E_coli", "K_pneumoniae": "K_pneumoniae", "Salmonella": "Salmonella"}

def pick(df, names):
    out = pd.Series([None] * len(df), index=df.index, dtype="object")
    for n in names:
        for c in df.columns:
            if c.lower() == n.lower():
                v = df[c].astype("object")
                v = v.where(~v.astype(str).str.strip().str.lower().isin(MISSING), None)
                out = out.where(out.notna(), v)
    return out

files = sorted(glob.glob("Malawi_*.csv") + glob.glob("SouthAfrica_*.csv") + glob.glob("Zambia_*.csv"))
frames = []
for f in files:
    country, sp = f[:-4].split("_", 1)
    df = pd.read_csv(f, dtype=str)
    m = pd.DataFrame({"sample": SPECIES[sp] + "_" + df["Run"].str.strip()})
    m["species"] = SPECIES[sp]
    m["country_file"] = country
    m["country_meta"] = pick(df, ["geo_loc_name_country"])
    m["collection_date"] = pick(df, ["Collection_Date"])
    m["isolation_source"] = pick(df, ["isolation_source"])
    m["sample_type"] = pick(df, ["sample_type"])
    m["host"] = pick(df, ["HOST", "host_scientific_name"])
    m["env_medium"] = pick(df, ["env_medium"])
    m["bioproject"] = pick(df, ["BioProject"])
    m["biosample"] = pick(df, ["BioSample"])
    m["lat_lon"] = pick(df, ["lat_lon"])
    frames.append(m)

meta = pd.concat(frames, ignore_index=True)
meta["year"] = meta["collection_date"].str.extract(r"((?:19|20)\d{2})")[0]
meta.to_csv("metadata_all.csv", index=False)

print(len(meta), "isolates;", meta["sample"].nunique(), "unique sample names")
print(meta.groupby(["country_file", "species"]).size().unstack(fill_value=0), "\n")
src = ["isolation_source", "sample_type", "host", "env_medium"]
meta["has_source"] = meta[src].notna().any(axis=1)
print("missing year / no source info, by country and species:")
print(meta.groupby(["country_file", "species"]).agg(
    n=("sample", "size"), no_year=("year", lambda s: s.isna().sum()),
    no_source=("has_source", lambda s: (~s).sum())), "\n")
print("country as written in the metadata vs. file name:")
print(pd.crosstab(meta["country_file"], meta["country_meta"].fillna("(blank)")))
