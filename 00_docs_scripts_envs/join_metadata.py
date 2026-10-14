import re
import pandas as pd

meta = pd.read_csv("metadata_all.csv", dtype=str)
lk = pd.read_csv("bioproject_lookup.csv", dtype=str)
m = meta.merge(lk, on="bioproject", how="left", suffixes=("", "_lk"))

ENV = r"water|river|borehole|tap|dam\b|spring|septic|sewage|wastewater|drain|soil|effluent"
ANI = r"chicken|poultry|fowl|pig|porcine|lamb|sheep|cattle|cow|beef|goat|meat|food|carcass|lobster|egg|airsac|farm|animal"
CLI = r"human|patient|blood|urine|stool|clinical|hospital|sputum"
text = (m[["isolation_source","sample_type","host","env_medium"]].fillna("")
          .agg(" ".join, axis=1).str.lower())

def classify(t):
    if re.search(ENV, t): return "environmental"
    if re.search(ANI, t): return "animal_food"
    if re.search(CLI, t): return "clinical"
    return None

m["source_text_class"] = text.map(classify)
m["source_category"] = m["source_category"].where(
    ~m["source_category"].isin(["mixed","check"]) & m["source_category"].notna(),
    m["source_text_class"])
m["source_category"] = m["source_category"].fillna(m["source_text_class"])
m["source_category"] = m["source_category"].fillna("needs_review")

m["year_iso"] = pd.to_numeric(m["year"], errors="coerce")
ymin = pd.to_numeric(m["year_min"], errors="coerce")
ymax = pd.to_numeric(m["year_max"], errors="coerce")
m["year_final"] = m["year_iso"].fillna(ymin.where(ymin == ymax))
m["year_range"] = (ymin.astype("Int64").astype(str) + "-" + ymax.astype("Int64").astype(str))
m["year_basis"] = m["year_iso"].notna().map({True: "isolate", False: "bioproject"})
m["country"] = m["country_file"]

out = m[["sample","species","country","bioproject","source_category","year_final",
         "year_range","year_basis","isolation_source","sample_type","host","env_medium"]]
out.to_csv("metadata_final.csv", index=False)

print(len(out), "isolates")
print(pd.crosstab([out["country"], out["species"]], out["source_category"]), "\n")
print("year from bioproject range only (no single year):",
      (out["year_final"].isna()).sum(), "isolates")
print("needs review:"); print(out[out["source_category"]=="needs_review"][["sample","bioproject"]].to_string(index=False))
