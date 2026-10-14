import glob, os
import pandas as pd

def merge(filename, out, split_id=False):
    frames = []
    for path in sorted(glob.glob(f"results/*/{filename}")):
        sample = os.path.basename(os.path.dirname(path))
        if os.path.getsize(path) == 0:
            continue
        df = pd.read_csv(path, sep="\t", dtype=str)
        if df.empty:
            continue
        if split_id and "sample_id" in df.columns:
            # mobtyper sample_id looks like "<sample>:<plasmid cluster>"
            df["plasmid_cluster"] = df["sample_id"].str.split(":").str[-1]
            df = df.drop(columns="sample_id")
        else:
            df = df.drop(columns="sample_id", errors="ignore")
        df.insert(0, "sample", sample)
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(out, sep="\t", index=False)
    print(f"{filename}: {merged['sample'].nunique()} samples, {len(merged)} rows -> {out}")

merge("contig_report.txt", "merged/contig_report_all.tsv")
merge("mobtyper_results.txt", "merged/mobtyper_all.tsv", split_id=True)
merge("mge.report.txt", "merged/mge_report_all.tsv")
