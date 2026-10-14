# Comparative Genomics Workflow: Command Record

Record of the commands behind the comparative genomics analysis of 94 bacterial genomes (*E. coli*, *K. pneumoniae*, *Salmonella*; Malawi, Zambia, South Africa), from read download to plasmid and phylogeny analysis.

| Stage | Source |
| --- | --- |
| Batch 1 download, read QC, trimming, SPAdes (original 90 genomes) | Workflow Commands document |
| Batch 2 download to assembly (25 replacement genomes) | Shell history, WSL workstation |
| CheckM2, QUAST, MultiQC, Prokka, AMRFinderPlus, input prep, merge scripts | Shell history, WSL workstation |
| Flye assembly, IntegronFinder | Galaxy server (not in shell history) |
| MOB-suite, PLSDB BLASTn, plasmid comparison, Snippy, IQ-TREE | Shell history, HPC virtual machine |

Paths are relative to the home directory on the WSL workstation unless stated. Assembly work is under `~/project1/My_Assemblies`.

Placeholders: `<port>` = SSH port for the HPC VM; `<user>@<hpc>` = your HPC login; `<accessions>` = accession list for the loop.

## Sample set

90 genomes selected. 21 failed CheckM2 QC and were replaced by 25, giving 94 analysed (32 *E. coli*, 32 *K. pneumoniae*, 30 *Salmonella*).

| Group | Failed QC | Replacements | Replacement accessions |
| --- | --- | --- | --- |
| South Africa *E. coli* | 10 | 10 | SRR29213288, SRR29213289, SRR29213310, SRR29213325-29, SRR29213331, SRR29213334 |
| Zambia *E. coli* | 1 | 3 | SRR31865353, SRR31865355, SRR31865377 |
| Zambia *K. pneumoniae* | 10 | 12 | ERR15165589-600 (excluding ERR15165596, ERR15165601) |
| Total | 21 | 25 | |

Full accession lists are in `02_download_inputs/Acc_list/`.

## Software environments

Each tool group has its own conda environment. Exports are in `00_docs_scripts_envs/`.

| Environment | Tools | Used for |
| --- | --- | --- |
| qc-tools | FastQC, fastp, MultiQC, NanoPlot, Filtlong, NanoFilt, pandas | Read QC, trimming, reports, table merging |
| assembly | SPAdes, QUAST | Illumina assembly and assembly stats |
| checkm2 | CheckM2 | Completeness and contamination |
| annotation | Prokka | Genome annotation |
| genomic_profile | AMRFinderPlus | Resistance gene detection |
| mobilome (HPC) | MOB-suite | Plasmid reconstruction and typing |
| phylo (HPC) | Snippy, IQ-TREE, mlst | Core SNPs, trees, sequence typing |
| netviz (HPC) | networkx, matplotlib | Plasmid-sharing network |

Creation commands:

```bash
conda env create -f ~/qc-tools.yml
conda create -n assembly -c bioconda -c conda-forge spades
conda install -n assembly -c bioconda -c conda-forge quast
mamba create -n checkm2 python=3.8 -y
mamba install -n checkm2 -c bioconda -c conda-forge "checkm2=1.0.2" -y
conda create -n annotation -c bioconda -c conda-forge prokka
conda create -n genomic_profile -c bioconda ncbi-amrfinderplus

# HPC virtual machine
conda create -n mobilome --override-channels -c conda-forge -c bioconda mob_suite -y
conda create -n phylo --override-channels -c conda-forge -c bioconda snippy iqtree mlst -y
conda create -n netviz --override-channels -c conda-forge python=3.11 pandas networkx matplotlib -y
```

MOB-suite and IntegronFinder would not install cleanly alongside the other tools on the workstation (pandas/NumPy conflicts), so MOB-suite ran on the HPC VM and IntegronFinder on Galaxy.

## Download and raw-read QC

Reads were downloaded from the NCBI SRA, one accession at a time, into one folder per species and country:

```bash
for acc in <accessions>; do
  echo -e "\n== Processing ${acc} =="
  fasterq-dump --split-3 -p -e 4 "$acc" -O ./<species>_<country>_fastq/
  gzip -v ./<species>_<country>_fastq/${acc}*.fastq
done
```

Raw-read QC (Illumina):

```bash
conda activate qc-tools
mkdir -p QC_ILLUMINA QC_ILLUMINA_REPORT
for d in Ecoli_*fastq Salmonella_*fastq K_pneumoniae_SouthAfrica_fastq K_pneumoniae_Zambia_fastq; do
  fastqc "$d"/*.fastq.gz -o QC_ILLUMINA
done
multiqc QC_ILLUMINA -o QC_ILLUMINA_REPORT
```

Batch 2 used the same commands with `QC_ILLUMINA2` and `QC_ILLUMINA_REPORT2`.

PacBio reads (Malawi *K. pneumoniae*) checked with NanoPlot:

```bash
mkdir QC_PACBIO
for f in K_pneumoniae_Malawi_fastq/*.fastq.gz; do
  NanoPlot --fastq "$f" -o QC_PACBIO/$(basename "$f" .fastq.gz)
done
```

## Read trimming

Illumina pairs trimmed with fastp (paired-end adapter detection, sliding-window Q20, min length 50 bp):

```bash
mkdir -p QC_ILLUMINA_TRIM reports_fastp
for fq1 in */*_1.fastq.gz; do
  fq2=${fq1/_1/_2}
  acc=$(basename "$fq1" _1.fastq.gz)
  dir=$(basename $(dirname "$fq1"))
  fastp -i "$fq1" -I "$fq2" \
    -o QC_ILLUMINA_TRIM/${acc}_1.trimmed.fastq.gz \
    -o QC_ILLUMINA_TRIM/${acc}_2.trimmed.fastq.gz \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 50 --thread 8 \
    --html reports_fastp/${dir}_${acc}.html \
    --json reports_fastp/${dir}_${acc}.json
done
fastqc QC_ILLUMINA_TRIM/*.fastq.gz -o reports_fastp/
multiqc reports_fastp/ -o QC_ILLUMINA_REPORT/
```

Batch 2: same commands, `QC_ILLUMINA_TRIM2` and `reports_fastp2`.

PacBio: Filtlong (min 1,000 bp, keep 90 %) then NanoFilt (Q20), then NanoPlot:

```bash
mkdir -p QC_PACBIO_TRIM reports_pacbio
for fq in K_pneumoniae_Malawi_fastq/*.fastq.gz; do
  acc=$(basename "$fq" .fastq.gz)
  filtlong --min_length 1000 --keep_percent 90 "$fq" | NanoFilt -q 20 | gzip \
    > QC_PACBIO_TRIM/${acc}_filtered_trimmed.fastq.gz
  NanoPlot --fastq QC_PACBIO_TRIM/${acc}_filtered_trimmed.fastq.gz \
    -o reports_pacbio/${acc}_NanoPlot --threads 8
done
multiqc reports_pacbio/ -o reports_pacbio/
```

## Assembly

Illumina genomes assembled with SPAdes in isolate mode. Batch 1: 4 threads, 4 GB. Batch 2: 2 threads, 2 GB. Six SA *E. coli* genomes re-run at 4 threads, 5 GB after OOM errors.

```bash
conda activate assembly
for R1 in *1.trimmed.fastq.gz; do
  R2=${R1/1.trimmed.fastq.gz/2.trimmed.fastq.gz}
  acc=$(basename "$R1" _1.trimmed.fastq.gz)
  spades.py -1 "$R1" -2 "$R2" -o spades_results2/$acc -t 2 -m 2 --isolate \
    | tee spades_results2/$acc/run.log
done

for a in SRR29213325 SRR29213326 SRR29213327 SRR29213328 SRR29213331 SRR29213334; do
  spades.py --isolate \
    -1 QC_ILLUMINA_TRIM2/${a}_1.trimmed.fastq.gz \
    -2 QC_ILLUMINA_TRIM2/${a}_2.trimmed.fastq.gz \
    -o spades_results2_rerun/$a -t 4 -m 5 || echo "FAILED: $a" >> spades_failed.txt
done
```

Batches merged, then sorted per species:

```bash
rsync -av --ignore-existing spades_results2/ spades_results1/
rsync -av --ignore-existing spades_results2_rerun/ spades_results1/
while read acc; do
  [ -d "spades_results/$acc" ] && mv "spades_results/$acc" grouped/E_coli/
done < E_coli_accessions_clean.txt
```

The same sorting loop was run for *K. pneumoniae* and *Salmonella* with their accession lists.

The 10 PacBio *K. pneumoniae* genomes were assembled with Flye on Galaxy.

## Assembly QC

QUAST per species; CheckM2 per genome; MultiQC combined.

```bash
conda activate assembly
for org in E_coli Salmonella K_pneumoniae; do
  quast grouped/$org/*/contigs.fasta -o ../../quast_results/$org
done
```

CheckM2 database downloaded from Zenodo into `~/project1/checkm2_db`; run one genome at a time (single thread, low memory):

```bash
for asm in "${INPUT_DIR}"/*/contigs.fasta; do
  sample=$(basename "$(dirname "$asm")")
  checkm2 predict --input "$asm" --output-directory "$OUTPUT_DIR/$sample" -x fasta \
    --database_path ~/project1/checkm2_db/CheckM2_database/uniref100.K0.1.dmnd \
    --threads 1 --lowmem --force
done

multiqc quast_results CheckM2_all -o ../multiqc_report         # batch 1
multiqc quast_results1 checkm2_results1 -o multiqc_report_25samples --force  # batch 2
```

## Annotation

```bash
conda activate annotation
for org in E_coli Salmonella K_pneumoniae; do
  mkdir -p ~/project1/My_Assemblies/prokka_results/$org
  for asm in grouped/$org/*/scaffolds.fasta; do
    acc=$(basename "$(dirname "$asm")")
    prokka --outdir ~/project1/My_Assemblies/prokka_results/$org/$acc --prefix $acc "$asm"
  done
done
```

## Resistome profiling

```bash
conda activate genomic_profile
amrfinder -u
mkdir -p ~/project1/My_Assemblies/amr_results
declare -A ORG
ORG[E_coli]="Escherichia"
ORG[K_pneumoniae]="Klebsiella_pneumoniae"
ORG[Salmonella]="Salmonella"

for species in E_coli K_pneumoniae Salmonella; do
  for acc in ~/project1/My_Assemblies/spades_results/grouped/$species/*; do
    acc_name=$(basename "$acc")
    asm="$acc/contigs.fasta"
    [ -f "$asm" ] || asm="$acc/scaffolds.fasta"
    [ -f "$asm" ] || { echo "WARNING: no assembly for $species/$acc_name"; continue; }
    amrfinder -n "$asm" -O "${ORG[$species]}" --plus \
      -o ~/project1/My_Assemblies/amr_results/${species}_${acc_name}_amr.tsv
  done
done
```

94 per-genome files merged by `merge_amr_integron.py` into `amr_all.tsv`.

## Mobilome and integrons

Species-prefixed FASTA copies for consistent sample names:

```bash
mkdir -p integron_input
for d in spades_results/grouped/*/*; do
  acc=$(basename "$d")
  sp=$(basename "$(dirname "$d")")
  if [ -f "$d/contigs.fasta" ]; then
    cp "$d/contigs.fasta" integron_input/${sp}_${acc}.fasta
  else
    echo "MISSING contigs.fasta: $sp $acc"
  fi
done
```

IntegronFinder on the 94 FASTA files via Galaxy.

MOB-suite on the HPC VM, 4 genomes in parallel, 4 threads each:

```bash
conda activate mobilome

mob_recon --infile input/E_coli_SRR29213288.fasta \
          --outdir results/E_coli_SRR29213288 --num_threads 4

nohup sh -c "ls input/*.fasta | xargs -P 4 -I{} sh -c 'S=\$(basename {} .fasta); mob_recon --infile {} --outdir results/\$S --num_threads 4 > logs/\$S.log 2>&1'" > logs/run_all.log 2>&1 &

python merge_mobsuite.py
```

PLSDB BLASTn search:

```bash
makeblastdb -in plsdb.fasta -dbtype nucl -out plsdb_db -title PLSDB

nohup blastn -query plsdb/plasmid_queries.fasta -db plsdb/plsdb_db \
  -evalue 1e-20 -max_target_seqs 50 -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore stitle" \
  -out plsdb/plsdb_hits.tsv > logs/plsdb_blast.log 2>&1 &

python plsdb_best_hits.py
```

## Metadata and merged tables

```bash
python merge_metadata.py
python join_metadata.py
python merge_amr_integron.py
```

Outputs: `metadata_all.csv`, `metadata_final.csv`, `amr_all.tsv`, `integron_elements_all.tsv`, `integron_contigs_all.tsv`, `combined_summary.tsv`.

The 10 Malawi *K. pneumoniae* genomes (BioProject PRJEB64289) have no recorded isolation source in the run metadata, so their source category is set to `community_unspecified`.

## Plasmid comparison, typing and phylogeny

Plasmid comparison (all vs all, 99 % identity, ≥80 % coverage of the smaller plasmid; 40 cross-species pairs collapsed to 6 plasmid types):

```bash
makeblastdb -in plasmid_queries.fasta -dbtype nucl -out plasmid_self_db

nohup blastn -query plasmid_queries.fasta -db plasmid_self_db \
  -evalue 1e-20 -perc_identity 99 -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  -out plasmid_allvsall.tsv > ~/ruva/logs/allvsall.log 2>&1 &

python plasmid_pairs.py
python tiers.py
python plasmid_amr_integron.py
python shared_region_genes.py
```

Sequence typing:

```bash
conda activate phylo
nohup ./run_mlst_forced.sh > logs/mlst_forced.log 2>&1 &
```

Core-genome SNPs (Snippy on assemblies, `--ctgs`, one reference per species):

```bash
conda activate phylo
nohup ./run_snippy_all.sh > logs/snippy_all.log 2>&1 &

snippy-core --ref refs/Ecoli_MG1655.fasta --prefix core/E_coli snippy/E_coli_*
snippy-core --ref refs/Salmonella_LT2.fasta --prefix core/Salmonella snippy/Salmonella_*
snippy-core --ref refs/K_pneumoniae_ERR13744124_chr.fasta \
  --prefix core/K_pneumoniae_no119 \
  $(ls -d snippy/K_pneumoniae_* | grep -v ERR13744119)
```

ERR13744119 was left out of the *K. pneumoniae* core alignment.

Trees (GTR+ASC+G, 1,000 ultrafast bootstraps):

```bash
nohup sh -c 'for sp in E_coli K_pneumoniae_no119 Salmonella; do
  iqtree -s core/$sp.aln -m GTR+ASC+G -B 1000 -T 4 --prefix tree/$sp
done' > logs/iqtree.log 2>&1 &
```

Annotation and network:

```bash
conda activate mobilome
python make_itol.py
python make_cytoscape.py

conda activate netviz
python draw_network.py
```

Trees copied back with rsync:

```bash
rsync -avP -e "ssh -p <port>" \
  --include='*.treefile' --include='*.iqtree' --exclude='*' \
  <user>@<hpc>:~/ruva/tree/ ~/project1/hpc_results/tree/
```

## Still to attach

See `STILL_TO_ADD.md` at the repository root.
