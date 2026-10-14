# Resistome and mobilome of co-circulating pathogens (Malawi, Zambia, South Africa)

Comparative genomics of 94 bacterial genomes (32 *E. coli*, 32 *K. pneumoniae*, 30 *Salmonella*) from Malawi, Zambia and South Africa, to compare resistance genes and mobile genetic elements (plasmids, integrons) across species within a One Health framework.

**Author:** Ruva Mateta · **Supervisor:** Dr Sapkota Manoj · **Institution:** Midlands State University · **Status:** Work in progress, October 2026

## What is here

| Folder | Contents |
| --- | --- |
| `00_docs_scripts_envs` | Workflow document, commands record, environment files, analysis scripts |
| `02_download_inputs` | SRA accession lists and isolate metadata |
| `03_read_QC` to `06_assembly_QC` | Read QC, trimming, assembly and CheckM2/QUAST summary reports |
| `07_resistome_AMR` | AMRFinderPlus results (merged) |
| `08_mobilome_integrons` | MOB-suite, IntegronFinder and PLSDB result tables |
| `09_HGT_phylogeny` | Plasmid comparison tables and IQ-TREE trees |

Folder numbers follow the Workflow Commands document. Every command is in `COMMANDS.md`, in the order it was run.

## Not included

Raw reads and assemblies are large and public. Download the reads with the accession lists in `02_download_inputs` (`fasterq-dump`), then follow `COMMANDS.md`.

## Sample set

90 genomes were selected. 21 were flagged by CheckM2 as incomplete or contaminated and replaced by 25 new genomes, giving 94. [Add the CheckM2 completeness and contamination thresholds used.]

## Methods in one line each

SRA Toolkit, FastQC, MultiQC, fastp, SPAdes, Flye (Galaxy), QUAST, CheckM2, Prokka, AMRFinderPlus, MOB-suite, IntegronFinder (Galaxy), BLASTn with PLSDB, Snippy, IQ-TREE, iTOL. Versions are in `00_docs_scripts_envs`.

## Reproducing

```bash
conda env create -f 00_docs_scripts_envs/qc-tools.yml
```

Then follow `COMMANDS.md`. Scripts that were run on an HPC virtual machine are in `00_docs_scripts_envs`.

## Licence and citation

Code: MIT. Cite as: Mateta R. (2026). Resistome and mobilome of co-circulating pathogens (Malawi, Zambia, South Africa). GitHub repository.
