#!/bin/bash
cd ~/ruva
mkdir -p mlst_forced
run() {
  for f in input/$1_*.fasta; do
    b=$(basename "$f" .fasta)
    [ -s mlst_forced/$b.tsv ] || echo "$f"
  done | xargs -P 2 -I{} sh -c 'b=$(basename {} .fasta); mlst --scheme '"$2"' --threads 1 {} > mlst_forced/$b.tsv 2> mlst_forced/$b.err'
}
run E_coli ecoli_achtman_4
run K_pneumoniae klebsiella
run Salmonella salmonella
