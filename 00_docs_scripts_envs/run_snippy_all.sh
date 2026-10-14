#!/bin/bash
cd ~/ruva
for sp in E_coli K_pneumoniae Salmonella; do
  for f in input/${sp}_*.fasta; do
    b=$(basename $f .fasta)
    [ "$b" = "K_pneumoniae_ERR13744124" ] && continue
    [ -s snippy/$b/snps.tab ] && continue
    echo "$sp $f"
  done
done | xargs -P 3 -L 1 ./snippy_one.sh
