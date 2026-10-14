#!/bin/bash
cd ~/ruva
sp=$1; f=$2; b=$(basename $f .fasta)
case $sp in
  E_coli) ref=refs/Ecoli_MG1655.fasta;;
  K_pneumoniae) ref=refs/K_pneumoniae_ERR13744124_chr.fasta;;
  Salmonella) ref=refs/Salmonella_LT2.fasta;;
esac
snippy --cpus 2 --force --outdir snippy/$b --ref $ref --ctgs $f > logs/snippy_$b.log 2>&1
