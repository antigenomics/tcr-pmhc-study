#!/bin/bash

command -v blastp >/dev/null 2>&1 || { echo >&2 "BLASTP not installed. Aborting."; exit 1; }
command -v igblastp >/dev/null 2>&1 || { echo >&2 "IGBLASTP not installed. Aborting."; exit 1; }

mkdir logs/
mkdir output/
mkdir annotation/tmp/

export PDB_DIR="`pwd`/annotation/tmp/"

#cd annotation/
#groovy -cp . AnnotatePdb.groovy ../input/good_pdb_ids.txt ../output/annotations.txt 2>&1 | tee ../logs/annotate_pdb.log
#groovy -cp . MapMhc.groovy ../output/annotations.txt ../output/mhc.annotations.txt 2>&1 | tee ../logs/mapmhc.log
#groovy -cp . MapTcr.groovy ../output/annotations.txt ../output/tcr.annotations.txt 2>&1 | tee ../logs/maptcr.log
#groovy -cp . CombineAnnotation.groovy ../output/annotations.txt ../output/mhc.annotations.txt ../output/tcr.annotations.txt ../output/final.annotations.txt 2>&1 | tee ../logs/combine.log
#cd ..

python structure/src/structure.py -i output/final.annotations.txt -o output/structure.txt