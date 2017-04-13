#!/bin/bash

mkdir logs/
mkdir output/
mkdir annotation/tmp/

export PDB_DIR="`pwd`/annotation/tmp/"

cd annotation/
groovy -cp . AnnotatePdb.groovy ../input/good_pdb_ids.txt ../output/annotations.txt 2>&1 | tee ../logs/annotate_pdb.log
groovy -cp . MapMhc.groovy ../output/annotations.txt ../output/mhc.annotations.txt 2>&1 | tee ../logs/mapmhc.log
#groovy -cp . MapTcr.groovy ../output/annotations.txt ../output/tcr.annotations.txt 2>&1 | tee ../logs/maptcr.log
#groovy -cp . CombineAnnotation.groovy ../output/annotations.txt ../output/mhc.annotations.txt ../output/tcr.annotations.txt ../output/final.annotations.txt 2>&1 | tee ../logs/combine.log
cd ..

#STRUCTURE="python structure/src"
#$STRUCTURE/structure.py -i result/final.annotations.txt -o result/structure.txt