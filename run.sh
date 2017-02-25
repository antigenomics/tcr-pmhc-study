#!/bin/bash

ANNOTATION="groovy src"

#cd annotation/src/

#rm -r ../tmp
#mkdir ../tmp

#groovy AnnotatePdb.groovy ../../input/extended_pdb_ids.txt ../../result/annotations.txt 2>&1 | tee ../tmp/annotate_pdb.log
#groovy MapMhc.groovy ../../result/annotations.txt ../../result/mhc.annotations.txt 2>&1 | tee ../tmp/mapmhc.log
#groovy MapTcr.groovy ../../result/annotations.txt ../../result/tcr.annotations.txt 2>&1 | tee ../tmp/maptcr.log
#groovy CombineAnnotation.groovy ../../result/annotations.txt ../../result/mhc.annotations.txt ../../result/tcr.annotations.txt ../../result/final.annotations.txt 2>&1 | tee ../tmp/combine.log
#cd ..

STRUCTURE="python structure/src"
$STRUCTURE/structure.py -i result/final.annotations.txt -o result/structure.txt