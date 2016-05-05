#!/bin/bash

cd annotation/src
#groovy MapMhc.groovy ../../result/annotations.txt ../../result/mhc.annotations.txt
#groovy MapTcr.groovy ../../result/annotations.txt ../../result/tcr.annotations.txt
#groovy MapJreg.groovy ../tmp/for_J_reg_annot.txt ../../result/tcr.jreg.annotations.txt
#groovy CombineAnnotation.groovy ../../result/annotations.txt ../../result/mhc.annotations.txt ../../result/final.annotations.txt
#python merge_j_with_final.py
python find_refs.py
cd -

#STRUCTURE="python structure/src"
#$STRUCTURE/pdbstat.py result/final.annotations.txt result/structure.txt
