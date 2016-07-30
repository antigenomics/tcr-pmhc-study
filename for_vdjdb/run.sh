#!/bin/bash

groovy MapJreg.groovy ../tmp/for_J_reg_annot.txt ../../result/tcr.jreg.annotations.txt
python merge_j_with_final.py
python find_refs.py