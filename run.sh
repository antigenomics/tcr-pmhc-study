#!/bin/bash

ANNOTATION="groovy annotation/src"
$ANNOTATION/AnnotatePdb.groovy result/extended_pdb_ids.txt result/annotations.txt
$ANNOTATION/MapMhc.groovy result/annotations.txt result/mhc.annotaitons.txt
$ANNOTATION/MapTcr.groovy result/annotations.txt result/tcr.annotaitons.txt
$ANNOTATION/CombineAnnotation.groovy result/annotations.txt result/mhc.annotations.txt result/final.annotations.txt

STRUCTURE="python structure/src"
$STRUCTURE/structure.py -i result/final.annotations.txt -o result/structure.txt