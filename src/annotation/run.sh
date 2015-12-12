AnnotatePdb.groovy result/extended_pdb_ids.txt result/annotations.txt
MapMhc.groovy result/annotations.txt result/mhc.annotaitons.txt
MapTcr.groovy result/annotations.txt result/tcr.annotaitons.txt
CombineAnnotation.groovy result/annotations.txt result/mhc.annotations.txt result/final.annotations.txt