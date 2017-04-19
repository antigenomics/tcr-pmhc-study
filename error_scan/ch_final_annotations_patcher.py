import os
import glob
import pandas as pd
import copy
import Bio.PDB
import time
from antigenomics_util import *
import argparse

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will use annotations table patch to fix errors in final annotations table.")
parser.add_argument("-i", nargs=1, type=str, default="../result/final.annotations.patch.txt",
                    help="Annotation table patch path.")
parser.add_argument("-o", nargs=1, type=str, default="../result/final.annotations.txt",
                    help="Annotation table path.")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
patch_file = os.path.abspath(args.i)

if type(args.o) is list:
    args.o = args.o[0]
annotation_file = os.path.abspath(args.o)

patchdata = pd.read_csv(patch_file, sep='\t', dtype=str)
inpdata = pd.read_csv(annotation_file, sep='\t', dtype=str)

for pdb in patchdata.groupby('pdb_id'):
    inpdata = inpdata[inpdata['pdb_id'] != pdb[0]]
    inpdata = inpdata.append(patchdata[patchdata['pdb_id'] == pdb[0]])

inpdata.sort_values(by='pdb_id').to_csv(annotation_file, sep='\t', index=None)