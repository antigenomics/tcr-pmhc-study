from __future__ import print_function

import argparse
import warnings
import pandas as pd
import time

from Bio.PDB import *
from shutil import rmtree
from util import *
from os import chdir, path, makedirs

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="backbone_coordinates.py\nCompute X/Y/Z coordinates of Calpha atoms of CDR loops and antigen. C terminal residue is used as origin, coordinates are rotated so that C-N terminus direction is aligned with X axis and YZ projection of the center of mass is aligned with Z axis.")
parser.add_argument("-i", nargs=1, type=str, default="../../output/final.annotations.txt",
                    help="Annotation table path.")
parser.add_argument("-o", nargs=1, type=str, default="../../output/backbone.txt",
                    help="Output table path.")
parser.add_argument("-o2", nargs=1, type=str, default="../../output/backbone_ag.txt",
                    help="Output table path for antigen.")
parser.add_argument("-t", nargs=1, type=str, default="../tmp/",
                    help="Temporary folder path")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_file = path.abspath(args.i)

if type(args.o) is list:
    args.o = args.o[0]
output_file = path.abspath(args.o)

if type(args.o2) is list:
    args.o2 = args.o[0]
output_file2 = path.abspath(args.o2)

chdir(path.dirname(path.realpath(__file__)))

# Create temporary folders
pdb_dir = path.abspath(args.t + '/pdb/')

if not path.exists(pdb_dir):
    makedirs(pdb_dir)

param_template_path = path.abspath("../res/")

# Writing output file header
col_names = ['pdb_id', 'species',
             'mhc_type', 'mhc_a_allele', 'mhc_b_allele',
             'tcr_v_allele', 'tcr_region', 'tcr_region_seq',
             'aa_tcr', 'len_tcr', 'pos_tcr', 
             'x', 'y', 'z',
             'x_cb', 'y_cb', 'z_cb']

with open(output_file, 'w') as f:
    f.write('\t'.join(col_names) + '\n')

# Main loop

pdb_list = PDBList()
pdb_parser = PDBParser()

table = pd.read_table(input_file)
table = table[table['tcr_region'].str.contains("CDR")]

bypdb = table.groupby("pdb_id")

i = 0

results_ag = []

for pdb_id, pdb_group in bypdb:
    # Load PDB file
    pdb_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=pdb_dir)

    print("[", time.strftime("%c"), "|", i, "/", table.shape[0], "]")
    print(pdb_id, "- preparing for computation")

    # Load model from original pdb file
    model = pdb_parser.get_structure(pdb_id, pdb_file)[0]

    # Store annotation for entire complex
    pdb_annot = pdb_group.iloc[0]

    # Get and check antigen residues
    antigen_chain = model[pdb_annot['chain_antigen']]
    antigen_seq = pdb_annot['antigen_seq']
    antigen_range = range(len(antigen_seq))
    antigen_residues = get_residues(antigen_chain, antigen_range)
    antigen_seq_obs = get_seq(antigen_residues)

    if antigen_seq != antigen_seq_obs:
        warning("Antigen sequence mismatch (expected observed): ", antigen_seq, antigen_seq_obs,
                ". Replacing with one from PDB.")
        pdb_annot['antigen_seq'] = antigen_seq_obs

    print(pdb_id, "- computing positions for", pdb_annot['mhc_type'], ":", antigen_seq_obs)

    # Compute positions and add them to results
    positions = calc_coords_ag(pdb_id, 
        pdb_annot['mhc_type'],
        antigen_residues)

    results_ag.extend(positions)

    # Iterate by TCR chain/region (CDR1,2,3)
    byregion = pdb_group.groupby(['chain_tcr', 'tcr_region'])

    results_by_pdb = []
    for tcr_region_id, tcr_region_group in byregion:
        # Get and check tcr region residues
        tcr_annot = tcr_region_group.iloc[0]

        tcr_chain = model[tcr_annot['chain_tcr']]
        tcr_region_seq = tcr_annot['tcr_region_seq']
        tcr_region_range = range(tcr_annot['tcr_region_start'], tcr_annot['tcr_region_end'])
        tcr_region_residues = get_residues(tcr_chain, tcr_region_range)
        tcr_region_seq_obs = get_seq(tcr_region_residues)

        if tcr_region_seq != tcr_region_seq_obs:
            warning("TCR:", tcr_region_id, " sequence mismatch (expected observed): ", tcr_region_seq,
                    tcr_region_seq_obs, ". Replacing with one from PDB.")
            tcr_annot['tcr_region_seq'] = tcr_region_seq_obs

        print(pdb_id, "- computing positions for", tcr_annot['tcr_v_allele'], ':', tcr_region_id[1])

        # Compute positions and add them to results
        positions = calc_coords(tcr_annot['tcr_v_allele'][0:2], tcr_annot['tcr_region'],
            tcr_region_residues)

        for row in positions:
            row.update(tcr_annot.to_dict())

        results_by_pdb.extend(positions)

        i += 1

    # Write selected columns
    res = pd.DataFrame(results_by_pdb)[col_names].drop_duplicates()
    res.to_csv(output_file, sep='\t', header=False, index=False, mode='a')
    print("Done")

pd.DataFrame(results_ag).to_csv(output_file2, sep='\t', index=False)

print("Finished processing", table.shape[0], " entries.")
