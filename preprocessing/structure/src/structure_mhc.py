from __future__ import print_function

import argparse
import warnings
import pandas as pd

from Bio.PDB import *
from shutil import rmtree
from util import *
from os import chdir, path, makedirs

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="structure.py")
parser.add_argument("-i", nargs=1, type=str, default="../../result/final.annotations.txt",
                    help="annotation table ptah")
parser.add_argument("-o", nargs=1, type=str, default="../../result/structure.mhc.txt",
                    help="output table path")
parser.add_argument("-d", "--distance", nargs=1, type=int, default=25,
                    help="max allowed distance")
parser.add_argument("-t", nargs=1, type=str, default="../tmp/",
                    help="temporary folder path")
parser.add_argument("-m", nargs=1, type=bool, default=False,
                    help="perform an additional minimization round prior to estimating energies")
parser.add_argument("-e", nargs=1, type=str, default="total",
                    help="energy term to compute (Coul/LJ-SR/LR, total or none)")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_file = path.abspath(args.i)

if type(args.o) is list:
    args.o = args.o[0]
output_file = path.abspath(args.o)

if type(args.distance) is list:
    args.d = args.d[0]
max_distance = args.distance

chdir(path.dirname(path.realpath(__file__)))

# Create temporary folders
pdb_dir = path.abspath(args.t + '/pdb/')

if not path.exists(pdb_dir):
    makedirs(pdb_dir)

param_template_path = path.abspath("../res/")

# Writing output file header
col_names = ['pdb_id', 'species',
             'mhc_type',
             'tcr_v_allele', 'tcr_region', 'tcr_region_seq',
             'aa_tcr', 'aa_mhc', 'len_tcr', 'len_mhc',
             'pos_tcr', 'pos_mhc',
             'distance', "energy"]

with open(output_file, 'w') as f:
    f.write('\t'.join(col_names) + '\n')

pdb_list = PDBList()
pdb_parser = PDBParser()

table = pd.read_table(input_file)

bypdb = table.groupby("pdb_id")

i = 0

for pdb_id, pdb_group in bypdb:

    # Load PDB file
    pdb_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=pdb_dir)

    print("[", i, "/", table.shape[0], "]")
    print(pdb_id, "- preparing for computation")

    # Load model from original pdb file
    model = pdb_parser.get_structure(pdb_id, pdb_file)[0]

    # Store annotation for entire complex
    pdb_annot = pdb_group.iloc[0]

    # Get and check antigen residues
    mhc_a_chain = model[pdb_annot['chain_mhc_a']]
    mhc_a_range = range(len(model[pdb_annot['chain_mhc_a']]))
    mhc_a_residues = get_residues(mhc_a_chain, mhc_a_range)

    mhc_b_chain = model[pdb_annot['chain_mhc_b']]
    mhc_b_range = range(len(model[pdb_annot['chain_mhc_b']]))
    mhc_b_residues = get_residues(mhc_b_chain, mhc_b_range)

    # Iterate by TCR chain/region (CDR1,2,3 + FR1,2,3)
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

        print(pdb_id, "- computing energies for", tcr_annot['tcr_v_allele'], ':', tcr_region_id[1])

        # Compute distances and add them to results
        distances = calc_distances_mhc(tcr_chain.get_id(), mhc_a_chain.get_id(),
                                   tcr_annot['tcr_v_allele'], tcr_annot['tcr_region'],
                                   tcr_region_residues, mhc_a_residues, tcr_region_range, mhc_a_range)
        for row in distances:
            if row["distance"] <= max_distance:
                row.update(tcr_annot.to_dict())
                row['energy'] = "NA"
                results_by_pdb.append(row)


        distances = calc_distances_mhc(tcr_chain.get_id(), mhc_b_chain.get_id(),
                                   tcr_annot['tcr_v_allele'], tcr_annot['tcr_region'],
                                   tcr_region_residues, mhc_b_residues, tcr_region_range, mhc_b_range)
        for row in distances:
            if row["distance"] <= max_distance:
                row.update(tcr_annot.to_dict())
                row['energy'] = "NA"
                results_by_pdb.append(row)


        i += 1

    # Write selected columns and delete gmx/ content
    res = pd.DataFrame(results_by_pdb)[col_names]
    res.to_csv(output_file, sep='\t', header=False, index=False, mode='a')
    print("Done")

print("Finished processing", table.shape[0], "entries.")
