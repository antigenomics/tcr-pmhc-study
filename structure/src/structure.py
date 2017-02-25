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

parser = argparse.ArgumentParser(description="structure.py")
parser.add_argument("-i", nargs=1, type=str, default="../../result/final.annotations.txt",
                    help="Annotation table path.")
parser.add_argument("-o", nargs=1, type=str, default="../../result/structure.txt",
                    help="Output table path.")
parser.add_argument("-t", nargs=1, type=str, default="../tmp/",
                    help="Temporary folder path")
parser.add_argument("--gmx", action="store_true",
                    help="Compute point energies using GROMACS (very time-consuming!).")
parser.add_argument("-m", nargs=1, type=bool, default=False,
                    help="Perform an additional minimization round prior to estimating energies. Only takes effect if --gmx is specified.")
parser.add_argument("-e", nargs=1, type=str, default="total",
                    help="Energy term to compute: Coul/LJ-SR/LR, total or none. Only takes effect if --gmx is specified.")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_file = path.abspath(args.i)

if type(args.o) is list:
    args.o = args.o[0]
output_file = path.abspath(args.o)

minimized = args.m
energy_term_name = args.e
use_gmx = args.gmx

chdir(path.dirname(path.realpath(__file__)))

# Create temporary folders
(pdb_dir, gmx_dir) = path.abspath(args.t + '/pdb/'), path.abspath(args.t + '/gmx/')

if not path.exists(pdb_dir):
    makedirs(pdb_dir)

if use_gmx:
    print("Using GROMACS...")

    if path.exists(gmx_dir):
        rmtree(gmx_dir)

    makedirs(gmx_dir)

param_template_path = path.abspath("../res/")

# Writing output file header
col_names = ['pdb_id', 'species',
             'mhc_type', 'mhc_a_allele', 'mhc_b_allele',
             'antigen_seq',
             'tcr_v_allele', 'tcr_region', 'tcr_region_seq',
             'aa_tcr', 'aa_antigen', 'len_tcr', 'len_antigen',
             'pos_tcr', 'pos_antigen',
             'distance', 'distance_CA', 'distance_CB', 
             'energy']

with open(output_file, 'w') as f:
    f.write('\t'.join(col_names) + '\n')

# Main loop
# Work in GROMACS path, as inevitably there will be lots of mess created in work dir
# and many paths are specified relative to it

chdir(gmx_dir)

pdb_list = PDBList()
pdb_parser = PDBParser()

table = pd.read_table(input_file)

bypdb = table.groupby("pdb_id")

i = 0

for pdb_id, pdb_group in bypdb:
    # Load PDB file
    pdb_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=pdb_dir)

    print("[", time.strftime("%c"), i, "/", table.shape[0], "]")
    print(pdb_id, "- preparing for computation")

    # Load model from original pdb file
    model = pdb_parser.get_structure(pdb_id, pdb_file)[0]
    
    if use_gmx:
        # Fix PDB structure and make GROMACS files, we'll use it later
        print(pdb_id, "-- fixing PDB")
        pdb_file = fix_pdb(pdb_id, pdb_file, pdb_group)

        print(pdb_id, "-- making topology")
        prepare_gmx(pdb_id, pdb_file, gmx_dir, param_template_path)
        if minimized:
            print(pdb_id, "-- running energy minimization")
            run_minimization(pdb_id, param_template_path)

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
        distances = calc_distances(tcr_chain.get_id(), antigen_chain.get_id(),
                                   tcr_annot['tcr_v_allele'], tcr_annot['tcr_region'],
                                   tcr_region_residues, antigen_residues, tcr_region_range, antigen_range)

        if use_gmx:
            # Run the following separately for each region as GROMACS will "not support >64 energy groups on NxN kernels"
            pdb_sub_id = pdb_id + '.' + '.'.join(tcr_region_id)

            # Assign energy groups, create indexes
            keep_idxs = set()

            for row in distances:
                keep_idxs.add(row['idx_tcr'])
                keep_idxs.add(row['idx_antigen'])

            # Create index and store energy group (aka residue id) order, as GROMACS routines will fail otherwise
            idxs = create_index(pdb_id, pdb_sub_id, keep_idxs)

            # Compute energies
            run_single_point(pdb_id, pdb_sub_id, param_template_path, minimized, idxs)

            # Append general annotation, extract and append energies
            energies = read_enemat(pdb_sub_id, energy_term_name, idxs)

        for row in distances:
            row.update(tcr_annot.to_dict())
            if use_gmx:
                row['energy'] = energies.get((row['idx_tcr'], row['idx_antigen']), float('nan'))
            else:
                row['energy'] = 'NA'

        results_by_pdb.extend(distances)

        i += 1

    # Write selected columns and delete gmx/ content
    res = pd.DataFrame(results_by_pdb)[col_names]
    res.to_csv(output_file, sep='\t', header=False, index=False, mode='a')
    if use_gmx:
        clear_folder(gmx_dir)
    print("Done")

print("Finished processing", table.shape[0], " entries.")
