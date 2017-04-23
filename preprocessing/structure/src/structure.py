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
parser.add_argument("-m", nargs=1, type=str, default="../../result/structure_mock.txt",
                    help="Mock complexes output table path.")
parser.add_argument("-t", nargs=1, type=str, default="../tmp/",
                    help="Temporary folder path")
parser.add_argument("--pdbfixer", action="store_true",
                    help="Fix PDB files with PDBFixer (time-consuming!).")
parser.add_argument("--mock", action="store_true",
                    help="Generate mock complexes (negative controls).")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_file = path.abspath(args.i)
if type(args.o) is list:
    args.o = args.o[0]
output_file = path.abspath(args.o)
if type(args.m) is list:
    args.m = args.m[0]
output_file_mock = path.abspath(args.m)
use_pdbfixer = args.pdbfixer
generate_mock = args.mock

chdir(path.dirname(path.realpath(__file__)))

# Create temporary folders
pdb_dir = path.abspath(args.t + '/pdb/')

if not path.exists(pdb_dir):
    makedirs(pdb_dir)

# Writing output file header, we'll drop some original annotation
# columns such as chain identifiers, etc
col_names = ['pdb_id', 'species',
             'mhc_type', 'mhc_a_allele', 'mhc_b_allele',
             'antigen_seq', 'tcr_gene',
             'tcr_v_allele', 'tcr_region', 'tcr_region_seq',
             'aa_tcr', 'aa_antigen', 'len_tcr', 'len_antigen',
             'pos_tcr', 'pos_antigen',
             'distance', 'distance_CA']

with open(output_file, 'w') as f:
    f.write('\t'.join(col_names) + '\n')

# Main loop

pdb_list = PDBList()
pdb_parser = PDBParser()

table = pd.read_table(input_file)

bypdb = table.groupby("pdb_id")

i = 0

df_residues = []

for pdb_id, pdb_group in bypdb:
    # Load PDB file
    pdb_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=pdb_dir)

    print("[", time.strftime("%c"), i, "/", table.shape[0], "]")
    print(pdb_id, "- preparing for computation")

    # Load model from original pdb file
    # In case of PDBFixer option it will be used to find
    # region ranges by overlapping Calpha atom positions
    model_original = pdb_parser.get_structure(pdb_id, pdb_file)[0]

    if use_pdbfixer:
        print(pdb_id, "-- fixing PDB")
        pdb_file_fixed = fix_pdb(pdb_id, pdb_file, pdb_group)
        # fixed model
        model = pdb_parser.get_structure(pdb_id, pdb_file_fixed)[0]
    else:
        model = model_original

    # Store annotation for entire complex
    pdb_annot = pdb_group.iloc[0]

    # Get and check antigen residues
    antigen_chain_id = pdb_annot['chain_antigen']
    antigen_chain = model[antigen_chain_id]
    antigen_seq = pdb_annot['antigen_seq']
    antigen_range = range(len(antigen_seq))

    if use_pdbfixer:
        antigen_residues = get_residues_pdbfixer(antigen_chain, antigen_range,
                                                 model_original[antigen_chain_id])
    else:
        antigen_residues = get_residues(antigen_chain, antigen_range)

    antigen_seq_obs = get_seq(antigen_residues)

    if antigen_seq != antigen_seq_obs:
        warning(pdb_id, "Antigen sequence mismatch (expected observed):", antigen_seq,
                antigen_seq_obs, ". Replacing with one from PDB.")
        pdb_annot['antigen_seq'] = antigen_seq_obs

    mhc_type = pdb_annot['mhc_type']

    # Iterate by TCR chain/region (CDR1,2,3 )
    byregion = pdb_group.groupby(['chain_tcr', 'tcr_region'])
    results_by_pdb = []
    for tcr_region_id, tcr_region_group in byregion:
        # Get and check tcr region residues
        tcr_annot = tcr_region_group.iloc[0]
        tcr_region_name = tcr_annot['tcr_region']

        if tcr_region_name.startswith("FR"):
            continue

        tcr_v_name = tcr_annot['tcr_v_allele']
        tcr_gene = tcr_v_name[0:3]
        tcr_chain_id = tcr_annot['chain_tcr']
        tcr_chain = model[tcr_chain_id]
        tcr_region_seq = tcr_annot['tcr_region_seq']
        tcr_region_range = range(tcr_annot['tcr_region_start'],
                                 tcr_annot['tcr_region_end'])

        if use_pdbfixer:
            tcr_region_residues = get_residues_pdbfixer(tcr_chain, tcr_region_range,
                                                        model_original[tcr_chain_id])
        else:
            tcr_region_residues = get_residues(tcr_chain, tcr_region_range)

        tcr_region_seq_obs = get_seq(tcr_region_residues)

        if tcr_region_seq != tcr_region_seq_obs:
            warning(pdb_id, "TCR:", tcr_region_id, " sequence mismatch (expected observed): ",
                    tcr_region_seq, tcr_region_seq_obs,
                    ". Replacing with one from PDB.")
            tcr_annot['tcr_region_seq'] = tcr_region_seq_obs

        # we'll need this for mock complexes
        if generate_mock and tcr_region_name == "CDR3":
            df_residues.append({'pdb_id': pdb_id,
                                'mhc_type': mhc_type,
                                'ag_seq': antigen_seq_obs,
                                'tcr_seq': tcr_region_seq_obs,
                                'ag_residues': antigen_residues,
                                'tcr_gene': tcr_gene,
                                'tcr_region': tcr_region_name,
                                'tcr_residues': tcr_region_residues})

        # Compute distances and add them to results

        print(pdb_id, "- computing distances for", tcr_v_name,
              ':', tcr_region_id[1])

        distances = calc_distances(tcr_region_residues, antigen_residues)

        # Append annotation
        for row in distances:
            row.update(tcr_annot.to_dict())
            row.update({'tcr_gene': tcr_gene,
                        'tcr_v_allele': tcr_v_name,
                        'tcr_region': tcr_region_name})

        results_by_pdb.extend(distances)

        i += 1

    # Write selected columns
    res = pd.DataFrame(results_by_pdb)[col_names]
    res.to_csv(output_file, sep='\t', header=False, index=False, mode='a')
    print("Done")

print("Finished processing", table.shape[0], "entries.")

if generate_mock:
    print("Generating mock complexes.")

    col_names_mock = ['mhc_type', 'tcr_gene', 'tcr_region',
                      'pdb_id_a', 'pdb_id_t',
                      'ag_seq_a', 'ag_seq_t',
                      'tcr_seq_a', 'tcr_seq_t',
                      'aa_tcr', 'aa_antigen',
                      'len_tcr', 'len_antigen',
                      'pos_tcr', 'pos_antigen',
                      'distance', 'distance_CA']

    with open(output_file_mock, 'w') as f:
        f.write('\t'.join(col_names_mock) + '\n')

    results_mock = []

    df_residues = pd.DataFrame(df_residues)

    # split in two tables

    df_res_tcr = df_residues[['pdb_id', 'ag_seq', 'tcr_seq',
                             'tcr_gene', 'tcr_region', 'tcr_residues']]

    df_res_tcr.columns = ['pdb_id_t', 'ag_seq_t', 'tcr_seq_t',
                          'tcr_gene', 'tcr_region', 'tcr_residues_t']

    df_res_ag = df_residues[['pdb_id', 'mhc_type', 'ag_seq', 'tcr_seq',
                            'ag_residues', 'tcr_gene', 'tcr_region',
                            'tcr_residues']]

    df_res_ag.columns = ['pdb_id_a', 'mhc_type', 'ag_seq_a', 'tcr_seq_a',
                         'ag_residues', 'tcr_gene', 'tcr_region',
                         'tcr_residues_a']

    # merge (+ cartesian product)

    df_res = pd.merge(df_res_tcr, df_res_ag,
                      on = ['tcr_gene', 'tcr_region'])

    # filter

    df_res = df_res[(df_res.pdb_id_t != df_res.pdb_id_a) &
                    (df_res.ag_seq_t != df_res.ag_seq_a) &
                    (df_res.tcr_seq_t != df_res.tcr_seq_a)]

    # compute distances for all pairs

    for index, row in df_res.iterrows():
        if index % 100 == 0:
            print("[", time.strftime("%c"), index, "/", df_res.shape[0], "]",
                  "Processing mock entry", row['pdb_id_a'], row['pdb_id_t'],
                  row['ag_seq_a'], row['tcr_seq_t'])

        # transfer CDRs to expected position in original pMHC
        tcr_residues_new = superimpose(row['tcr_residues_a'],
                                       row['tcr_residues_t'],
                                       row['tcr_region'])

        distances = calc_distances(tcr_residues_new, row['ag_residues'])

        for row2 in distances:
            row2.update(row.to_dict())

        res = pd.DataFrame(distances)[col_names_mock]
        res.to_csv(output_file_mock, sep='\t', header=False, index=False, mode='a')

    print("Finished processing", df_res.shape[0], "mock entries.")
