import argparse
import warnings
from shutil import rmtree

import pandas as pd
from Bio.PDB import *

from util import *

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="structure.py annotation_table output_table")

args = parser.parse_args()

# df = pd.read_table(args[0])

table = pd.read_table("../../result/final.annotations.txt")  # TODO

bypdb = table.groupby("pdb_id")

# Create temporary folders
(pdb_dir, gmx_dir) = os.path.abspath('../tmp/pdb/'), \
                     os.path.abspath('../tmp/gmx/')

if os.path.exists(gmx_dir):
    rmtree(gmx_dir)

if not os.path.exists(pdb_dir):
    os.makedirs(pdb_dir)

os.makedirs(gmx_dir)

param_template_path = os.path.abspath("../res/sp.mdp")

os.chdir(gmx_dir)

pdb_list = PDBList()
pdb_parser = PDBParser()

results = pd.DataFrame()

# Main loop
for pdb_id, pdb_group in itertools.islice(bypdb, 2, 4):  # TODO
    # pdb_file = "{}pdb{}.ent".format(pdb_dir, pdb_id)

    # Load PDB file
    pdb_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=pdb_dir)

    # Load model from original pdb file
    model = pdb_parser.get_structure(pdb_id, pdb_file)[0]

    # Fix PDB structure and make GROMACS files, we'll use it later
    pdb_file = fix_pdb(pdb_id, pdb_file, pdb_group)
    prepare_gmx(pdb_id, pdb_file, gmx_dir)

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

        # Compute distances and add them to results
        distances = calc_distances(tcr_chain.get_id(), antigen_chain.get_id(),
                                   tcr_annot['tcr_v_allele'], tcr_annot['tcr_region'],
                                   tcr_region_residues, antigen_residues, tcr_region_range, antigen_range)

        # Run the following separately for each region as GROMACS will not support >64 energy groups on NxN kernels
        pdb_sub_id = pdb_id + "." + ".".join(tcr_region_id)

        # Assign energy groups, create indexes
        keep_idxs = set()

        for row in distances:
            keep_idxs.add(row['idx_tcr'])
            keep_idxs.add(row['idx_ag'])

        # won't work if the order is not the same
        idxs = create_index(pdb_id, pdb_sub_id, [chain.get_id() for chain in model], keep_idxs)

        # Run molecular dynamics
        run_single_point(pdb_id, pdb_sub_id, param_template_path, idxs)

        # Append annotation, extract and append energies
        energies = read_enemat(pdb_sub_id, idxs)

        for row in distances:
            row['energy'] = energies.get((row['idx_tcr'], row['idx_ag']), float('nan'))
            row.update(tcr_annot.to_dict())

        results_by_pdb.extend(distances)

    results = results.append(pd.DataFrame(results_by_pdb))

results.to_csv("structure.txt", sep='\t', index=False)
