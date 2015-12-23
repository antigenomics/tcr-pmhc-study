import argparse

import pandas as pd
from Bio.PDB import *

from util import *

parser = argparse.ArgumentParser(description='structure.py annotation_table output_table')

args = parser.parse_args()

# df = pd.read_table(args[0])

table = pd.read_table('../../result/final.annotations.txt')  # TODO

bypdb = table.groupby('pdb_id')

# Create temporary folders
(pdb_dir, gmx_dir) = '../tmp/pdb/', '../tmp/gmx/'
for d in [pdb_dir, gmx_dir]:
    if not os.path.exists(d):
        os.makedirs(d)
    else:
        clear_folder(d)

# Main loop
for pdb_id, pdb_group in bypdb:
    pdb_file = '{}pdb{}.ent'.format(pdb_dir, pdb_id)
    gmx_prefix = '{}{}'.format(gmx_dir, pdb_id)

    # Load PDB file
    PDBList().retrieve_pdb_file(pdb_id, pdir=pdb_dir)

    # Load model from fixed pdb file
    model = PDBParser().get_structure(pdb_id, pdb_file)[0]

    # Store annotation for entire complex
    pdb_annot = pdb_group.iloc[0]

    # Get and check antigen residues
    antigen_chain = model[pdb_annot['chain_antigen']]
    antigen_seq = pdb_annot['antigen_seq']
    antigen_range = range(len(antigen_seq))
    antigen_residues = get_residues(antigen_chain, antigen_range)

    assert antigen_seq == get_seq(antigen_residues)

    # Iterate by TCR chain/region (CDR1,2,3)
    byregion = pdb_group.groupby(['chain_tcr', 'tcr_region'])
    results = []
    for tcr_region_id, tcr_region_group in byregion:
        # Get and check tcr region residues
        tcr_annot = tcr_region_group.iloc[0]
        tcr_chain = model[tcr_annot['chain_tcr']]
        tcr_region_seq = tcr_annot['tcr_region_seq']
        tcr_region_range = range(tcr_annot['tcr_region_start'], tcr_annot['tcr_region_end'])
        tcr_region_residues = get_residues(tcr_chain, tcr_region_range)

        assert tcr_region_seq == get_seq(tcr_region_residues)

        # Compute distances and add them to results
        results.extend(calc_distances(tcr_annot['tcr_v_allele'], tcr_annot['tcr_region'],
                                      tcr_region_residues, antigen_residues, tcr_region_range, antigen_range))

    print(results)

