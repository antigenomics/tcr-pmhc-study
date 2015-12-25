from __future__ import print_function

import itertools
import os
import sys
from collections import OrderedDict
from os.path import dirname
from shutil import copyfile
from subprocess import check_call

from Bio.PDB import PDBParser
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile


def clear_folder(d):
    for the_file in os.listdir(d):
        file_path = os.path.join(d, the_file)
        if os.path.isfile(file_path):
            os.unlink(file_path)


def warning(*objs):
    print("WARNING: ", *objs, file=sys.stderr)


# (Bio)PDB methods

_trans = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
          'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
          'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
          'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
          # modified residues go below
          'PFF': 'F', 'F2F': 'F'}


def get_aa_code(res):
    return _trans[res.get_resname()]


def get_seq(residues):
    return "".join(get_aa_code(res) for res in residues)


def get_residues(residues, range):
    return [x for i, x in enumerate(residues) if i in range]


def calc_distance(aa1, aa2):
    return min([abs(atom1 - atom2) for atom1 in aa1 for atom2 in aa2])


def calc_distances(tcr_chain, antigen_chain, tcr_v_allele, tcr_region, tcr_residues, ag_residues, tcr_range, ag_range):
    # indexes are required to access gromacs data
    return [{'tcr_v_allele': tcr_v_allele,
             'tcr_region': tcr_region,
             'idx_tcr': tcr_chain + '_' + str(aa_tcr.get_id()[1]),
             'idx_ag': antigen_chain + '_' + str(aa_ag.get_id()[1]),
             'aa_tcr': get_aa_code(aa_tcr),
             'aa_ag': get_aa_code(aa_ag),
             'len_tcr': len(tcr_range),
             'len_ag': len(ag_range),
             'pos_tcr': tcr_range[i],
             'pos_ag': ag_range[j],
             'dist': calc_distance(aa_tcr, aa_ag)}
            for i, aa_tcr in enumerate(tcr_residues)
            for j, aa_ag in enumerate(ag_residues)]


# PDBfixer


def get_required_chains(pdb_group):
    return set(itertools.chain(*[pdb_group[c].tolist() for c in
                                 ['chain_tcr', 'chain_antigen', 'chain_mhc_a', 'chain_mhc_b']]))


def fix_pdb(pdb_id, pdb_file, pdb_group):
    chains_to_retain = get_required_chains(pdb_group)
    chains_to_remove = []

    for chain in PDBParser().get_structure(pdb_id, pdb_file)[0]:
        if chain.get_id() not in chains_to_retain:
            chains_to_remove.append(chain.get_id())

    fixer = PDBFixer(filename=pdb_file)

    fixer.removeChains(chainIds=chains_to_remove)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(True)
    fixer.addMissingHydrogens(7.0)
    # fixer.addSolvent(fixer.topology.getUnitCellDimensions())

    # KeepIds flag is critical here, otherwise we loose all information binding
    pdb_file = dirname(pdb_file) + '/' + pdb_id + '.pdb'
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'), keepIds=True)

    return pdb_file


# GROMACS

os.environ['GMX_MAXBACKUP'] = '0'  # suppress shitload of gromacs backups


def create_index(pdb_id, gmx_dir, chains, keep_idxs):
    chains_iter = iter(chains)
    i = 0
    index = OrderedDict()

    with open(gmx_dir + '/' + pdb_id + '.gro') as g_file:
        for _ in range(2):
            next(g_file)

        prev_res_id = -1
        chain = next(chains_iter)

        for line in g_file:
            if len(line) >= 44:
                line = line.replace('\t', '    ')
                res_id = int(line[0:5].strip())
                atom_id = int(line[15:20].strip())

                if res_id < prev_res_id:
                    i += 1
                    if i >= len(chains):
                        break
                    chain = next(chains_iter)

                idx = chain + "_" + str(res_id)

                if idx in keep_idxs:
                    index.setdefault(idx, []).append(atom_id)

                prev_res_id = res_id

    res_ids = []

    with open(gmx_dir + '/' + pdb_id + '.ndx', 'w') as i_file:
        with open(gmx_dir + '/' + pdb_id + '.dat', 'w') as d_file:
            d_file.write(str(len(index)) + '\n')
            for res_id, atom_ids in index.items():
                res_ids.append(res_id)
                i_file.write('[ ' + res_id + ' ]\n')
                d_file.write(res_id + '\n')
                for atom_id in atom_ids:
                    i_file.write(str(atom_id) + '\n')

    return res_ids


def prepare_gmx(pdb_id, pdb_file, gmx_dir):
    os.chdir(gmx_dir)

    check_call('gmx pdb2gmx -f {0} -o {1}.gro -p {1}.top -water spce -ff oplsaa; '
               'gmx editconf -f {1}.gro -o {1}.gro -c -d 1.0 -bt cubic; '.
               format(pdb_file, pdb_id), shell=True
               )


def run_md(pdb_id, gmx_dir, keep_idxs):
    param_template_path = os.path.abspath("../res/sp.mdp")  # TODO: abs path from script

    os.chdir(gmx_dir)

    copyfile(param_template_path, pdb_id + '.mdp')

    with open(pdb_id + '.mdp', 'a') as f:
        f.write('\nenergygrps\t= ' + ' '.join(keep_idxs))

    # Writing the following lines take several days owing to best-practices appliend in writing GROMACS documentation
    # and the intuitiveness and diversity of GROMACS mailing lists Funny quotes written by GROMACS to stderr and
    # mailing list replies containing solely (broken) documentation links to pages containing 5-10 words were really
    # helpful during this process.

    # Single-point energies are calculated with -rerun
    # ndx/groups files are used to mark residues of interest and calculate interaction energies
    # NOTE: without -nb CPU energies are not stored.

    check_call('gmx grompp -c {}.gro -p {}.top -n {}.ndx -f {}.mdp -o {}.tpr; '
               'gmx mdrun -s {}.tpr -rerun {}.gro -e {}.edr -g {}.log -nb cpu; '
               'gmx enemat -groups {}.dat -nlevels 10000 -f {}.edr -emat .{}.xpm'
               .format(pdb_id), shell=True
               )