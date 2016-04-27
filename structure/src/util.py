from __future__ import print_function

import itertools
import os
import sys
import re

from collections import OrderedDict
from os.path import dirname
from shutil import copyfile
from subprocess import check_call
from Bio.PDB import PDBParser
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
from Bio.PDB.Vector import rotaxis


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


#def calc_distance(aa1, aa2):
#    return min([abs(atom1 - atom2) for atom1 in aa1 for atom2 in aa2])

def get_expected_cb_pos(aa):
	n = aa['N'].get_vector() 
	c = aa['C'].get_vector() 
	ca = aa['CA'].get_vector()
	# center at origin
	n = n - ca
	c = c - ca
	# find rotation matrix that rotates n -120 degrees along the ca-c vector
	rot = rotaxis(-3.14*120.0/180.0, c)
	# apply rotation to ca-n vector
	cb_at_origin = n.left_multiply(rot)
	# put on top of ca atom
	return (cb_at_origin + ca)

def get_cbeta_pos(aa):
	if aa.get_resname() == 'GLY':
		return get_expected_cb_pos(aa)
	else:
		ca = aa['CA'].get_vector()
		try:
			cb = aa['CB'].get_vector()
		except KeyError:
			cb = get_expected_cb_pos(aa)
		#print ((ca + (((cb - ca).normalized()) ** 2.4)))
		return (ca + (((cb - ca).normalized()) ** 2.4))

def calc_distance(aa1, aa2):
    #return min([abs(atom1 - atom2) for atom1 in aa1 for atom2 in aa2])
    return (get_cbeta_pos(aa1) - get_cbeta_pos(aa2)).norm()


def calc_distances(tcr_chain, antigen_chain, tcr_v_allele, tcr_region, tcr_residues, ag_residues, tcr_range, ag_range):
    # indexes are required to access gromacs data
    return [{'tcr_v_allele': tcr_v_allele,
             'tcr_region': tcr_region,
             'idx_tcr': tcr_chain + '_' + str(aa_tcr.get_id()[1]),
             'idx_antigen': antigen_chain + '_' + str(aa_ag.get_id()[1]),
             'aa_tcr': get_aa_code(aa_tcr),
             'aa_antigen': get_aa_code(aa_ag),
             'len_tcr': len(tcr_range),
             'len_antigen': len(ag_range),
             'pos_tcr': i,
             'pos_antigen': j,
             'distance': calc_distance(aa_tcr, aa_ag)}
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

    # KeepIds flag is critical here, otherwise we loose all information binding
    pdb_file = dirname(pdb_file) + '/' + pdb_id + '.pdb'
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'), keepIds=True)

    return pdb_file


# GROMACS

os.environ['GMX_MAXBACKUP'] = '0'  # suppress shitload of gromacs backups


def create_index(pdb_id, pdb_sub_id, keep_idxs):
    # TODO: check AAs

    index = OrderedDict()

    with open(pdb_id + '.pdb') as g_file:
        for _ in range(4):
            next(g_file)

        for line in g_file:
            if line[0:3] == 'TER':
                continue
            if line[0:6] == 'ENDMDL':
                break

            res_id = int(line[22:26].strip())
            atom_id = int(line[4:11].strip())
            chain_id = line[20:22].strip()

            idx = chain_id + "_" + str(res_id)

            if idx in keep_idxs:
                index.setdefault(idx, []).append(atom_id)

    res_ids = []

    with open(pdb_sub_id + '.ndx', 'w') as i_file:
        with open(pdb_sub_id + '.dat', 'w') as d_file:
            d_file.write(str(len(index)) + '\n')
            for res_id, atom_ids in index.items():
                res_ids.append(res_id)
                i_file.write('[ ' + res_id + ' ]\n')
                d_file.write(res_id + '\n')
                for atom_id in atom_ids:
                    i_file.write(str(atom_id) + '\n')

    missing_res = [x for x in keep_idxs if x not in res_ids]
    if missing_res:
        warning("The following residues are missing in processed structure: ",
                "; ".join(missing_res))

    return res_ids


def prepare_gmx(pdb_id, pdb_file, gmx_dir, param_template_path):
    param_template_path += "/ions.mdp"

    check_call('gmx pdb2gmx -f {0} -o {1}.pdb -p {1}.top -water spce -ff oplsaa >> gmx.log 2>&1; '
               'gmx editconf -f {1}.pdb -o {1}.pdb -c -d 1.2 -bt cubic >> gmx.log 2>&1; '
               'gmx solvate -cp {1}.pdb -cs spc216.gro -o {1}.pdb -p {1}.top >> gmx.log 2>&1; '.
               # 'gmx grompp -f {2} -c {1}.pdb -p {1}.top -o {1}.ions.tpr >> gmx.log 2>&1; '
               # 'gmx genion -nice 0 -neutral -s {1}.ions.tpr -o {1}.pdb -p {1}.top '
               # '-pname NA -nname CL -conc 0.154 -rmin 0.4 >> gmx.log 2>&1'.
               format(pdb_file, pdb_id, param_template_path), shell=True
               )


def run_minimization(pdb_id, param_template_path):
    param_template_path += '/minim.mdp'

    check_call('gmx grompp -c {0}.pdb -p {0}.top -f {1} -o {0}.m.tpr >> gmx.log 2>&1; '
               'gmx mdrun -s {0}.m.tpr -c {0}.m.pdb -nb cpu >> gmx.log 2>&1; '
               .format(pdb_id, param_template_path), shell=True
               )


def run_single_point(pdb_id, pdb_sub_id, param_template_path, minimized, keep_idxs):
    copyfile(param_template_path + '/sp.mdp', pdb_sub_id + '.sp.mdp')

    with open(pdb_sub_id + '.sp.mdp', 'a') as f:
        f.write('\nenergygrps\t= ' + ' '.join(keep_idxs))

    suffix = '.m' if minimized else ''

    # Writing the following lines take several days. Funny quotes written by GROMACS to stderr and
    # mailing list replies containing not a single word except a (broken) link to extremely concise
    # documentation were really helpful during this process.

    # Single-point energies are calculated with -rerun
    # ndx/groups files are used to mark residues of interest and calculate interaction energies
    # NOTE: without -nb CPU energies are not stored.

    check_call('gmx grompp -c {1}{2}.pdb -p {1}.top -n {0}.ndx -f {0}.sp.mdp -o {0}.tpr >> gmx.log 2>&1; '
               'gmx mdrun -s {0}.tpr -rerun {1}{2}.pdb -e {0}.edr -nb cpu >> gmx.log 2>&1; '
               'gmx enemat -groups {0}.dat -nlevels 10000 -f {0}.edr -emat .{0}.xpm >> gmx.log 2>&1'
               .format(pdb_sub_id, pdb_id, suffix), shell=True
               )


def read_enemat(pdb_sub_id, ene_term, idxs):
    with open(ene_term + '.' + pdb_sub_id + '.xpm') as x_file:
        line = ""
        while "static char *gromacs_xpm[] = {" not in line:
            line = next(x_file)

        next(x_file)

        levels = {}
        while "x-axis" not in line:
            line = next(x_file)
            tokens = re.split('[ "]', line)
            levels[tokens[1]] = float(tokens[8])

        next(x_file)

        word_size = len(list(levels.keys())[0])

        results = {}
        i = 0
        for line in x_file:
            for j in range(len(idxs)):
                results[(idxs[-i - 1], idxs[j])] = levels[line[(1 + word_size * j):(1 + word_size * (j + 1))]]
            i += 1

    return results
