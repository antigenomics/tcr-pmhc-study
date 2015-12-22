import itertools
import os
import subprocess
import sys

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

def get_required_chains(pdb_group):
    return set(itertools.chain(*[pdb_group[c].tolist() for c in
                                 ['chain_tcr', 'chain_antigen', 'chain_mhc_a', 'chain_mhc_b']]))


'''
def select_model(pdb_id, pdb_file, pdb_group):
    required_chains = get_required_chains(pdb_group)

    print(PDBParser().get_structure(pdb_id, pdb_file))

    for i, model in enumerate(PDBParser().get_structure(pdb_id, pdb_file)):
        print(i, set([chain.get_id() for chain in model]), required_chains, sep="\t")
        if set([chain.get_id() for chain in model]).issuperset(required_chains):
            return i

    return -1
'''


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
    # fixer.addMissingHydrogens(7.0)
    # fixer.addSolvent(fixer.topology.getUnitCellDimensions())

    # KeepIds flag is critical here, otherwise we loose all information binding
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'), keepIds=True)


_trans = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
          'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
          'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
          'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def get_aa_code(res):
    return _trans[res.get_resname()]


def get_seq(residues):
    return "".join(get_aa_code(res) for res in residues)


def get_residues(residues, range):
    return [x for i, x in enumerate(residues) if i in range]


def calc_distance(aa1, aa2):
    return min([abs(atom1 - atom2) for atom1 in aa1 for atom2 in aa2])


def calc_distances(tcr_v_allele, tcr_region, tcr_residues, ag_residues, tcr_range, ag_range):
    # Here we also store the 'resseq' codes, the only way to traceback to
    # our residues when working with structures fixed by PDBFixer
    return [{'tcr_v_allele': tcr_v_allele,
             'tcr_region': tcr_region,
             'resseq_tcr': aa_tcr.get_id()[1],
             'resseq_ag': aa_ag.get_id()[1],
             'aa_tcr': get_aa_code(aa_tcr),
             'aa_ag': get_aa_code(aa_ag),
             'len_tcr': len(tcr_range),
             'len_ag': len(ag_range),
             'pos_tcr': tcr_range[i],
             'pos_ag': ag_range[j],
             'dist': calc_distance(aa_tcr, aa_ag)}
            for i, aa_tcr in enumerate(tcr_residues)
            for j, aa_ag in enumerate(ag_residues)]


# GROMACS methods

def prepare_gmx(pdb_file, gmx_prefix):
    pdb_file = os.path.abspath(pdb_file)
    gmx_prefix = os.path.abspath(gmx_prefix)
    param_path = os.path.abspath("../res/minim.mdp")

    if os.name == 'nt':  # should run via cygwin on windows
        pdb_file = pdb_file.replace('\\', '\\\\\\\\')
        gmx_prefix = gmx_prefix.replace('\\', '\\\\\\\\')
        param_path = param_path.replace('\\', '\\\\\\\\')

    print(param_path)

    subprocess.check_call(
            'gmx pdb2gmx -f "{}" -o "{}.gro" -p "{}.top" -n "{}" -water spce -ff oplsaa'.
                format(pdb_file, gmx_prefix, gmx_prefix, gmx_prefix),
            shell=True
    )
    subprocess.check_call(
            'gmx grompp -maxwarn 10 -f "{}" -c "{}.gro" -p "{}.top" -n "{}.ndx" -o "{}.tpr"'.
                format(param_path, gmx_prefix, gmx_prefix, gmx_prefix, gmx_prefix),
            shell=True
    )
    subprocess.check_call(
            'gmx editconf -f "{}.gro" -o "{}.gro" -c -d 1.0 -bt cubic'.
                format(gmx_prefix, gmx_prefix),
            shell=True
    )
    subprocess.check_call(
            'gmx mdrun -v -nt 8 -s "{}.tpr" -e "{}.edr"'.
                format(gmx_prefix, gmx_prefix, gmx_prefix),
            shell=True
    )
    # mdrun -v -deffnm em$NAME -nb cpu