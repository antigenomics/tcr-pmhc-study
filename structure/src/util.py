import itertools
import os
import sys


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
    return [{'tcr_v_allele': tcr_v_allele,
             'tcr_region': tcr_region,
             'aa_tcr': get_aa_code(aa_tcr),
             'aa_ag': get_aa_code(aa_ag),
             'len_tcr': len(tcr_range),
             'len_ag': len(ag_range),
             'pos_tcr': tcr_range[i],
             'pos_ag': ag_range[j],
             'dist': calc_distance(aa_tcr, aa_ag)}
            for i, aa_tcr in enumerate(tcr_residues)
            for j, aa_ag in enumerate(ag_residues)]
