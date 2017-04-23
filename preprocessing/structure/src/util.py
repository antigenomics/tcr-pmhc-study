from __future__ import print_function

import pandas as pd
import itertools
import os
import sys
import re
import numpy as np
from numpy.linalg import inv
import math

from collections import OrderedDict, namedtuple
from os.path import dirname
from shutil import copyfile
from subprocess import check_call
from Bio.PDB import PDBParser, Vector
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
from simtk.unit import is_quantity, angstroms, norm
from Bio.PDB.Vector import rotaxis, rotmat


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

_skip = set(["HOH", "EMC", "NDG", "NAG", "GOL", "IOD"])


def get_aa_code(res):
    return _trans[res.get_resname()]


def is_hydrogen(atom):
    return atom.get_name().startswith("H")


def get_seq(residues):
    return "".join(get_aa_code(res) for res in residues)


def get_residues(residues, query_range):
    return [x for i, x in enumerate(residues) if i in query_range]


def unpack_residues(residues):
    return pd.DataFrame(list(itertools.chain(*[[(i, res.get_resname(), j, atom.get_name(), atom.get_vector()) for
             j, atom in enumerate(res)] for i, res in enumerate(residues)])))


def get_residues_pdbfixer(residues, query_range, original_residues):
    # If we match everything is OK
    exact = get_residues(residues, query_range)
    if exact == get_residues(original_residues, query_range):
        return exact

    # Find range in original chain data
    original_residue_list = get_residues(original_residues, query_range)

    first = original_residue_list[0]
    best_dist_first = 1000
    best_first = None
    last = original_residue_list[-1]
    best_dist_last = 1000
    best_last = None

    # Find residues in new (fixed) chain data that have nearest CA coordinates
    for i, x in enumerate(residues):
        dist_first = calc_distance(first, x, 'CA')
        if dist_first < best_dist_first:
            best_dist_first = dist_first
            best_first = (i, x)
        dist_last = calc_distance(last, x, 'CA')
        if dist_last < best_dist_last:
            best_dist_last = dist_last
            best_last = (i, x)

    if best_dist_first > 0.1:
        warning("Failed to find CA match between", first,
                "and residues from fixed PDB, best match is", best_first)
    if best_dist_last > 0.1:
        warning("Failed to find CA match between", last,
                "and residues from fixed PDB, best match is", best_last)

    new_range = range(best_first[0], best_last[0] + 1)
    return get_residues(residues, new_range)


def get_calpha_pos(aa):
    if 'CA' in aa:
        return aa['CA'].get_vector()
    else:
        return Vector(float('nan'), float('nan'), float('nan'))


# Distance & energy


def calc_distance(aa1, aa2, dist_type = 'closest_atom'):
    if dist_type == 'CA':
        return (get_calpha_pos(aa1) - get_calpha_pos(aa2)).norm()
    else: # min dist
        return min([abs(atom1 - atom2) for atom1 in aa1 if not is_hydrogen(atom1)
            for atom2 in aa2 if not is_hydrogen(atom2)])


def calc_distances(tcr_residues, ag_residues):
    len_tcr = len(tcr_residues)
    len_antigen = len(ag_residues)

    return [{'aa_tcr': get_aa_code(aa_tcr),
             'aa_antigen': get_aa_code(aa_ag),
             'len_tcr': len_tcr,
             'len_antigen': len_antigen,
             'pos_tcr': i,
             'pos_antigen': j,
             'distance': calc_distance(aa_tcr, aa_ag),
             'distance_CA': calc_distance(aa_tcr, aa_ag, 'CA')}
            for i, aa_tcr in enumerate(tcr_residues)
            for j, aa_ag in enumerate(ag_residues) if aa_ag.get_resname() not in _skip]


# Rotations and superimpositions

CanonicalTransformation = namedtuple('CanonicalTransformation',
                                     'offset Rx Rz')


def find_canonical_transform(residues):
    start = get_calpha_pos(residues[0])
    end = get_calpha_pos(residues[-1])

    # C terminus -> N terminus axis
    direction = end - start

    # Matrix that rotates 'direction' to align with X axis
    # i.e. rotation along Z axis
    Rz = rotmat(direction, Vector(1, 0, 0))

    res = [(i, aa, (get_calpha_pos(aa) - start).left_multiply(Rz))
        for i, aa in enumerate(residues)]

    # Find 'center of mass' coord in YZ plane
    cm = Vector(0, 0, 0)

    for _, _, vec in res:
        cm = cm + Vector(0, vec[1], vec[2])

    # Matrix that rotates 'center of mass' to align with Z axis in YZ plane
    # i.e. rotation along X axis
    Rx = rotmat(cm, Vector(0, 0, 1))

    Rx[0,:] = np.array([1.0, 0.0, 0.0])
    Rx[:,0] = np.array([1.0, 0.0, 0.0])

    return CanonicalTransformation(offset = start,
                                   Rz = Rz, Rx = Rx)


def move_residue(res, sfun):
    res1 = res.copy()

    for atom in res1:
        atom.set_coord(sfun(atom.get_vector()).get_array())

    return res1


def to_canonical(residues, ct):
    ctfun = lambda x: (x - ct.offset).left_multiply(ct.Rz).left_multiply(ct.Rx)

    return [move_residue(res, ctfun) for res in residues]


def from_canonical(residues, ct):
    Rxi = inv(ct.Rx)
    Rzi = inv(ct.Rz)
    ctifun = lambda x: ct.offset + x.left_multiply(Rxi).left_multiply(Rzi)

    return [move_residue(res, ctifun) for res in residues]


def superimpose(residues_ref, residues_query, tcr_region):
    # Find canonical transformation for query and apply it
    ctq = find_canonical_transform(residues_query)
    residues_query = to_canonical(residues_query, ctq)

    # Find canonical transformation for reference
    ctr = find_canonical_transform(residues_ref)

    # Move canonicalized query atoms to reference
    residues_query = from_canonical(residues_query, ctr)

    return residues_query



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
    fixer.addMissingHydrogens(7.4)

    # KeepIds flag is critical here, otherwise we loose all information binding
    pdb_file = dirname(pdb_file) + '/' + pdb_id + '.pdb'

    f = open(pdb_file, 'w')

    PDBFile.writeHeader(fixer.topology, f)
    PDBFile_writeModel(fixer.topology, fixer.positions, f)
    PDBFile.writeFooter(fixer.topology, f)

    return pdb_file

# Below goes modified OPENMM code
# from https://github.com/pandegroup/openmm/blob/f21bb6a5434eafe1fbf315cf686fd2d198ddbde0/wrappers/python/simtk/openmm/app/pdbfile.py
# Preserving chain IDs but renumberig residues, this solves https://github.com/pandegroup/pdbfixer/issues/153 for us

def PDBFile_writeModel(topology, positions, file=sys.stdout, modelIndex=None):
    """Write out a model to a PDB file.
    Parameters
    ----------
    topology : Topology
        The Topology defining the model to write
    positions : list
        The list of atomic positions to write
    file : file=stdout
        A file to write the model to
    modelIndex : int=None
        If not None, the model will be surrounded by MODEL/ENDMDL records
        with this index
    """

    if len(list(topology.atoms())) != len(positions):
        raise ValueError('The number of positions must match the number of atoms')
    if is_quantity(positions):
        positions = positions.value_in_unit(angstroms)
    if any(math.isnan(norm(pos)) for pos in positions):
        raise ValueError('Particle position is NaN')
    if any(math.isinf(norm(pos)) for pos in positions):
        raise ValueError('Particle position is infinite')
    nonHeterogens = PDBFile._standardResidues[:]
    nonHeterogens.remove('HOH')
    atomIndex = 1
    posIndex = 0
    if modelIndex is not None:
        print("MODEL     %4d" % modelIndex, file=file)
    for (chainIndex, chain) in enumerate(topology.chains()):
        chainName = chain.id

        residues = list(chain.residues())
        for (resIndex, res) in enumerate(residues):
            if len(res.name) > 3:
                resName = res.name[:3]
            else:
                resName = res.name

            resId = "%4d" % ((resIndex+1)%10000)

            if res.name in nonHeterogens:
                recordName = "ATOM  "
            else:
                recordName = "HETATM"
            for atom in res.atoms():
                if atom.element is not None:
                    symbol = atom.element.symbol
                else:
                    symbol = ' '
                if len(atom.name) < 4 and atom.name[:1].isalpha() and len(symbol) < 2:
                    atomName = ' '+atom.name
                elif len(atom.name) > 4:
                    atomName = atom.name[:4]
                else:
                    atomName = atom.name
                coords = positions[posIndex]
                line = "%s%5d %-4s %3s %s%4s    %s%s%s  1.00  0.00          %2s  " % (
                    recordName, atomIndex%100000, atomName, resName, chainName, resId, _format_83(coords[0]),
                    _format_83(coords[1]), _format_83(coords[2]), symbol)
                assert len(line) == 80, 'Fixed width overflow detected'
                print(line, file=file)
                posIndex += 1
                atomIndex += 1
            if resIndex == len(residues)-1:
                print("TER   %5d      %3s %s%4s" % (atomIndex, resName, chainName, resId), file=file)
                atomIndex += 1
    if modelIndex is not None:
        print("ENDMDL", file=file)


def _format_83(f):
    """Format a single float into a string of width 8, with ideally 3 decimal
    places of precision. If the number is a little too large, we can
    gracefully degrade the precision by lopping off some of the decimal
    places. If it's much too large, we throw a ValueError"""
    if -999.999 < f < 9999.999:
        return '%8.3f' % f
    if -9999999 < f < 99999999:
        return ('%8.3f' % f)[:8]
    raise ValueError('coordinate "%s" could not be represented '
                     'in a width-8 field' % f)
