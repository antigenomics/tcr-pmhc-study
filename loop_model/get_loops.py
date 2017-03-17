from __future__ import print_function, division

import warnings
import argparse
import time
import gzip
import sys
from os import listdir, path
from os.path import isfile, join

import multiprocessing

import pandas as pd
import numpy as np
from Bio.PDB import *

warnings.filterwarnings("ignore")

### CLI

parser = argparse.ArgumentParser(description="structure.py")
parser.add_argument("-i", nargs=1, type=str, default="backbone.txt",
                    help="Path to file containing CDR1-3 coordinates.")
parser.add_argument("-o", nargs=1, type=str, default="loops.txt",
                    help="Output table path.")
parser.add_argument("-d", nargs=1, type=str, default="data/",
                    help="Path to data folder containing PDB files.")
parser.add_argument("-l", nargs=2, type=int, default=[13,13],
                    help="Min and max loop length.")

args = parser.parse_args()
backbone_filename = path.abspath(args.i)
output_filename = path.abspath(args.o)
data_path = path.abspath(args.d)
cdr3_lens = range(args.l[0], args.l[1] + 1)

### UTILS

def get_file_handle(file_name):
    if file_name.endswith('.gz'):
        return gzip.open(file_name, 'rt')
    else:
        return open(file_name)


def get_pdb_id(file_name):
    return file_name.split('.')[0][3:]


pdb_parser = PDBParser()


def load_pdb_model(file_name):
    pdb_id = get_pdb_id(file_name)
    with get_file_handle(file_name) as handle:
        return pdb_parser.get_structure(pdb_id, handle)[0]
    

_trans = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
          'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
          'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
          'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def get_aa_code(res):
    return _trans.get(res.get_resname(), "X")


def get_calpha_pos(res):
    return res['CA'].get_vector()


def check_residues(residues):
    return all(get_aa_code(res) != "X" and 'CA' in res for res in residues)


def rotate(residues):
    start = get_calpha_pos(residues[0])
    end = get_calpha_pos(residues[-1])

    # C terminus -> N terminus axis
    direction = end - start

    # Matrix that rotates 'direction' to align with X axis
    r = rotmat(direction, Vector(1, 0, 0))
    
    res = [(i, aa, (get_calpha_pos(aa) - start).left_multiply(r))
        for i, aa in enumerate(residues)]

    # Find 'center of mass' coord in YZ plane
    cm = Vector(0, 0, 0)
    
    for _, _, vec in res:
        cm = cm + Vector(0, vec[1], vec[2])
    
    # Matrix that rotates 'center of mass' to align with Z axis in YZ plane
    r = rotmat(cm, Vector(0, 0, 1))

    res = [(i, aa, vec[0], Vector(0, vec[1], vec[2]).left_multiply(r))
        for i, aa, vec in res]

    return [(i, len(res), aa, Vector(x, vecYZ[1], vecYZ[2])) 
        for i, aa, x, vecYZ in res]


def calc_coords(pdb_id_kmer, chain_id_kmer, start_kmer, residues):
    return [{'pdb_id_kmer': pdb_id_kmer,
             'chain_id_kmer': chain_id_kmer,
             'start_kmer': start_kmer,
             'aa_kmer': get_aa_code(aa),
             'len_tcr': l,
             'pos_tcr': i,
             'x': vec[0],
             'y': vec[1],
             'z': vec[2]
             }
            for i, l, aa, vec in rotate(residues)]


def compute_rmsd(df, kmer, rmsd_threshold = 1.5):
    result = pd.merge(df, kmer, on=['len_tcr', 'pos_tcr'], suffixes=['', '_kmer'])
    result['delta'] = (result['x'] - result['x_kmer'])**2 + (result['y'] - result['y_kmer'])**2 + (result['z'] - result['z_kmer'])**2
    rmsd = np.sqrt(result.groupby(['pdb_id', 'tcr_v_allele', 'tcr_region', 'pdb_id_kmer', 'chain_id_kmer', 'start_kmer'])['delta'].mean())
    rmsd = rmsd[rmsd <= rmsd_threshold]

    return pd.merge(result, rmsd.to_frame('rmsd').reset_index())

### MAIN

#### Load and filter backbone data

df = pd.read_table(backbone_filename)
df = df.loc[df['len_tcr'].isin(cdr3_lens)][['pdb_id', 'tcr_v_allele', 'tcr_region', 'aa_tcr', 'len_tcr', 'pos_tcr', 'x', 'y', 'z']]

print("[", time.strftime("%c"), "] Loaded backbone data")

#### List PDB files

pdb_files = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.startswith('pdb')]

num_tasks = len(pdb_files)

print("[", time.strftime("%c"), "] Running for", num_tasks, "PDB files")

#### Main loop

def compute_rmsd_pdb(pdb_file):
    pdb_id = get_pdb_id(pdb_file)
    df_rmsd = pd.DataFrame()
    for chain in load_pdb_model(join(data_path, pdb_file)):
        residues = [r for r in chain.get_residues()]
        for cdr3_len in cdr3_lens:
            df_kmer = []
            for i in range(0, len(residues) - cdr3_len):
                kmer = residues[i:(i+cdr3_len)]
                if check_residues(kmer):
                    coords = calc_coords(pdb_id, chain.get_id(), i, kmer)
                    df_kmer.extend(coords)

            if df_kmer:
                df_kmer = pd.DataFrame(df_kmer)
                df_rmsd = df_rmsd.append(compute_rmsd(df, df_kmer))

    progress_string = '\r[ ' + time.strftime("%c") + ' ] Finished processing ' + pdb_id
    sys.stderr.write(progress_string)

    return df_rmsd

rmsd_filtered_all = multiprocessing.Pool().map(compute_rmsd_pdb, pdb_files)

print("\n[", time.strftime("%c"), "] Concatinating results and writing")

pd.concat(rmsd_filtered_all).to_csv(output_filename, sep='\t', index=False)
