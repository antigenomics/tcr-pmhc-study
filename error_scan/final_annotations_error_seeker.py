import os
import glob
import pandas as pd
import copy
import Bio.PDB
import time
from antigenomics_util import *
import argparse

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will try to find errors in annotations of complexes TCR+MHC+antigen.")
parser.add_argument("-i", nargs=1, type=str, default="../result/final.annotations.txt",
                    help="Annotation table path.")
parser.add_argument("-o", nargs=1, type=str, default="annotation_errors/",
                    help="Path to output files.")
parser.add_argument("-o2", nargs=1, type=str, default="../result/final.annotations.patch.txt",
                    help="Output path for annotation table patch.")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_file = os.path.abspath(args.i)

if type(args.o) is list:
    args.o = args.o[0]
outp = os.path.abspath(args.o)

if type(args.o2) is list:
    args.o2 = args.o2[0]
output_file = os.path.abspath(args.o2)


tto = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def crdir(pathto):
    if os.path.exists(pathto):
        return 0
    else:
        os.mkdir(pathto)

crdir(outp)
crdir(os.path.join(outp, 'pdb'))

inpdata = pd.read_csv(input_file, sep='\t', dtype=str)
ok_pdbs = pd.unique(inpdata['pdb_id'])


def get_pdbs(inplist, outpath):
    pdbl = Bio.PDB.PDBList()
    iti = 0
    print(len(inplist), iti)
    for i in inplist:
        iti += 1
        if not os.path.isfile(os.path.join(outpath,'pdb'+i+'.ent')):
            print(len(inplist), iti)
            pdbl.retrieve_pdb_file(i, pdir=os.path.join(outpath))

bad_pdbs = {'chain_problem':[], 'seqpos_problem':[], 'bad_structure_problem':[], 'loss_of_data':[], 'possible_problems':[]}
possible_problem = {'chain_problem':[], 'seqpos_problem':[], 'bad_structure_problem':[], 'loss_of_data':[]}

def intersect(a, b):
    return(sorted(list(set(sorted(a))&set(sorted(b)))))

def is_in_list(a, b):
    if sorted(intersect(a, b)) == sorted(a):
        return True
    else:
        return False

def refine(a):
    b = a
    to_pop = set()
    if len(a) == 1:
        return b
    for i in range(1, len(a)+1):
        for j in range(i+1, len(a)+1):
            if is_in_list(a[str(i)], a[str(j)]):
                to_pop.add(str(i))
            elif is_in_list(a[str(j)], a[str(i)]):
                to_pop.add(str(j))
    for i in to_pop:
       b.pop(i)
    return b

pdb_models = {} #{pdb:[biomolecules, pdb_chains_interch]}
for i in ok_pdbs:
    if not os.path.isfile(os.path.join(outp, 'pdb', 'pdb' + i + '.ent')):
        get_pdbs(ok_pdbs, os.path.join(outp, 'pdb'))
    with open(os.path.join("annotation_errors/pdb/pdb"+i+".ent"), 'r') as inp:
        yep = 0
        pdbinfo = {} # for biomolecules
        pdb_chains_interch = []
        for line in inp:
            if line.startswith('COMPND'):
                if line.strip().split()[2] == 'CHAIN:':
                    pdb_chains_interch.append(list(map(lambda x: x.replace(',','').replace(';',''), line.strip().split()[3:])))
            if line.startswith("REMARK 350"):
                #work with biomolecule. Add all of them to pdbinfo.
                yep = 1
                for line in inp:
                    if line.startswith('REMARK 350 BIOMOLECULE:'):
                        biomol = line.strip().split()[3]
                        for line in inp:
                            if line.startswith('REMARK 350 APPLY THE FOLLOWING TO CHAINS:'):
                                pdbinfo[biomol] = list(map(lambda x: x.replace(',', '').replace(';', ''), line.strip().split()[7:]))
                                if line.strip().split()[-1].endswith(','): #still there is something..
                                    for line in inp:
                                        pdbinfo[biomol] += list(
                                            map(lambda x: x.replace(',', '').replace(';', ''), #add chains
                                                line.strip().split()[4:]))
                                        if line.strip().split()[-1].endswith(','): # continue cycle
                                            pass
                                        else:
                                            break
                                break
                    elif line.startswith('REMARK 350'):
                        pass
                    else:
                        break
            elif line.startswith("SEQRES"):
                if yep == 0:
                    possible_problem['chain_problem'].append(i)
                    possible_problem['seqpos_problem'].append(i)
                    possible_problem['bad_structure_problem'].append(i)
                    print('Possible problem found: ', i, '\n')
                else:
                    biomolecules = refine(pdbinfo)
                    pdb_models[i] = [biomolecules, pdb_chains_interch]
                break


def look_for_chain(keys, model):
    for key in keys:
        if key in model:
            return key
    return 0

colnames =  inpdata.columns.values.tolist()
print(colnames)
pdb_grouped = inpdata.groupby('pdb_id')
pdbiti = 1
fixed_annotations = pd.DataFrame(columns=colnames)
print(fixed_annotations)
for ipdb in pdb_grouped:
    print("[", time.strftime("%c"), "|", pdbiti, "/", len(pdb_grouped), "]")
    print(ipdb[0], "- start scanning")
    #PART I. CHAIN PROBLEM
    tcr = list(map(str, (list(pd.unique(ipdb[1]['chain_tcr'])))))
    mhc_a = list(map(str, (list(pd.unique(ipdb[1]['chain_mhc_a'])))))
    mhc_b = list(map(str, (list(pd.unique(ipdb[1]['chain_mhc_b'])))))
    antigen = list(map(str, list(pd.unique(ipdb[1]['chain_antigen']))))
    #chains_in_seq = list(map(str, (list(pd.unique(ipdb[1]['chain_tcr'])) + list(pd.unique(ipdb[1]['chain_antigen'])))))
    chains_in_seq = tcr+antigen+mhc_a+mhc_b
    mhcchains_in_seq = mhc_a+mhc_b
    #chain in seq = chains in annotations.
    #pdb_models[ipdb[0]] - pdb file
    #pdb_models[ipdb[0]][0] - pdbinfo2 - chains in each model. [1] - pdb_chains_interch - chains that can be changed into each other
    if 'nan' in chains_in_seq:
        bad_pdbs['loss_of_data'].append(ipdb[0])
    else:
        is_ok = 0
        #bad_pdbs = {'chain_problem':[], 'seqpos_problem':[], 'bad_structure_problem':[], 'loss_of_data':[]}
        for i in pdb_models[ipdb[0]][0]:
            #if chains_in_seq are in biomolecule i
            if intersect(pdb_models[ipdb[0]][0][i], chains_in_seq) == sorted(chains_in_seq):
                is_ok = 1
                if len(pdb_models[ipdb[0]][1]) != 5:
                    bad_pdbs['bad_structure_problem'].append([ipdb[0], 'Found '+str(len(pdb_models[ipdb[0]][1]))+' chains, while should be 5.'])
                if len(pdb_models[ipdb[0]][0][i]) != 5:
                    bad_pdbs['possible_problems'].append([ipdb[0], 'Found ' + str(len(pdb_models[ipdb[0]][0][i])) +
                                                          ' chains in biomolecule, while should be 5.'])
                break
            elif intersect(pdb_models[ipdb[0]][0][i], chains_in_seq) == [] or intersect(tcr, pdb_models[ipdb[0]][0][i]) == sorted(pdb_models[ipdb[0]][0][i]):
                pass
            else:
                if intersect(tcr, pdb_models[ipdb[0]][0][i]) == sorted(tcr):
                    is_ok = 2
                    #now we are looking for antigen and (mhc_a or/and mhc_b).
                    antigen_in = 0
                    mhc_a_in = 0
                    mhc_b_in = 0

                    if intersect(mhc_a, pdb_models[ipdb[0]][0][i]) == sorted(mhc_a):
                        mhc_a_in = 1

                    if intersect(mhc_b, pdb_models[ipdb[0]][0][i]) == sorted(mhc_b):
                        mhc_b_in = 1

                    if intersect(antigen, pdb_models[ipdb[0]][0][i]) == sorted(antigen):
                        antigen_in = 1

                    for chainclass in pdb_models[ipdb[0]][1]:
                        if antigen_in == 0:
                            if is_in_list(antigen, chainclass):
                                outant = look_for_chain(chainclass, pdb_models[ipdb[0]][0][i])
                                if outant != 0:
                                    antigen_in = 2
                        if mhc_a_in == 0:
                            if is_in_list(mhc_a, chainclass):
                                mhc_a_new = look_for_chain(chainclass, pdb_models[ipdb[0]][0][i])
                                if mhc_a_new != 0:
                                    mhc_a_in = 2
                        if mhc_b_in == 0:
                            if is_in_list(mhc_b, chainclass):
                                mhc_b_new = look_for_chain(chainclass, pdb_models[ipdb[0]][0][i])
                                if mhc_b_new != 0:
                                    mhc_b_in = 2

                    if antigen_in == 2:
                        bad_pdbs['chain_problem'].append([ipdb[0], 'Change antigen from ' +
                                                          str(antigen).replace('[', '').replace(']', '') +
                                                          ' to \'' + str(outant)+'\''])
                        inpdata.set_value(inpdata['pdb_id'] == ipdb[0], 'chain_antigen', str(outant))
                    elif antigen_in != 1:
                        bad_pdbs['chain_problem'].append([ipdb[0], 'No antigen found in biomolecule'])

                    if mhc_a_in == 2:
                        bad_pdbs['chain_problem'].append([ipdb[0], 'Change mhc_a from ' +
                                                          str(mhc_a).replace('[', '').replace(']', '') +
                                                          ' to \'' + str(mhc_a_new) + '\''])
                        inpdata.set_value(inpdata['pdb_id'] == ipdb[0], 'chain_mhc_a', str(mhc_a_new))
                    elif mhc_a_in != 1:
                        bad_pdbs['chain_problem'].append([ipdb[0], 'No mhc_a found in biomolecule'])

                    if mhc_b_in == 2:
                        bad_pdbs['chain_problem'].append([ipdb[0], 'Change mhc_b from ' +
                                                          str(mhc_b).replace('[', '').replace(']', '') +
                                                          ' to \'' + str(mhc_b_new) + '\''])
                        inpdata.set_value(inpdata['pdb_id'] == ipdb[0], 'chain_mhc_b', str(mhc_b_new))
                    elif mhc_b_in != 1:
                        bad_pdbs['chain_problem'].append([ipdb[0], 'No mhc_b found in biomolecule'])

                    if len(pdb_models[ipdb[0]][1]) != 5:
                        bad_pdbs['bad_structure_problem'].append(
                            [ipdb[0], 'Found ' + str(len(pdb_models[ipdb[0]][1])) + ' chains, while should be 5.'])

                    if len(pdb_models[ipdb[0]][0][i]) != 5:
                        bad_pdbs['possible_problems'].append([ipdb[0], 'Found ' + str(len(pdb_models[ipdb[0]][0][i])) +
                                                              ' chains in biomolecule, while should be 5.'])

                    if antigen_in == 2 or mhc_a_in == 2 or mhc_b_in == 2:
                        print('fixing')
                        fixed_annotations = fixed_annotations.append(inpdata[inpdata['pdb_id'] == ipdb[0]])

        if is_ok == 0:
            bad_pdbs['chain_problem'].append([ipdb[0], 'Chaos in tcr chains and/or in antigen'])
            is_ok = 0

    #PART II. ANNOTATION PROBLEM
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

    work_structure = pdb_parser.get_structure("first", os.path.join(outp,'pdb', "pdb" + ipdb[0] + ".ent"))
    work_model = work_structure[0]

    for seq in ipdb[1].groupby('tcr_region_seq'):
        chain = list(pd.unique(seq[1]['chain_tcr']))[0]
        work_range = range(int(seq[1]['tcr_region_start']), int(seq[1]['tcr_region_end']))
        work_residues = get_residues(work_model[chain], work_range)
        work_seq = ''
        for i in work_residues:
            if i.get_id()[0] != ' ':
                work_residues.remove(i)
            else:
                work_seq += tto[str(i.get_resname())]
        if seq[0] != work_seq:
            if ipdb[0] not in bad_pdbs['seqpos_problem']:
                bad_pdbs['seqpos_problem'] = {ipdb[0]:[[seq[0], work_seq]]}
            else:
                bad_pdbs['seqpos_problem'][ipdb[0]].append([seq[0], work_seq])
    pdbiti+=1
    print("Done")
print("Finished scanning", len(pdb_grouped), " entries.")
with open(os.path.join(outp, 'problems.txt'), 'w') as out:
    for i in sorted(bad_pdbs):
        if i == 'chain_problem':
            out.writelines('Problems with chain:\n')
            for k in bad_pdbs[i]:
                out.writelines('\t'.join(k)+'\n')
            out.writelines('//\n\n')
        elif i == 'loss_of_data':
            out.writelines('Loss of data:\n'+'\n'.join(bad_pdbs[i])+'\n//\n\n')
        elif i == 'seqpos_problem':
            out.writelines('Problem with data annotation:\npdb\tshould_be\twhat_is\n')
            for k in bad_pdbs[i]:
                out.writelines('\t'.join(k)+'\n')
            out.writelines('//\n\n')
        elif i == 'bad_structure_problem':
            out.writelines('Problem with structures:\n')
            for k in bad_pdbs[i]:
                out.writelines('\t'.join(k)+'\n')
            out.writelines('//\n\n')
        elif i == 'possible_problems':
            out.writelines('Possible problems:\n')
            for k in bad_pdbs[i]:
                out.writelines('\t'.join(k)+'\n')
            out.writelines('//\n\n')
fixed_annotations.to_csv(output_file, sep='\t', index=None)
