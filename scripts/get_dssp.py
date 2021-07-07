from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from pathlib import Path
from sys import argv

def get_id(dir):
    path = dir+'/var/pdbID.txt'
    pdb = ''
    with open(path,'r') as p:
       pdb=p.read()
    return pdb

if __name__ == '__main__':
    dir = argv[1]
    path = dir+'/model.dssp'
    pdb = get_id(dir)
    p = PDBParser()
    structure = p.get_structure(pdb,dir+'/model.pdb')
    model = structure[0]
    dssp = DSSP(model, dir+'/model.pdb',dssp='mkdssp')
    a_key = list(dssp.keys())[2]
    dssp_list = []
    with open(path,'w') as f:
        for key in dssp.keys():
            dssp_list.append(dssp[key][2])
        ss=''.join(dssp_list)
        f.write(ss)
