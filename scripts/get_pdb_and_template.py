from ost import io, seq
from promod3 import modelling, loop
from pathlib import Path
from sys import argv
# this will work to fetch the pdb file, with the name of the pdb passed in the dommand line
def get_id(dir):
    path = dir+'/var/pdbID.txt'
    pdb = ''
    with open(path,'r') as p:
       pdb=p.read()
    pdb = pdb.split('.')
    return pdb

def get_fasta(pdb,dir):
    #for res in pdb.residues:
    #print(''.join([res.one_letter_code for r in pdb.residues]))
    # for some reason, this does not work
    res_list=[]
    for res in pdb.residues:
        o_res = res.one_letter_code
        if o_res != '?':
            res_list.append(o_res)
    res_string = ''.join(res_list)
    path = Path(dir+'/var/template_from_PDB.fasta')
    f = open(path,'w')
    f.write(res_string)
    f.close()

def get_pdb(pdb_tup,dir):
    if len(pdb_tup)!=2:
        id = pdb_tup[0]
        chain = ''
    else:
        id = pdb_tup[0]
        chain = pdb_tup[1]

    p=io.LoadPDB(id,seqres=True,remote=True,remote_repo='pdb')
    pdb=p[0]

    pdb_path = dir+'/var/template.pdb'
    if len(chain) != 0:
        l = pdb.GetChainList()
        print('Saving {} chain {} to file'.format(id,chain))
        query = 'chain='+chain
        pdb_crn = pdb.Select(query) # returns entityview object
        io.SavePDB(pdb_crn,pdb_path)
        pdb=io.LoadPDB(pdb_path)
        #get_fasta(p,dir)

    else:
        io.SavePDB(pdb,pdb_path)
        #get_fasta(pdb,dir)
    #for res in pdb.residues:                
    #    print(res.one_letter_code)
    #print(''.join([res.one_letter_code for r in pdb.residues]))
    '''
    test part
    '''
    res_list=[]
    for res in pdb.residues:
        o_res = res.one_letter_code
        if o_res != '?':
            res_list.append(o_res)
    res_string = ''.join(res_list)
    path = Path(dir+'/var/template_from_PDB.fasta')
    f = open(path,'w') 
    f.write(res_string)
    f.close()
    print('pdb saved to file')


if __name__ == '__main__':

    dir = argv[1]
    pdb_id = get_id(dir)
    print('Getting pdb from database...')
    pdb = get_pdb(pdb_id,dir)
    print('extracting sequence from template...')
