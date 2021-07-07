from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from os import write              
from pathlib import Path
from sys import argv

#P63241tmpl_from_hhblits.fasta
#template_from_pdb.fasta
def get_paths(dir,fasta):
    p1String=fasta+'tmpl_from_hhblits.fasta'
    path=Path(dir+'/var/'+p1String)
    path2=Path(dir+'/var/template_from_PDB.fasta')
    return path,path2

def get_templates(path):
    fa=''
    with open(path, 'r') as fasta:
        fa=fasta.read()           
        fa = fa.replace('-','')
        
    return fa

def trim_left(og_fas,pdb_fas):
    og = list(og_fas)
    pdb = list(pdb_fas)
    co = 0
    while og[0] == '-':
        if len(og)==0:
           break 
        popped=og.pop(0)
        co += 1
    r = range(0,co)
    for i in r:
        pdb.pop(0)
    pdb = ''.join(pdb)
    og=''.join(og)
    return og,pdb,co

def trim_right(og_fas,pdb_fas):
    og_fas=og_fas.strip()
    pdb_fas=pdb_fas.strip()
    og = list(og_fas)
    pdb = list(pdb_fas)
    co = 0
    while og[-1] == '-':     
        if len(og)==0:
            break
        popped=og.pop(-1)    
        co += 1
    r = range(1,co)
    for i in r:
        pdb.pop(-1)
    pdb = ''.join(pdb)
    return pdb,co


def trim(og_fas,pdb_fas):
    print('trimming fasta files')
    left_results=trim_left(og_fas,pdb_fas)
    og = left_results[0]
    pdb=left_results[1]
    count=left_results[2]
    right_results = trim_right(og,pdb)
    pdb = right_results[0]
    right_deletions=right_results[1]
    return pdb,count,right_deletions

if __name__ == '__main__':
    dir=argv[1]
    fasta = argv[2]
    paths = get_paths(dir,fasta) 
    og_path=paths[0]
    pdb_path=paths[1]
    og_fasta = get_templates(og_path)
    pdb_fasta = get_templates(pdb_path)
    alignments = pairwise2.align.globalxx(og_fasta, pdb_fasta,one_alignment_only=True)
    #print(format_alignment(*alignments[0]))
    r = range(0,len(alignments))
    for i in r:
        x=''
    #    print(format_alignment(*alignments[i]))
#    print(alignments)
#    print(alignments[0].seqA)

# if the orignial is shorter than the pdb, the pdb fasta will need to be trimmed. We will use this format (it does not however, fix internal spaces.) to do that, we will remove all spaces, and then align it to the old file.
    og_aligned=''
    pdb_alligned=''
    # need to return otherwise
    pdb_aligned=alignments[0].seqB
    if len(og_fasta) <= len(pdb_fasta):
        og_aligned = alignments[0].seqA
        pdb_aligned = alignments[0].seqB
        trimmed = trim(og_aligned,pdb_aligned)
        pdb=trimmed[0]
        count=trimmed[1]
        right_count=trimmed[2]

    p=Path(dir+'/var/'+fasta+'_template.fasta')
    #remove - from pdb
    pdb_aligned=pdb_aligned.replace('-','')
    with open(p,'w') as file:
        file.write(pdb_aligned)
