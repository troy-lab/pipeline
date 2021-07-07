from os import write              
from pathlib import Path
import subprocess
from sys import argv 

def get_path(fasta,dir):
    pString=fasta+'_aln.fasta'
    path=Path(dir+'/'+pString)
    return path

def get_aln_string(path):
    with open(path, 'r') as fasta:
        tmpl_index = 0
        for line in enumerate(fasta):
            fa=fasta.readlines()
        r = range(0,len(fa))
        for int in r:
            if fa[int][0:1]=='>' and (fa[int][1:2] != 's' and fa[int][1:3]!='Co'):
                tmpl_index=int+1
        tmpList=fa[tmpl_index:]

        aln_fasta=''.join(tmpList)
        return aln_fasta

def save_aln_file(aln_string,fasta,dir):
    path=Path(dir+"/"+fasta+'tmpl.fasta')
    with open(path,'w') as aln:
        aln.write(aln_string)
    print('Alignment written to file')

if __name__ == "__main__":

    dir=argv[1]
    fasta = argv[2]   

    path = get_path(fasta,dir)
    aln=get_aln_string(path)
    save_aln_file(aln,fasta,dir)
