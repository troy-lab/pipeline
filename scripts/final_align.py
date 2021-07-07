from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from os import write              
from pathlib import Path
from sys import argv

# include offset here

def get_paths(dir,fasta):
    p1String=fasta+'_template.fasta'
    path=Path(dir+'/var/'+p1String)
    path2=Path(dir+'/var/'+fasta+'_target.fasta')
    print(path2)
    return path,path2

# Q8N137_target.fasta

if __name__ == '__main__':
    dir=argv[1]
    fasta=argv[2]

    paths=get_paths(dir,fasta)
    temp_path = paths[0]
    trg_path=paths[1]
    target = ''

    #get target fasta
    with open(trg_path,'r') as trg:
        target=trg.read()
    print(target)
    target=target.replace('\n','')

    print('target:{}'.format(target))
    target = target.replace('-','')
    

    #get template fasta
    template=''
    with open(temp_path, 'r') as tmp:
        template=tmp.read()
    print('template:{}'.format(template))
    aln = pairwise2.align.globalms(target, template, 5, -.5, -4, -2,one_alignment_only=True)

    print(format_alignment(*aln[0]))

    file_string='>target\n'+aln[0].seqA+'\n>template\n'+aln[0].seqB
    align_path = Path(dir+'/final_aln.fasta')
    with open(align_path, 'w') as f:
        f.write(file_string)
