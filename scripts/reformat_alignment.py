# updated file for extracting original template amino acid sequence
from os import write
from pathlib import Path
import subprocess
from sys import argv

def get_path(acc_code,dir):
    # gets the paths for the summary file (hhr)
    # gets the path for the original alignment
    p1String=acc_code+'.hhr'
    p2String=fasta+'_aln.fasta'
    path1=Path(dir+'/var/'+p1String)
    path2=Path(dir+'/var/'+p2String)
    return path1,path2

def get_hhr(path):
    # takes path as a parameter
    # returns the text from the summary
    with open(path,'r') as hhr:
        for line in enumerate(hhr):
            doc=hhr.readlines()
    return doc

def get_tmpl_header(hhr,dir):
    # taking the summary file of the input, it parses it
    # to extract the pdb chain for the top hit from the query
    x = range(0,len(hhr))
    top_template = ''
    for i in x:
        if hhr[i][0:3]==' No': # how the header just above the top template
        # begins    
            top_template=hhr[i+1]

    tt_list = top_template.split() # tokenized string

    pdb_id = tt_list[1].lower() # print id
    pdb_id = pdb_id.split('_')
    pdb = pdb_id[0]
    try:
        chain = pdb_id[1].upper()
    except:
        chain = ''
    # other potentially useful values
    #tpl_range = tt_list[-2] # = template range
    #tpl_offset = tpl_range.split('-')[0]
    #print(tt_list[-3]) # = target range
    #print(tt_list[-4]) # = cols
    #print(tt_list[-6]) # = scores
    #print(tt_list[-8]) # = e value
    with open(Path(dir+'/var/pdbID.txt'),'w') as f:
        id = ''.join([pdb,'.',chain])
        print(id)
        f.write(id)

    tmpl_header = '>'+pdb+'.'+chain.upper()+'\n'
    return tmpl_header

def get_aln_string(path,dir,target,header):
    # although we don't use this alignment,
    # we will keep it for testing purposes
    with open(path, 'r') as fasta:
        trg_index=0
        end_trg=0 
        tmpl_index = 0
        for line in enumerate(fasta):
            fa=fasta.readlines()
        r = range(0,len(fa))
        for int in r:
            # where target begins
            if fa[int][0:3]=='>sp' or fa[int][0:3]=='>tr':
                trg_index = int+1
            if fa[int][0:5]=='>ss_p':
                end_trg=int
            if fa[int][0:1]=='>' and (fa[int][1:2] != 's' and fa[int][1:3]!='Co'):
                tmpl_index=int+1
        trgList=fa[trg_index:end_trg]
        tmpList=fa[tmpl_index:]

        aln_string='>target\n'+''.join(trgList)+header+ ''.join(tmpList)
        with open(Path(dir+'/var/'+ target+'_target.fasta'),'w') as trg:
            trg_str = ''.join(trgList)
            trg.write(trg_str)

        # P63241tmpl_from_hhblits.fasta
        with open(Path(dir+'/var/'+target+'tmpl_from_hhblits.fasta'),'w') as f:
            tmp = ''.join(tmpList)
            f.write(tmp)
        return aln_string

def save_aln_file(aln_string,fasta,dir):
    # write out new alignment
    path=Path(dir+"/var/"+fasta+'_aln_2.fasta')
    with open(path,'w') as aln:
        aln.write(aln_string)
    print('Alignment written to file')

if __name__ == "__main__":
    dir=argv[1] # requires directory to write out file
    fasta = argv[2]
    pathTup=get_path(fasta,dir)
    path1=pathTup[0]
    path2=pathTup[1]

    hhr=get_hhr(path1)

    tmpl_header=get_tmpl_header(hhr,dir)
    aln_string = get_aln_string(path2,dir,fasta,tmpl_header)

    save_aln_file(aln_string,fasta,dir)
