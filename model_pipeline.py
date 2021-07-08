import subprocess
from pathlib import Path
import os
import urllib.request
import json
from ost import io, seq
from promod3 import modelling, loop

# 2. reformat the alignment for check
#python3 reformat_alignment.py $dataDir $target

# 3. get .pdb structure for template and extract amino acid sequence
# requires the promod build and singularity
#singularity run --app PM $HOME/pipeline/promod.sif get_pdb_and_template.py $dataDir

# 4. align orginal template to extracted template 
#python3 align_template1_2.py $dataDir $target

# 5. align target to new template
#python3 final_align.py $dataDir $target

#singularity run -- app PM $HOME/pipeline/promod.sif 
# 6. build model
#singularity run --app PM $HOME/pipeline/promod.sif build-model -f ${dataDir}/final_aln.fasta -p ${dataDir}/var/template.pdb -o ${dataDir}/model.pdb

# Calculate DSSP
#python3 get_dssp.py $dataDir


class GetFastaFailed(Exception):
    pass

class WriteFileFailed(Exception):
    pass

class protein:
    # class variables
    root_dir = '.'
    
    database = '/root/pipeline/database/pdb70' 
    
    def __init__(self, acc_code):
        self.accession_code = acc_code.strip()
        self.build_dict = {}
        # make directory for build
        self.output_dir = self.root_dir +'/'+ self.accession_code
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        # make directory for misc build files
        self.misc_dir = self.output_dir+'/var'
        if not os.path.exists(self.misc_dir):
            os.mkdir(self.misc_dir)

    def build_to_json(self):
        json_object = json.dumps(self.build_dict, indent = 4)
        file = Path(self.misc_dir,'{}.json'.format(self.accession_code))
        with open(file,'w+') as f:
            f.write(json_object)


    # use wget to retrieve the amino acid
    # sequence in fasta format
    # from uniprotKB
    def get_target_fasta(self):
        # url for fasta sequence
        fasta_url ='https://www.uniprot.org/uniprot/{}.fasta'.format(self.accession_code)
        trg_fasta = ''
        try:
            with urllib.request.urlopen(fasta_url) as url:
                trg_fasta = url.read().decode('utf-8')
            if len(trg_fasta) != 0:
                self.build_dict['target_fasta'] = trg_fasta
        except:
            raise GetFastaFailed
            print('get failed')
        try:
            # write target to path for querying
            fasta_path =self.output_dir+'/'+ self.accession_code+'.fasta' 
            self.build_dict['target_fasta_path'] = fasta_path
            with open(fasta_path,'w') as f:
                f.write(trg_fasta)
        except:
            raise WriteFileFailed


    # queries hhblits for template hit
    # saves query output file
    # adds template hit to build_dict
    def query_hhblits(self):
        #hhblits -id 100 -cov 50 -cpu $cpu -i $dataDir/${target}.fasta -hide_dssp -B 1 -b 1 -Ofas ${dataDir}/var/${target}_aln.fasta -o ${dataDir}/var/${target}.hhr -d $database

        aln_path = self.misc_dir + '/'+self.accession_code + '_aln.fasta'
        hhr_path = self.misc_dir + '/'+self.accession_code + '.hhr'
        print('alnpath{}\nhhrpath{}'.format(aln_path,hhr_path))
        cmd = ['/usr/local/share/hhblits/hh-suite/build/bin/hhblits', '-id','100','-cov','50','-i',self.build_dict['target_fasta_path'], '-hide_dssp', '-B', '1', '-b', '1', '-Ofas', aln_path, '-o', hhr_path, '-d', self.database]# path for HHBLITS for now  
        p = subprocess.Popen(cmd)
        (output, err) = p.communicate()  

        #This makes the wait possible
        p_status = p.wait()

        #This will give you the output of the command being executed
        print("Command output: ")
        print(output)

    def extract_from_query(self):
        hhr = ''
        # path to .hhr output file
        path = self.misc_dir + '/'+self.accession_code + '.hhr'
        with open(path,'r') as file:
            for line in enumerate(file):
                hhr=file.readlines()
        
        x = range(0,len(hhr))
        top_template = ''


        # taking the summary file of the input, it parses it
        # to extract the pdb chain for the top hit from the query
        for i in x:
            if hhr[i][0:3]==' No': # how the header just above the top template
            # begins    
                top_template=hhr[i+1]

        tt_list = top_template.split() # tokenized string

        pdb_id = tt_list[1].lower() # print id
        pdb_id = pdb_id.split('_')
        pdb = pdb_id[0]
        self.build_dict['pdb_id'] = pdb
        try:
            chain = pdb_id[1].upper()
        except:
            chain = ''

        self.build_dict['chain_id']= chain
    # other potentially useful values
        tpl_range = tt_list[-2] # = template range
        self.build_dict['template_range'] = tpl_range
    #tpl_offset = tpl_range.split('-')[0]
    #print(tt_list[-3]) # = target range
    #print(tt_list[-4]) # = cols
    #print(tt_list[-6]) # = scores
        e_value = tt_list[-8] # = e value
        self.build_dict['e_value']=e_value


    def get_pdb(self):
        p=io.LoadPDB(id,seqres=True,remote=True,remote_repo='pdb')
        pdb=p[0]
        pdb_path = self.misc_dir + '/template.pdb'
        chain = self.build_dict['chain_id']
        if chain != '':
            l = pdb.GetChainList()
            query = 'chain='+chain
            pdb_crn = pdb.Select(query) # returns entityview object
            io.SavePDB(pdb_crn,pdb_path)
            pdb=io.LoadPDB(pdb_path)
        
        else:
            io.SavePDB(pdb,pdb_path)

        res_list=[]
        for res in pdb.residues:
            o_res = res.one_letter_code
            if o_res != '?':
                res_list.append(o_res)
        res_string = ''.join(res_list)
        self.build_dict['template_from_pdb']=res_string



if __name__ == '__main__':
    x = protein('P06733')
    try:
        x.get_target_fasta()
    except (GetFastaFailed, WriteFileFailed) as e:
        pass

    x.query_hhblits()

    x.extract_from_query()

    x.build_to_json()
