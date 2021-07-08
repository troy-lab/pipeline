import subprocess
from pathlib import Path
import os
import urllib.request
import json

from ost import io, seq
from promod3 import modelling, loop
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

# Check 2. reformat the alignment for check
#python3 reformat_alignment.py $dataDir $target

# Check 3. get .pdb structure for template and extract amino acid sequence
# requires the promod build and singularity
#singularity run --app PM $HOME/pipeline/promod.sif get_pdb_and_template.py $dataDir

# check 4. align orginal template to extracted template 
#python3 align_template1_2.py $dataDir $target

# check 5. align target to new template
#python3 final_align.py $dataDir $target

#singularity run -- app PM $HOME/pipeline/promod.sif 
# check 6. build model
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
        pdb_id = self.build_dict['pdb_id']
        p=io.LoadPDB(pdb_id,seqres=True,remote=True,remote_repo='pdb')
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


    def extract_from_aln(self):
        chain = self.build_dict['chain_id']
        if chain != '':
            header = '>'+self.build_dict['pdb_id']+'.'+chain.upper()+'\n'
        else:
            header = '>'+self.build_dict['pdb_id'] 
        # path to .aln output file
        path = self.misc_dir + '/'+self.accession_code + '_aln.fasta'
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


            tmp = ''.join(tmpList)
            tmp = tmp.replace('\n','')
            tmp = tmp.replace('-','')
            self.build_dict['template_fasta_from_aln']=tmp

            trg = ''.join(trgList)
            trg = trg.replace('\n','')
            trg = trg.replace('-','')
            self.build_dict['target_fasta_from_aln']=trg

    # for trimming alignment fasta
    def trim_left(self,og_fas,pdb_fas):
        og = list(og_fas)
        pdb = list(pdb_fas)
        co = 0
        while og[0] == '-':
            if len(og)==0:
                break 
            popped=og.pop(0)
            co += 1
        r = range(0,co)
        raw_count = co
        co = 0
        for i in r:
            popped = pdb.pop(0)
            if popped != '-':
                co +=1
        pdb = ''.join(pdb)
        og=''.join(og)
        return og,pdb,co

    # for trimming alignment fasta
    def trim_right(self,og_fas,pdb_fas):
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


    # for trimming alignment fasta
    def trim(self,og_fas,pdb_fas):
        print('trimming fasta files')
        left_results=self.trim_left(og_fas,pdb_fas)
        og = left_results[0]
        pdb=left_results[1]
        count=left_results[2]
        right_results = self.trim_right(og,pdb)
        pdb = right_results[0]
        right_deletions=right_results[1]
        return pdb,count,right_deletions

    def fix_templates(self):
        og_template = self.build_dict['template_fasta_from_aln']
        og_template = og_template.replace('-','')
        pdb_template = self.build_dict['template_from_pdb']
        pdb_template = pdb_template.replace('-','')
        alignments = pairwise2.align.globalms(og_template, pdb_template,3,1,-3,-1,one_alignment_only=True)
        print(format_alignment(*alignments[0]))
        print(alignments)
        print(alignments[0])
        pdb_aligned = alignments[0][0]
        print(pdb_aligned)

        if len(og_template) < len(pdb_aligned):
            og_aligned = alignments[0][1]
            print('og aligned')
            print(og_aligned)
            trimmed = self.trim(og_aligned,pdb_aligned) # og template, then pdb_template
            pdb_aligned = trimmed[0]
            count = trimmed[1]
            self.build_dict['offset']=count

        pdb_aligned = pdb_aligned.replace('-','')

        self.build_dict['final_template_fasta'] = pdb_aligned

    def make_final_alignment(self):
        template = self.build_dict['final_template_fasta']
        target = self.build_dict['target_fasta_from_aln']
        aln = pairwise2.align.globalms(target, template, 5, -.5, -4, -2,one_alignment_only=True)
        print(format_alignment(*aln[0]))
        file_string='>target\n'+aln[0][0] + '\n>template\n'+aln[0][1]
        align_path = Path(self.misc_dir+'/final_aln.fasta')
        with open(align_path, 'w+') as f:
            f.write(file_string)

    def create_model(self):
        # get raw model
        tpl = io.LoadPDB(self.misc_dir+'/template.pdb')
        aln = io.LoadAlignment(self.misc_dir+'/final_aln.fasta')
        aln.AttachView(1, tpl.CreateFullView())
        mhandle = modelling.BuildRawModel(aln)

        # build final model
        final_model = modelling.BuildFromRawModel(mhandle)
        io.SavePDB(final_model, self.output_dir+'/model.pdb')

    def get_dssp(self):
        p = PDBParser()
        structure = p.get_structure(self.accession_code,self.output_dir+'/model.pdb')
        model = structure[0]
        dssp = DSSP(model, self.output_dir+'/model.pdb',dssp='mkdssp')
        a_key = list(dssp.keys())[2]
        dssp_list = []
        for key in dssp.keys():
            dssp_list.append(dssp[key][2])
        secondary_structure = ''.join(dssp_list)
        self.build_dict['dssp']=secondary_structure





if __name__ == '__main__':

    x = protein('P06733')
    try:
        x.get_target_fasta()
    except (GetFastaFailed, WriteFileFailed) as e:
        pass

    x.query_hhblits()

    x.extract_from_query()

    x.extract_from_aln()

    x.get_pdb()

    x.fix_templates()

    x.build_to_json()

    x.make_final_alignment()

    x.create_model()

    x.get_dssp()

    x.build_to_json()
