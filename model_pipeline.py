import subprocess
from pathlib import Path
import os
import urllib.request
import json
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
        s = subprocess.Popen(cmd)

if __name__ == '__main__':
    x = protein('P53621')
    try:
        x.get_target_fasta()
    except (GetFastaFailed, WriteFileFailed) as e:
        pass

    x.query_hhblits()

    x.build_to_json()
