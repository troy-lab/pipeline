cpu=4

# get fasta file from uniprot
target=$1
# P08133,,,,
# P07737,,,,
# P31946,,,,
# P09651,,,,
# base data
# Database base directory
database=$2 #/root/pipeline/database/pdb70
# base url for getting amino acid sequence for target protein
url=https://www.uniprot.org/uniprot
rootDir= $3 #/root/pipeline run files from root, store in data
dataDir=${rootDir}/data/${target}
wget $url/${target}.fasta -P ${dataDir}
mkdir ${dataDir}/var

# 1. Query database using hhblits
# first, hhblits uses the fasta sequence to query the database
# it outputs an alignment file, and a summary file (.hhr
hhblits -id 100 -cov 50 -cpu $cpu -i $dataDir/${target}.fasta -hide_dssp -B 1 -b 1 -Ofas ${dataDir}/var/${target}_aln.fasta -o ${dataDir}/var/${target}.hhr -d $database

# 2. reformat the alignment for check
python3 reformat_alignment.py $dataDir $target

# 3. get .pdb structure for template and extract amino acid sequence
# requires the promod build and singularity
singularity run --app PM $HOME/pipeline/promod.sif get_pdb_and_template.py $dataDir

# 4. align orginal template to extracted template
python3 align_template1_2.py $dataDir $target

# 5. align target to new template
python3 final_align.py $dataDir $target

#singularity run -- app PM $HOME/pipeline/promod.sif 
# 6. build model
singularity run --app PM $HOME/pipeline/promod.sif build-model -f ${dataDir}/final_aln.fasta -p ${dataDir}/var/template.pdb -o ${dataDir}/model.pdb

# Calculate DSSP
python3 get_dssp.py $dataDir

