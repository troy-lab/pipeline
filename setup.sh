# first load singularity - note, not necessary for certain systems
module load singularity
# pull proper singularity container
singularity pull library://ttimmerman/default/modeling_engine:sha256.87e8cec025b097430cf0ebda124a8eec1f640603657dbf6c2059cc05a3fbc4c3
# rename container
mv modeling_engine* promod.sif

# download hhblits compressed database
wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb70_from_mmcif_latest.tar.gz
# make directory for database
mkdir database
# move compressed file to the database
mv *.tar.gz ./database/
# change to database directory
cd database
# untar compressed database file
tar xvf *.tar.gz 
# if the build was successful, the compressed file can be deleted
rm *.tar.gz
