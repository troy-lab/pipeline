## Homology Modeling Pipeline

This homology modeling pipeline is designed to streamline and automate the task of creating 3-D models in a .pdb format for a given UniProtKB accession code. Other than remotely accessing UniProtKB in order to access the amino acid sequence and retrieving the template .pdb structure from the Protein Data Bank, the modeling process will take place entirely on your local machine.


### Requirements

This software is only compatible with Linux systems.
It has been tested on kernel 4.18.0-193.el8.x86_64 using Singularity 3.7+.

To download Singularity 3.8, visit their documentation:
https://sylabs.io/guides/3.8/user-guide/quick_start.html

### Getting Started

To setup the software, begin by cloning the git repository.

```markdown
git clone https://github.com/troy-lab/pipeline.git
```

Next, to setup the database and necessary directories, within the git repository run the setup.sh script.
This step involves downloading the database of amino acid sequences from the protein databank, and may take several hours.

```markdown
sh setup.sh
```

Lastly, load a .csv file into the same directory. This should contain the UniprotKB Accession codes for the proteins you wish to model. Make sure to follow the proper formatting for the .csv file, as described below.

### Running the Software

So long as the format of the .csv files is correct, simply issue this command from the cloned repository:

```markdown
singularity exec --app PM promod.sif ~/pipeline/model
```
Note that if you did not place your repository in "$HOME/pipeline", you will need to adjust the path "~/pipeline/model" accordingly.

If you wish to change the name of the output.csv file, edit the last line of model_pipeline.py

```markdown
df.to_csv(x.root_dir+'/output.csv')
```
to 
```markdown
df.to_csv(x.root_dir+'/desired_name.csv')
```

The model file is a bash script which calls the python script model_pipeline.py. This is to simplify running the Homology Pipeline on an HPC cluster.

You may also run the software from the base directory like so:

```markdown
singularity run --app PM promod.sif model_pipeline.py

```


#### Running on HPC Clusters
When running the Homology Pipeline on an HPC cluster using Singularity users must run "singularity exec ..." rather than "singularity run ..."

To use the pipeline with the job scheduler qsub, users must call a bash script which then calls the python script. Be sure to allow execution access to the relevant files in the qsub .pbs file. Here is an example for qsub below:

```markdown
#!/bin/bash
#PBS -q default
#PBS -N name of instance
#PBS -M user@email
#PBS -m abe
#PBS -l select=1:mem=64gb:ncpus=8
#PBS -l walltime=1000:00 # adjust to time needed for models (maximum is about 3:30 per model)

cd $PBS_O_WORKDIR

chmod +x path/to/model_pipeline.py

echo -------------

chmod 777 path/to/.csv/file

chmod +x path/to/model (bash file)

module load singularity # if needed for your HPC cluster

time singularity exec --app PM /path/to/promod.sif /path/to/model

exit 0
```

#### Input Formatting Requirements

In the base folder of the git repository, include a .csv file with the UniProtKB accession codes for the proteins you wish to model. They should be under a column with the header "Accession." (N.b. this is case sensitive).

If you include a column "Checked", it should contain boolean values which are False for unchecked proteins.

The following columns are optional, but if included must be of the proper type:

'helix' -- float
'beta' -- float
'c' -- float
'unknown' -- float
'error' -- float

If your .csv file is not named "input.csv", you will need to edit line 371 of model_pipeline.py accordingly.

I.e.
```markdown
    df = pd.read_csv(root_dir+'/input.csv')
```

should read
```markdown
    df = pd.read_csv(root_dir+'/your_file.csv')
```

This will contain how to format accession .csv files and where the place to change the defaults is found in the source code

### Output

After running the software, a .csv file called output.csv will be in the base git repository. It contains the columns 'helix', 'beta', 'c', and 'unknown' which track the number of each secondary structure identified in the final .pdb model.
The 'Checked' column for each protein will be True. In addition, if the error column holds a value, this indicates an error in the modeling process for that protein.

In addition to output from the .csv file, there is a directory called 'data' which contains directories for each protein which was modeled. These directories contin the final model ('model.pdb') and the amino acid sequence for the target protein in a .fasta file. In the subdirectory .var, there are two alignment files with the extension '.aln', which contain the original and formatted alignments. A .hhr file shows the results of the hhblits query, and the template.pdb file has the coordinates for the 3D template used to build the final model. Lastly, a .json file contains information used in the build as well as the calculated DSSP values. If there was an error in modeling, the details can be found in the .json file as well. (n.b. if there was no error modeling a protein, no error information will be present).

### References

Steinegger M, Meier M, Mirdita M, Vöhringer H, Haunsberger S J, and Söding J (2019) HH-suite3 for fast remote homology detection and deep protein annotation, BMC Bioinformatics, 473. doi: 10.1186/s12859-019-3019-7

Studer G, Tauriello G, Bienert S, Biasini M, Johner N, Schwede T (2021) ProMod3—A versatile homology modelling toolbox. PLoS Comput Biol 17(1): e1008667. https://doi.org/10.1371/journal.pcbi.1008667

A series of PDB related databases for everyday needs.
Wouter G Touw, Coos Baakman, Jon Black, Tim AH te Beek, E Krieger, Robbie P Joosten, Gert Vriend.
Nucleic Acids Research 2015 January; 43(Database issue): D364-D368.

Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features.
Kabsch W, Sander C,
Biopolymers. 1983 22 2577-2637.
PMID: 6667333; UI: 84128824.


This will be the places my software is from, such as HHBLITS, ProMod3, etc.

### Special Thanks

Thank you to Dr. Bakhtiyor Rasulev for supporting and guiding this project.

Thank you also to Dr. Dali Sun for your support of this project.

Thank you to the Electrical & Computer Engineering department and the Coatings and Polymeric Materials department at NDSU for supporting this project.

For more information, please visit the following:

https://rasulev.org/

https://www.ndsu.edu/cpm/

Any more information?
Include Grant information, thank you's, etc. here.
