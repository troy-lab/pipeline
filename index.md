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

If your .csv file is not named "input.csv", you will need to edit line 371 of model_pipeline.py accordingly.

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

This will contain how to format accession .csv files and where the place to change the defaults is found in the source code

### References

Steinegger M, Meier M, Mirdita M, Vöhringer H, Haunsberger S J, and Söding J (2019) HH-suite3 for fast remote homology detection and deep protein annotation, BMC Bioinformatics, 473. doi: 10.1186/s12859-019-3019-7

Studer G, Tauriello G, Bienert S, Biasini M, Johner N, Schwede T (2021) ProMod3—A versatile homology modelling toolbox. PLoS Comput Biol 17(1): e1008667. https://doi.org/10.1371/journal.pcbi.1008667


This will be the places my software is from, such as HHBLITS, ProMod3, etc.

### Special Thanks

Thank you to Dr. Bakhtiyor Rasulev for supporting and guiding this project.

Thank you also to Dr. Dali Sun for your support of this project.

Thank you to NDSU's Coatings and Polymeric Materials department.

For more information, please visit the following:

https://rasulev.org/

https://www.ndsu.edu/cpm/

Any more information?
Include Grant information, thank you's, etc. here.


### Github stuff

You can use the [editor on GitHub](https://github.com/troy-lab/pipeline/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.
Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/troy-lab/pipeline/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
