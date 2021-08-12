# pipeline
This is a project which takes input from a .csv file containing uniprot kb accession codes, and creates a model for each given code.

The modeling process begins by retrieving the .fasta sequence for the corresponding accession code, which is then queryed against the pdb70 database for homologs.
The top matching homolog is selected, and an alignment is generated for this homolog.
Next, the template structure for the template is collected from the Protein Databank. This is trimmed to match the alignment.
Finally, a model is generated using the ProMod3 modelling engine.
Using DSSP, the secondary structure for this model can be calculated.
This then can be added to a .csv file containing the accession codes.

For help running this software, please see the documentation:
https://troy-lab.github.io/pipeline/
