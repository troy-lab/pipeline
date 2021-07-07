from pathlib import Path
import subprocess
import pandas as pd
print('success')

#target=$1, find from .csv
# P08133,,,,
# P07737,,,,
# P31946,,,,
# P09651,,,,
# base data
# Database base directory
# base url for getting amino acid sequence for target protein

def accession_from_csv(filename):
    pass

def run_pipeline(acc_code, db, root):
    run_command = 'time sh pipeline.sh {} {} {}'.format(acc_code, db, root)
    print(run_command)
    subprocess.call(run_command,shell=True)

if __name__ == '__main__':
    database = './database/pdb70' # change to database path
    rootDir = '.' # change to root of pipeline
    path_to_csv = ''

    run_pipeline('P09651', database, rootDir)




