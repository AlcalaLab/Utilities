#!/usr/bin/env python3

import glob, os, shutil, subprocess, sys
from pathlib import Path

import pandas as pd


def prep_folder(folder_name: str):
    Path(folder_name).mkdir(exist_ok = True, parents = True)


def datasets_cli(accession: str, taxon_dir: str, taxon_name: str, files_to_grab: list):
    datasets_cmd = f'datasets download genome accession {accession} ' \
        f'--filename {taxon_dir}{taxon_name}.NCBI_Dataset.zip'

    datasets_cmd += f' --include {"".join(files_to_grab)}'

    print(datasets_cmd)

    try:
        dsets_rslt = subprocess.run(datasets_cmd.split(), check = True, capture_output = True, text = True)

    except subprocess.CalledProcessError as e:
        print(f'Failed to download acc: {accession} for {taxon_name}')

    return f'{taxon_dir}{taxon_name}.NCBI_Dataset.zip'


def unzip_prep(dset_zip: str, proj_dir: str, taxon_name: str, files_to_grab: list = ['cds']):
    unzip_cmd = f'unzip {dset_zip}'
    # print(unzip_cmd)
    unzip_rslt = subprocess.run(unzip_cmd.split(), check = True, capture_output = True, text = True)
    # print(unzip_rslt.stderr)
    for f in files_to_grab:
        if f == 'cds':
            try:
                query_file = glob.glob('ncbi_dataset/data/*/*cds_*.fna')[0]
                print(query_file)
            except IndexError:
                continue
            shutil.copy2(query_file, f'{proj_dir}/CDS_Files/{taxon_name}.WGS.CDS.fasta')
    os.system('rm -rf ncbi_dataset')
    os.system('rm md5sum.txt README.md')


def download_table_data(input_table: str, proj_dir: str, files_to_grab: list = ['cds']):

    df = pd.read_csv(input_table)

    for file_type in files_to_grab:
        Path(f'{proj_dir}{file_type.upper()}_Files/').mkdir(exist_ok = True, parents = True)

    for index, row in df.iterrows():
        taxon_name = row['Taxon'].rstrip().rstrip('.').replace(' ','_').replace('.','_')
        taxon_dir = f'{proj_dir}{taxon_name}_NCBI_Datasets/'

        prep_folder(taxon_dir)

        print(taxon_name)

        dset_zip = datasets_cli(row['Accession(s)'], taxon_dir, taxon_name, files_to_grab)
        # print(dset_zip, proj_dir, taxon_name, files_to_grab)
        unzip_prep(dset_zip, proj_dir, taxon_name, files_to_grab)



if __name__ == '__main__':
    try:
        input_table = sys.argv[1]
        proj_name = sys.argv[2]
    except:
        print('\nUsage:\n\n    python3 ncbi_genomes_prep.py [TABLE-OF-GENOMES] [OUTPUT-DIRECTORY]\n')
        sys.exit()

    files_to_grab = ['cds']

    proj_dir = f'{proj_name}_NCBI_Datasets/'


    keep_intermed = False

    download_table_data(input_table, proj_dir, files_to_grab)
