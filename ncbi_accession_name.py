#!/usr/bin/env python3

"""
Uses NCBI's `datasets-cli` tool to identify accession numbers for genome assemblies
while updating taxonomic names.

Dependencies:
- NCBI's datasets-cli

Author: Xyrus X. Maurer-Alcalá
Last Updated: 2026-09-17
"""

import subprocess, sys
import pandas as pd

from collections import defaultdict

def ncbi_dsets(
        taxon_name: str,
        acc: str = None,
        curr_attempt: int = 1,
        max_attempts: int = 5):
    dsets_cmd = ['datasets', 'summary', 'genome', 'taxon', taxon_name, '--reference']
    result = subprocess.run(dsets_cmd, capture_output = True, text = True)
    acc = result.stdout.partition('"accession":')[-1].partition(',')[0].replace('"','')
    if not acc and 'taxids:' in result.stderr:
        curr_attempt += 1
        if curr_attempt <= max_attempts:
            taxon_name = result.stderr.partition('taxids:\n')[-1].partition(' (')[0]
            taxon_name, acc, curr_attempt = ncbi_dsets(taxon_name, curr_attempt = curr_attempt)
    elif not acc:
        dsets_cmd = ['datasets', 'summary', 'genome', 'taxon', taxon_name]
        result = subprocess.run(dsets_cmd, capture_output = True, text = True)
        acc = result.stdout.partition('"accession":')[-1].partition(',')[0].replace('"','')
    return taxon_name, acc, curr_attempt


def check_name(ncbi_taxon_name: str):
    tname = ncbi_taxon_name.split(' ')
    if ncbi_taxon_name.startswith('Candidatus'):
        name = ' '.join(tname[:3])
        strain = ' '.join(tname[3:])
    else:
        name = ' '.join(tname[:2])
        strain = ' '.join(tname[2:])
    return name, strain

def parse_table(input_csv, curr_attempt = 1, max_attempts = 5):
    fin_dict = defaultdict(dict)

    for line in open(input_csv).readlines():
        taxon_name = line.split(',')[0].rstrip()
        print(f'Processing {taxon_name}')
        ncbi_taxon_name, acc, curr_attempt = ncbi_dsets(
                                            taxon_name,
                                            curr_attempt = 1,
                                            max_attempts = 5)
        ncbi_name, ncbi_strain = check_name(ncbi_taxon_name)
        fin_dict[taxon_name] = {'NCBI_Name': ncbi_name, 'Strain': ncbi_strain, 'NCBI_Accession':acc}


    acc_df = pd.DataFrame(fin_dict).T
    acc_df.index.name = 'Taxon'
    acc_df.to_csv(f'{input_csv.rpartition(".")[0]}.NCBI_Accessions.csv')

if __name__ == '__main__':
    try:
        input_csv = sys.argv[1]
    except:
        print('\n\nUsage:\n\n    python3 ncbi_accession_name.py [TAXON-TABLE]\n\n')
        sys.exit(1)

    parse_table(input_csv)
