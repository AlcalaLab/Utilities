#!/usr/bin/env python3

import glob
import os
import shutil
import sys

from collections import defaultdict
from pathlib import Path

def download_sra_file(
        sra_accession: str
        ) -> None:

    prefetch_cmd = f'prefetch -v -X 40G {sra_accession}'

    os.system(prefetch_cmd)


def sra_fastq_dump(
        taxon: str,
        sra_accession: str,
        threads: int = 8
        ) -> str:

    taxon = taxon.replace(' ','_').replace('.','').replace('-','_')

    fastq_out_dir = f'{taxon}_{sra_accession}'

    fastq_dump_cmd = f'parallel-fastq-dump ' \
                    '--split-files ' \
                    '--gzip ' \
                    f'--threads {threads} ' \
                    f'-s {sra_accession} ' \
                    f'-O {fastq_out_dir}'

    os.system(fastq_dump_cmd)

    try:
        shutil.rmtree(f'{sra_accession}/')
    except:
        pass

    return fastq_out_dir


def trim_reads(
        fastq_out_dir: str,
        common_adapters_path: str,
        min_read_len: int,
        min_quality: int = 18
        ) -> None:

    read_files = glob.glob(f'{fastq_out_dir}/*_1.fastq.gz')

    read_info = []

    bbduk_cmd = 'bbduk.sh ' \
                '-Xmx20g ' \
                'ktrim=r ' \
                'k=23 ' \
                'mink=11 ' \
                'hdist=1 ' \
                'qtrim=rl ' \
                f'trimq={min_quality} ' \
                f'ref={common_adapters_path} ' \
                f'minlen={min_read_len} '

    if not read_files:
        print('Unable to find compressed reads!')
        sys.exit(1)

    else:
        for fwd in read_files:
            out_fwd = f'Trimmed_Reads/{fastq_out_dir.rpartition("/")[-1]}.FPE.{min_quality}Q.{min_read_len}nt.fastq.gz'
            out_rev = out_fwd.replace(".FPE.",".RPE.")

            rev = fwd.replace("_1.fastq.gz","_2.fastq.gz")

            if os.path.isfile(rev):
                bbduk_cmd += f'in1={fwd} in2={rev} out1={out_fwd} out2={out_rev}'
                read_info.append((out_fwd, out_rev))

            else:
                bbduk_cmd += f'in={fwd} out={out_fwd.replace(".FPE.",".SE.")}'
                read_info.append(tuple([out_fwd.replace(".FPE.",".SE.")]))

            os.system(bbduk_cmd)

        shutil.rmtree(fastq_out_dir)

        return read_info



def prepare_sra_dataset(
        sra_lines: str,
        common_adapters_path: str,
        min_quality: int = 18,
        threads: int = 8
        ) -> None:

    reads_for_assembly = defaultdict(list)

    trim_lengths = {
        50: 40, 75: 60, 90: 75, 100: 75, 125: 100, 150: 125, 175:150, 200: 150, 225: 175,
        250: 200, 275: 225, 300: 250
        }

    Path('Trimmed_Reads').mkdir(exist_ok = True)

    for line in sra_lines:

        taxon, sra_accession, read_len = line.rstrip().split(",")

        taxon = taxon.replace(" ","_").replace(".","").replace('-','_')

        sra_accession = [sra_accession]

        if ';' in sra_accession[0]:
            sra_accession = [i.strip() for i in sra_accession[0].split(";")]

        min_read_len = trim_lengths[int(read_len)]

        if not (taxon or sra_accession or read_len):

            print('Table may be improperly formatted!')
            sys.exit()

        else:
            print(f'---- {taxon}, {sra_accession}, {read_len} ----')

            for acc in sra_accession:
                download_sra_file(
                    acc
                    )

                fastq_out_dir = sra_fastq_dump(
                                    taxon,
                                    acc,
                                    threads
                                    )


                trimmed_reads = trim_reads(
                                    fastq_out_dir,
                                    common_adapters_path,
                                    min_read_len,
                                    min_quality
                                    )

                reads_for_assembly[taxon] += trimmed_reads

    return reads_for_assembly


def assemble_transcriptome(
        reads_for_assembly: dict,
        threads: int = 8
        ) -> None:

    Path('Assembled_Transcriptomes/').mkdir(exist_ok = True)

    print(reads_for_assembly)

    for k, v in reads_for_assembly.items():
        read_prep = []
        taxon = k.replace(".","").replace('-','_').replace(' ','_')
        n = 1
        for i in v:
            if len(i) > 1:
                read_prep.append(f'--pe-1 {n} {i[0]} --pe-2 {n} {i[1]}')
            else:
                read_prep.append(f'--s {n} {i[0]}')
            n += 1

        read_cmd = ' '.join(read_prep)

        assembly_cmd = f'spades.py --rna -t {threads} -m 48 {read_cmd} -o {taxon}_rnaSPAdes'

        os.system(assembly_cmd)

        os.system(f'rm Trimmed_Reads/{taxon}*fastq.gz')

        shutil.copy(f'{taxon}_rnaSPAdes/transcripts.fasta',f'Assembled_Transcriptomes/{taxon}.rnaSPAdes.fasta')
        shutil.rmtree(f'{taxon}_rnaSPAdes')


if __name__ == '__main__':

    try:
        sra_table = sys.argv[1]

    except:
        print('\nUsage:\n\n    python3 sra_download_prep.py [SRA-ACCESSION-TABLE]\n\n')
        sys.exit()

    assemble = True

    common_adapters_path = 'common_adapters.fasta'

    if not common_adapters_path:
        print('Update the "common_adapters_path" with the path to the "common_adapters.fasta" file')
        print('This is around line 203')
        sys.exit()

    min_quality = 18

    threads = 8

    batch_size = 1

    all_accession_info = [line for line in open(sra_table).readlines()]

    for i in range(0, len(all_accession_info), batch_size):
        sra_lines = all_accession_info[i:i+batch_size]

        reads_for_assembly = prepare_sra_dataset(
                                sra_lines,
                                common_adapters_path,
                                min_quality,
                                threads
                                )

        if assemble:
            assemble_transcriptome(
                reads_for_assembly,
                threads)
