from fasta import FASTA
import os
import random as rd
import re
import sqlite3
from parse import parse

MFASeq_folder = 'ReDyMo-CPP/data'
SQLite_DB = 'ReDyMo-CPP/data/database.sqlite'

probabilities = [25, 50, 75, 100]

emeraldo_like = {
    'genome_filename': 'TTDB/TriTrypDB-46_TcruziCLBrenerEsmeraldo-like_Genome.fasta',
    'regions_filename': 'TTDB/TriTrypDB-46_TcruziCLBrenerEsmeraldo-like_AnnotatedCDSs.fasta',
    'organism': 'TcruziCLBrenerEsmeraldo-like'
}
non_emeraldo = {
    'genome_filename': 'TTDB/TriTrypDB-46_TcruziCLBrenerNon-Esmeraldo-like_Genome.fasta',
    'regions_filename': 'TTDB/TriTrypDB-46_TcruziCLBrenerNon-Esmeraldo-like_AnnotatedCDSs.fasta',
    'organism': 'TcruziCLBrenerNon-Esmeraldo-like'
}

organism = emeraldo_like

if __name__ == "__main__":
    # Load FASTA files
    genome = FASTA(organism['genome_filename'])
    genome.load()

    regions = FASTA(organism['regions_filename'])
    regions.load()

    # Load database file
    sqlite = sqlite3.connect(SQLite_DB)

    # Create MFASeq Folder
    Organism_MFASeq_folder = f"{MFASeq_folder}/MFA-Seq_{organism['organism']}"
    if not os.path.isdir(Organism_MFASeq_folder):
        os.mkdir(Organism_MFASeq_folder)

    # Create MFASeq Files
    for chromosome_id in genome.data.keys():
        Chromosome_file = f"{Organism_MFASeq_folder}/{chromosome_id}.txt"
        Chromosome_Prob_file = f"{Organism_MFASeq_folder}/{chromosome_id}_probability.txt"

        if not os.path.isfile(Chromosome_file):
            print(f'Creating MFASeq Scan for {chromosome_id}')
            mfaseq_f = open(
                Chromosome_file, 'w')

            for _ in range(1000):
                mfaseq_f.write(f'{rd.random()}\n')
            mfaseq_f.close()

        if not os.path.isfile(Chromosome_Prob_file):
            print(f'Creating MFASeq Probabilities for {chromosome_id}')
            mfaseq_p = open(
                Chromosome_Prob_file, 'w')

            for _ in range(1000):
                mfaseq_p.write('0.25\n')

            mfaseq_p.close()

    # Insert data into DB
    for chromosome_id in genome.data.keys():
        length = int(genome.data[chromosome_id]['params']['length'])

        sqlite.execute(
            f"INSERT OR REPLACE INTO Chromosome VALUES (\"{chromosome_id}\", {length}, 0.55, 62, 40, 9e999, 0, \"{organism['organism']}\")")

    for region_id in regions.data.keys():
        region = regions.data[region_id]
        location = region['params']['location']
        chromosome_id, start_loc, end_loc, sign = parse(
            '{}:{}-{}({})', location)

        start, end = 0, 0
        if sign == '+':
            start = start_loc
            end = end_loc
        if sign == '-':
            start = end_loc
            end = start_loc

        sqlite.execute(
            f"INSERT OR REPLACE INTO TranscriptionRegion VALUES ({start}, {end}, \"{chromosome_id}\")")

    sqlite.commit()
    sqlite.close()
