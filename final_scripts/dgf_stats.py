import numpy as np
import matplotlib.pyplot as plt
import time
import math
import json
import datetime
import argparse
from fasta import FASTA

import argparse

# Parse arguments
parser = argparse.ArgumentParser(
    description="Compresses and Decompresses number sequence files")

parser.add_argument("simulation", help="Simulation Folder")

parser.add_argument("-c", "--count", type=int, help="Number of simulations")
parser.add_argument("-C", "--chromosomes", type=int,
                    help="Number of chromosomes")
parser.add_argument("-f", "--fasta", type=str,
                    help="Fasta Annotated Transcripts File")
parser.add_argument("-s", "--search",
                    type=str, help="Search String")
parser.add_argument("-b", "--basepairs", action='store_true',
                    help="Use base pairs instead of genome counts")

args = parser.parse_args()

protein_fasta = args.fasta or '/home/seijihariki/Documents/TCC/TTDB/TriTrypDB-46_TcruziCLBrenerEsmeraldo-like_AnnotatedTranscripts.fasta'
simulation_folder = args.simulation
search = args.search or 'DGF-1'
base_pairs = args.basepairs

simulation_cnt = args.count or 50
chromosomes_cnt = args.chromosomes or 41


print('Loading annotations:')
transcripts = FASTA(protein_fasta)
transcripts.load()

collisions = {}

print('Detecting collisions:')
for chromosome in range(chromosomes_cnt):
    chromosome_name = f"TcChr{chromosome + 1}-S"
    collisions[chromosome_name] = []

    for simulation in range(simulation_cnt):
        with open(f"{simulation_folder}simulation_{simulation}/{chromosome_name}.cseq") as times:
            start, end = -2, -2
            current_location = 0

            for entry in times:
                xsplit = entry.split('x')
                dsplit = xsplit[0].split('-')

                amount = int(xsplit[1]) if len(xsplit) == 2 else 1
                nstart, nend = int(dsplit[0]), int(
                    dsplit[1] if len(dsplit) == 2 else dsplit[0])

                if end >= 0 and abs(end - nstart) < 2:
                    collisions[chromosome_name].append(current_location)

                current_location += amount * (abs(nend - nstart) + 1)

                start = nstart
                end = nend

statistics = {}

genes_regions = {}
dgf_1_regions = {}

dgf_collisions = {}
other_collisions = {}

print('Calculating statistics:')
for key in collisions.keys():
    collision_locs = collisions[key]
    dgf_1_regions[key] = []
    genes_regions[key] = []
    dgf_collisions[key] = []
    other_collisions[key] = []

    statistics[key] = {
        'in_dgf': 0,
        'out_dgf': 0,
        'dgf_cnt': 0,
        'other_cnt': 0,
        'dgf_total_bp': 0,
        'other_total_bp': 0
    }

    for context in transcripts.data.keys():
        try:
            params = transcripts.data[context]['params']
            chromosome = params['location'].split(':')[0]

            if params['gene_product'] and params['gene_product'].find(search) >= 0 and chromosome == key:
                statistics[key]['dgf_cnt'] += 1
                ch_range = list(map(int, params['location'].split(':')[
                    1].split('(')[0].split('-')))

                start, end = min(ch_range), max(ch_range)

                dgf_1_regions[key].append((start, end))

                statistics[key]['dgf_total_bp'] += end - start

                statistics[key]['in_dgf'] += len(
                    list(filter(lambda x: start < x and x < end, collision_locs)))
                dgf_collisions[key] += list(filter(lambda x: start <
                                                   x and x < end, collision_locs))

            elif params['location'].split(':')[
                    1].split('(')[0].split('-')[0] != '' and chromosome == key:
                statistics[key]['other_cnt'] += 1
                ch_range = list(map(int, params['location'].split(':')[
                    1].split('(')[0].split('-')))

                start, end = min(ch_range), max(ch_range)

                genes_regions[key].append((start, end))

                statistics[key]['other_total_bp'] += end - start
                other_collisions[key] += list(filter(lambda x: start < x and x < end,
                                                     collision_locs))

        except Exception as e:
            print('Error: ', e)
            pass

    statistics[key]['out_dgf'] = len(
        collisions[key]) - statistics[key]['in_dgf']

print()
print('Collisions:')


def bp_normalized():
    global_dgf = 0
    global_dgf_bp = 0

    global_other = 0
    global_other_bp = 0

    for key in collisions.keys():
        stats = statistics[key]

        normalized_dgf = 0
        normalized_other = 0

        global_dgf += stats['in_dgf']
        global_dgf_bp += stats['dgf_total_bp']

        global_other += stats['out_dgf']
        global_other_bp += stats['other_total_bp']

        if stats['dgf_total_bp'] > 0:
            normalized_dgf = stats['in_dgf'] / stats['dgf_total_bp']

        normalized_other = stats['out_dgf'] / stats['other_total_bp']

        print(f"{key}:".ljust(15) + f"Base Pairs: {stats['dgf_total_bp']}/{stats['other_total_bp']}".ljust(
            30) + f"DGF/OTHER={normalized_dgf/normalized_other}".ljust(30))
        print(stats['in_dgf'] + stats['out_dgf'])

    global_normalized_dgf = global_dgf / global_dgf_bp
    global_normalized_other = global_other / global_other_bp

    print(
        f"GLOBAL: DGF/OTHER={global_normalized_dgf / global_normalized_other}")


def gc_normalized():
    global_dgf = 0
    global_dgf_gc = 0

    global_other = 0
    global_other_gc = 0

    for key in collisions.keys():
        stats = statistics[key]

        normalized_dgf = 0
        normalized_other = 0

        global_dgf += stats['in_dgf']
        global_dgf_gc += stats['dgf_cnt']

        global_other += stats['out_dgf']
        global_other_gc += stats['other_cnt']

        if stats['dgf_cnt'] > 0:
            normalized_dgf = stats['in_dgf'] / stats['dgf_cnt']

        normalized_other = stats['out_dgf'] / stats['other_cnt']

        print(f"{key}:".ljust(15) + f"Gene Counts: {stats['dgf_cnt']}/{stats['other_cnt']}".ljust(
            30) + f"DGF/OTHER={normalized_dgf/normalized_other}".ljust(30))
        print(stats['in_dgf'] + stats['out_dgf'])

    global_normalized_dgf = global_dgf / global_dgf_gc
    global_normalized_other = global_other / global_other_gc

    print(
        f"GLOBAL: DGF/OTHER={global_normalized_dgf / global_normalized_other}")


if base_pairs:
    bp_normalized()
else:
    gc_normalized()

global_dgf_gc = 0
global_other_gc = 0

global_dgf_bp = 0
global_other_bp = 0

for key in collisions.keys():
    stats = statistics[key]

    global_dgf_gc += stats['dgf_cnt']
    global_other_gc += stats['other_cnt']

    global_dgf_bp += stats['dgf_total_bp']
    global_other_bp += stats['other_total_bp']

print(
    f"AVERAGE GENE LENGTH: DGF/OTHER={global_dgf_bp / global_dgf_gc}/{global_other_bp / global_other_gc}")
