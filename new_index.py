#!/usr/bin/env python3.6

import argparse
import itertools
from collections import Counter


def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    ntds = list(seq)
    ntds = [complement[ntd] for ntd in ntds]
    return ''.join(ntds)


def reverse_complement(s):
        return complement(s[::-1])


def get_arguments():
    parser = argparse.ArgumentParser(prog='index_hopping', description= 'Demultiplexing of paired end sequences by index paring')
    parser.add_argument('-1', '--R1', help='Read1 file path', required=True, type=argparse.FileType('rt', encoding='UTF-8 '))
    parser.add_argument('-2', '--R2', help='Read2 file path', required=True, type=argparse.FileType('rt', encoding='UTF-8 '))
    parser.add_argument('-3', '--R3', help='Read3 file path', required=True, type=argparse.FileType('rt', encoding='UTF-8 '))
    parser.add_argument('-4', '--R4', help='Read4 file path', required=True, type=argparse.FileType('rt', encoding='UTF-8 '))
    parser.add_argument('-c', '--qcutoff', help='Quality score cutoff', required=True, type=int)
    parser.add_argument('-i', '--ind', help='Index file path', required=True, type=argparse.FileType('rt', encoding='UTF-8 '))

    return parser.parse_args()


args = get_arguments()

min_score = args.qcutoff

R1 = args.R1.name
R2 = args.R2.name
R3 = args.R3.name
R4 = args.R4.name
index = args.ind.name


index_dict = {}

with open(index, 'r') as ind:
    for line in ind:
        line = line.strip('\n').split('\t')
        index_dict[line[3]] = line[4]


all_indices = []

for key in index_dict:
    all_indices.append(index_dict[key])


all_pairs_list = list(
    itertools.combinations_with_replacement(all_indices, 2)
    )

all_pairs_dict = {}

for pair in all_pairs_list:
    all_pairs_dict[pair] = 0

undetermined_dict = Counter()

with open(R1 ,'r') as r1, open(R2, 'r') as r2, open(R3, 'r') as r3, open(R4, 'r') as r4:
    NL = 0
    
    for line in zip(r1, r2, r3, r4):
        line = [l.strip() for l in line]
        if NL % 4 == 0:
            head = line
        if NL % 4 == 1:
            seq = line
        if NL % 4 == 2:
            plus = line
        if NL % 4 == 3:
            qline = line
            for j in qline[1:2]:
                qlist = []
                qtotal = 0
                for qual in j:
                    qtotal += ord(qual) - 33
                    qavg = qtotal/len(j)
                    if qavg >= min_score:
                        qlist.append(qual)
                        pair = seq[1], reverse_complement(seq[2])
                        if pair in all_pairs_dict:
                            all_pairs_dict[pair] += 1
                        else:
                            undetermined_dict[pair] += 1
        NL += 1


match = open('matched.tsv', 'w')
match.write('Matched Index Pair\tCounts\n')

swap = open('swapped.tsv', 'w')
swap.write('Swapped Index Pair\tCounts\n')

und = open('undetermined.tsv', 'w')
und.write('Undetermined Index Pair\tCounts\n')


swapped = 0
for pair in all_pairs_dict:
    if pair[0] == pair[1]:
        match.write(str(pair[0]) + '_' + str(pair[1]) + '\t' + str(all_pairs_dict[pair]) + '\n')
    else:
        swapped += 1
        swap.write(str(pair[0]) + '_' + str(pair[1]) + '\t' + str(all_pairs_dict[pair]) + '\n')

for pair in undetermined_dict:
    und.write((str(pair[0]) + '_' + str(pair[1]) + '\t' + str(undetermined_dict[pair]) + '\n'))


match = 0
for pair in all_pairs_dict:
    match += all_pairs_dict[pair]

undetermined = 0
for com in undetermined_dict:
    undetermined += undetermined_dict[com]

stats = open('metrics.tsv', 'w')
stats.write('Filtered' + '\t' + 'Matched' + '\t' + 'Swapped' + '\t' + 'Undetermined' + '\n')
stats.write(str(match+swapped+undetermined) + '\t' + str(match) + '\t' + str(swapped) + '\t' + str(undetermined)+'\n')
