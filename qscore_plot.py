#!/usr/bin/env python3.6

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import argparse


def get_arguments():
    parser = argparse.ArgumentParser(prog='qscore_plot', description= 'Outputs a distribution of quality scores as a function of the base position of sequencing reads.')

    parser.add_argument('-f', '--infile', help='fastq filepath', required=True, type=argparse.FileType('rt', encoding='UTF-8 '))

    return parser.parse_args()

def convert_phred(letter):
    score = ord(letter) - 33    # phred + 33
    return score


args = get_arguments()
file = args.infile.name

all_scores = np.zeros(101)
mean_read_scores = Counter()


with open(file, 'r') as fh:
    count = 0
    NR = 1
    for line in fh:
        line = str(line).strip('\n')
        count = 2
        if NR % 4 == 0:  # quality line
            ntd_pos = 0
            read_score = 0
            for ntd in line:
                all_scores[ntd_pos] = ((int(convert_phred(ntd)) + all_scores[ntd_pos]) / (count))

                read_score += int(convert_phred(ntd))
                read_avg = int(read_score // len(line))  # floor div
                ntd_pos += 1
            count += 1

            mean_read_scores[read_avg] += 1
        NR += 1 

# quality score distribution of reads
xdata = np.arange(0, len(all_scores[all_scores > 0]), 1)
plt.figure()
plt.bar(xdata, all_scores[all_scores > 0], width = 0.5)

print(all_scores[all_scores > 0])

plt.title(str(file) + '\n' + 'Mean Quality Score by Base Position')
plt.xlabel('Base Position')
plt.ylabel('Mean Quality Score')
plt.savefig("/home/kkinning/hop/plots/" + file + '_distribReads.png')



# quality score distribution across all reads
x = []
y = []
for key in mean_read_scores:
    x.append(key)
    y.append(mean_read_scores[key])


plt.figure()
plt.bar(x, y, width = 0.5)
plt.title(str(file) + '\n' + 'Frequency of Average Quality Score Across All Reads')
plt.xlabel('Avg score per read')
plt.ylabel('Frequency')
plt.savefig("/home/kkinning/hop/plots/" + file + '_distribAll.png')
