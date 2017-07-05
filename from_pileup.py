#!/usr/env/py36
#
# -----
# NGS MiSeq Analysis Pipeline
# Analysis of processed sequencing data
# Version 1.0
# Author: Lina Kim, klkim [at] mit.edu
# July 5, 2017
# -----
#
# Take pileup mapping data and convert to Excel
# Make information easily accessible
# Cell, barcode, lesion, position (6244-6280), refnt, coverage,
#   A%, G%, C%, T%, N%, D%, I%, *%,
# NEW ADDITION: allow for multiple input pileup files
#

import argparse
import re

def start(output):
    with open(output, 'w') as fo:
        info = ['cell','run','barcode','lesion','position','refnt','coverage']
        muts = ['A','G','C','T','N','-','+','*']
        fo.write(','.join(info + muts))
        fo.write('\n')

def extract(inputs, output, sample_codes):
    for input in inputs:
        muts = ['A','G','C','T','N','-','+','*']
        with open(input,'r') as fi, open(output,'a') as fo:
            nextLines = []
            while True:
                line = fi.readline().split("\t")
                if line == [""]: break

                sample = str(input)[input.rfind('/')+1:-15]
                cell = sample_codes[sample][:-2]
                run = sample_codes[sample][-1]
                barcode = str(input)[-14:-11]
                lesion = str(input)[-10:-7]
                pos = str(line[1])
                refnt = str(line[2])
                coverage = str(line[3])

                readbase = str(line[4])

                # [,.] are matches
                # [*] is del
                # [ACGTNacgtn] are mismatches
                # remove number and nt following +/- and record for insertions and deletions
                # remove char following ^ as quality score
                eor = r'(\^)[\x00-\x7F]' # just remove/ignore; ASCII char
                eol = r'\$' # just remove/ignore
                indel = r'[+-][0-9]{1,}' # greedy match, for most possible matches in case indel is more than 1 digit number long
                # starting from the NEXT position...? so remove and record for next few positions?
                delnt = r'\*'
                match = r'[\,\.]'
                mismt = r'[ACGTNacgtn]'

                percentages = {}
                for key in muts:
                    percentages[key] = 0

                num_indel = [(x.start(),x.end()) for x in re.finditer(indel,readbase)] # indices of beginnings
                to_remove = []
                for s,e in num_indel:
                    nt_counts = int(readbase[s+1:e]) # number of nt affected by indel, starting from the character after +/-
                    to_remove.append((s+1,e+nt_counts)) # to account for nt affected by indel

                new_bases = ''
                last_pos = 0
                for s,e in to_remove:
                    new_bases = new_bases + readbase[last_pos:s]
                    last_pos = e
                new_bases = new_bases + readbase[last_pos:]
                new_bases.upper()
                for nt in new_bases:
                    if (nt == '.') | (nt == ','):
                        percentages[refnt] = percentages[refnt] + 1
                    elif nt in muts:
                        percentages[nt] = percentages[nt] + 1

                total = sum([percentages[i] for i in ['A','C','G','T','N','*']])
                print(pos,readbase,[percentages[key] for key in muts])
                for key in muts:
                    percentages[key] = percentages[key] / float(total)

                l_out = ','.join([cell, run, barcode, lesion, pos, refnt, coverage] + [str(percentages[i]) for i in muts])

                # Write to file
                fo.write(l_out+'\n')

def main():
    parser = argparse.ArgumentParser(
        description = ("Takes in pileup files with mutation data and  "
                        "extracts necessary information regarding bases, "
                        "coverage, and percentage of mutations.")
    )
    parser.add_argument("-i", "--inputs",
                        action = "store",
                        type = str,
                        dest = "inputs",
                        help = "File path to pileup files with mutation data.",
                        required = True)
    parser.add_argument("-o", "--output",
                        action = "store",
                        type = str,
                        dest = "output",
                        help = "Destination file for output.",
                        required = True)

    args = parser.parse_args()
    inp = args.inputs.split(',') # is a list of input files
    out = args.output

    scs = dict([('D15-10283','AB1157-1'), ('D15-10284','AB1157-2'), ('D15-10285','AB1157-3'),
                ('D15-10286','mutY-1'), ('D15-10287','mutY-2'), ('D15-10288','mutY-3'),
                ('D15-10289','mutYSOS-1'), ('D15-10290','mutYSOS-2'), ('D15-10291','mutYSOS-3'),])

    start(out)
    extract(inp, out, scs)

if __name__ == '__main__':
    main()
