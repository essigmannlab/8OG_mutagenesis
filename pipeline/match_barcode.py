#!/usr/env/py36
#
# -----
# NGS MiSeq Analysis Pipeline
# Alignment and selection-matching of relevant reads
# Version 1.0
# Author: Lina Kim, klkim [at] mit.edu
# July 5, 2017
# -----
#
# Select reads matching 13-mers around barcode and +-5nt around lesion (6235-6247)
#   Replace low-quality bases with N
#   Extract downstream sequences
#   Split reads based on barcode, then lesion context
#

import re
import argparse

def reverse(read):
    comp = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    backwards = ''
    for nt in read[::-1]:
        backwards = backwards + comp[nt]
    return backwards

def extract(inputs, runid, sam, bcs):
    refs = ['AXA', 'AXC', 'AXG', 'AXT',
            'CXA', 'CXC', 'CXG', 'CXT',
            'GXA', 'GXC', 'GXG', 'GXT',
            'TXA', 'TXC', 'TXG', 'TXT']
    barcodes = bcs
    out_files = {}
    out_id = runid
    out_sam = sam
    for bc in barcodes:
        for r in refs:
            out_files[(bc,r)] = out_id+'/03/'+sam+'/'+bc+'/'+sam+'_'+bc+'_'+r+'.matched.fastq'
    match_mer = r'ACGGT([ACTGN]{3})TGCTC([ACTGN]{6})AGACC([ACTGN]{3})GCGTC' # stringent approach
    lowQ = ['!','\"','#','$','%','&','\'','(',')','*','+'] # Q <= 10

    for input in inputs:
        failed = 0
        short = 0
        with open(input,'r') as fi:
            while True:
                line1 = fi.readline().strip() # @ sequence ID
                if line1 == "": break
                line2 = fi.readline().strip() # sequence
                line3 = fi.readline().strip() # + sequence ID
                line4 = fi.readline().strip() # quality score

                # Replace low-quality bases with N
                for i in range(len(line4)):
                    if line4[i] in lowQ:
                        line2 = line2[:i] + 'N' + line2[i+1:]

                # Select reads with perfect match 13mers around lesion barcode
                matchF = re.search(match_mer, line2)
                if matchF:
                    # Only extract downstream sequences of barcode, inclusive
                    start = matchF.start(0)

                    # Write to different files based on lesion barcode and lesion context
                    barcode = line2[start+5:start+8]
                    if barcode in barcodes:
                        if len(line2) < (start+26+1): # discard reads that are too short to match barcode but exclude lesion
                            short = short + 1
                            continue
                        else:
                            con = line2[start+24]+'X'+line2[start+26]
                            if con in refs: # to ensure no N in lesion context
                                with open(out_files[(barcode,con)], 'a') as fo:
                                    fo.write(line1 + '_F' + '\n')
                                    fo.write(line2[start+8:] + '\n') # excluding lesion barcode
                                    fo.write(line3 + '\n')
                                    fo.write(line4[start+8:] + '\n')
                                    continue
                    with open(out_id+'/03/'+sam+'/'+sam+'_NNN.matched.fastq', 'a') as fo: # if dif barcode or N in lesion context
                            fo.write(line1 + '_F'+'\n')
                            fo.write(line2[start+8:] + '\n') # excluding lesion barcode
                            fo.write(line3 + '\n')
                            fo.write(line4[start+8:] + '\n')

                match_merR = r'GACGC([ACTGN]{3})GGTCT([ACTGN]{6})GAGCA([ACTGN]{3})ACCGT'
                matchR = re.search(match_merR, line2)
                if matchR:
                    line2R = reverse(line2)
                    line4R = line4[::-1] # since quality scores need only be read in reverse order

                    start = re.search(match_mer, line2R).start(0)

                    # Replace low-quality bases with N
                    for i in range(len(line4)):
                        if line4R[i] in lowQ:
                            line2R = line2R[:i] + 'N' + line2R[i+1:]
                    # Only extract downstream sequences of barcode, inclusive

                    # Write to different files based on lesion barcode and lesion context
                    barcode = line2R[start+5:start+8]
                    if barcode in barcodes:
                        if len(line2R) < (start+26+1): # discard reads that are too short to match barcode but exclude lesion
                            short = short + 1
                            continue
                        else:
                            con = line2R[start+24]+'X'+line2R[start+26]
                            if con in refs: # to ensure no N in lesion context
                                with open(out_files[(barcode,con)], 'a') as fo:
                                    fo.write(line1 + '_R' + '\n')
                                    fo.write(line2R[start+8:] + '\n') # excluding lesion barcode
                                    fo.write(line3 + '\n')
                                    fo.write(line4R[start+8:] + '\n')
                                    continue
                    with open(out_id+'/03/'+sam+'/'+sam+'_NNN.matched.fastq', 'a') as fo: # if dif barcode or N in lesion context
                            fo.write(line1 + '_R'+'\n')
                            fo.write(line2R[start+8:] + '\n') # excluding lesion barcode
                            fo.write(line3 + '\n')
                            fo.write(line4R[start+8:] + '\n')

                else:
                    failed = failed + 1

        print("Number of reads that failed match filter for " + str(input) + ": " + str(failed))
        print("Number of reads that failed length filter for " + str(input) + ": " + str(short))

def main():
    parser = argparse.ArgumentParser(
        description = ("Performs an match extract for provided files "
                        "by extracting the 13-mers surrounding the lesion barcode "
                        "and replacing low-quality bases with N.")
    )
    parser.add_argument("-i", "--inputs",
                        action = "store",
                        type = str,
                        dest = "input",
                        help = "File paths to input samples.",
                        required = True)
    parser.add_argument("-r", "--run-id",
                        action = "store",
                        type = str,
                        dest = "runid",
                        help = "ID number for general sample run.",
                        required = True)
    parser.add_argument("-s", "--sample",
                        action = "store",
                        type = str,
                        dest = "sample",
                        help = "ID number for sample of run.",
                        required = True)
    parser.add_argument("-b", "--barcodes",
                        action = "store",
                        type = str,
                        dest = "barcodes",
                        help = "List of lesion barcodes.",
                        required = True)

    args = parser.parse_args()

    inp = args.input.split(',')
    r = args.runid
    s = args.sample
    bc = args.barcodes.split(',')
    extract(inp, r, s, bc)

if __name__ == '__main__':
    main()
