#!/usr/env/py36
#
# -----
# NGS MiSeq Analysis Pipeline
# Visual analysis and comparison of mutation data
# Version 1.4
# Author: Lina Kim, klkim [at] mit.edu
# May 7, 2017
# -----
#
# Functions for analysis of COSMIC signatures and additional spectra
# Reads in text file of COSMIC mutational signatures
# Imports and stores in OrderedDict, the preferred method for signatures
#

from figures import spectrum_map
from collections import OrderedDict
from scipy.cluster.hierarchy import dendrogram
from fastcluster import linkage
from sklearn.metrics.pairwise import cosine_similarity

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hac
import scipy.stats as sc

def refine(mut_sig): # ensures OrderedDicts are in same order w.r.t keys
    # rewritten to allow for COSMIC dict input
    # cleans up dict so there's no impossibilities (ex: ACA, T>A)

    mutations = ['C>A','C>G','C>T','T>A','T>C','T>G']
    contexts = ['ACA','ACC','ACG','ACT',
                'CCA','CCC','CCG','CCT',
                'GCA','GCC','GCG','GCT',
                'TCA','TCC','TCG','TCT',
                'ATA','ATC','ATG','ATT',
                'CTA','CTC','CTG','CTT',
                'GTA','GTC','GTG','GTT',
                'TTA','TTC','TTG','TTT']
    new_sig = OrderedDict()
    if isinstance(list(mut_sig.keys())[0],str):
        for sig in list(mut_sig.keys()):
            inner_sig = OrderedDict()
            new_sig[sig] = inner_sig
            for m in mutations:
                for c in contexts:
                    if m[0] == c[1]:
                        inner_sig[(m,c)] = mut_sig[sig].setdefault((m,c),0)
        return new_sig
    for m in mutations:
        for c in contexts:
            if m[0] == c[1]:
                new_sig[(m,c)] = mut_sig.setdefault((m,c),0)
    return new_sig

def convert_cosmic_sigs(in_file):
    spect_dict = OrderedDict()
    firstLine = True
    f = open(in_file, 'r')
    for line in f:
        line = line.strip().split('\t')
        if firstLine:
            for num in line[3:]: # for each signature, past first three columns
                sig = "signature_" + str(num[10:]) # signature number
                spect_dict[sig] = OrderedDict()
            firstLine = False
            continue
        for col in range(1,31):
            sig = "signature_" + str(col)
            spect_dict[sig][(line[0],line[1])] = float(line[col+2]) # key is (mutation, context)
    return refine(spect_dict)

def convert_mutations(in_file): # for all mutations, normalized
    spect_dict = OrderedDict()
    mut_opt = []
    firstLine = True
    f = open(in_file, 'r')
    for line in f:
        line = line.strip().split(',')
        if firstLine:
            for mut in line: # each mutation, including empty first to align with index
                mut_opt.append(mut)
            firstLine = False
            continue
        for col in range(1,len(line)):
            spect_dict[(mut_opt[col],line[0])] = float(line[col]) # key is (mutation, context)
    return refine(spect_dict)

def take_CtA(spectrum): # make accessible to COSMIC nested dictionaries too
    if isinstance(list(spectrum.keys())[0],str):
        copy_dict = OrderedDict()
        for sig in list(spectrum.keys()):
            copy = spectrum[sig].copy()
            for mut, con in list(copy.keys()):
                if mut != 'C>A':
                    del copy[(mut,con)]
            copy_dict[sig] = copy
        return copy_dict
    copy = spectrum.copy()
    for mut, con in list(copy.keys()):
        if mut != 'C>A':
            del copy[(mut,con)]
    return copy

def uhc_cluster(cosmic_list, ref_sig):
    spectra = [list(ref_sig.values())] # so ref signature is value 0
    for sig in cosmic_list:
        spectra.append(list(cosmic_list[sig].values()))
    return linkage(spectra, method='ward', metric='cosine')

def plot_uhc_heatmap(cosmic_list, ref_sig, cluster_names, isText):
    # clean up to figure axes properly
    # also to account for different types of inputs
    # also to account for different parameters: annotations, methods, output format, etc.
    spectra = [list(ref_sig.values())]
    for sig in cosmic_list:
        spectra.append(list(cosmic_list[sig].values()))
    linkages = uhc_cluster(cosmic_list, ref_sig)
    with plt.rc_context({'lines.linewidth': 1.25}):
        columns = cluster_names
        grid = sns.clustermap(pd.DataFrame(cosine_similarity(spectra), [columns, ['']*len(columns)]),
                            method='weighted',
                            row_cluster=True, col_cluster=True,
                            row_linkage=linkages, col_linkage=linkages,
                            annot=isText, fmt='.2f', # change this line to show values in heatmap
                            cmap=sns.cubehelix_palette(start=9, rot=-0.2, dark=0.15, light=.85, as_cmap=True))

    _ = grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), size=9)
    _ = grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), size=9, rotation=0)
    plt.show()
#    plt.savefig(heatmap_file, dpi=450, bbox_inches='tight')

def plot_uhc_dendrogram(cosmic_list, ref_sig, cluster_names):
    spectra = [list(ref_sig.values())]
    for sig in cosmic_list:
        spectra.append(list(cosmic_list[sig].values()))
    data = spectra
    labels = cluster_names
    linkages = uhc_cluster(cosmic_list, ref_sig)
    with plt.rc_context({'lines.linewidth':0.5}):
        dendrogram(linkage(pd.DataFrame(data, [labels, ['']*len(labels)]),
                                 method='ward',
                                 metric='cosine'),
                         no_labels=False,
                         labels=labels,
                         leaf_rotation=90,
                         leaf_font_size=4,
                         link_color_func=lambda x:'k')
    plt.show()
#    plt.savefig('file.pdf', dpi=450)

def plot_uhc_corr_heatmap(cosmic_list, ref_sig, cluster_names, isText):
    # clean up to figure axes properly
    # also to account for different types of inputs
    # also to account for different parameters: annotations, methods, output format, etc.
    spectra = [list(ref_sig.values())]
    for sig in cosmic_list:
        spectra.append(list(cosmic_list[sig].values()))

    arr = [None]*len(spectra)
    for i in range(0,31):
        arr[i] = [sc.spearmanr(spectra[i],spectra[j])[0] for j in range(0,31)]

    # Use 'single' as UHC method
    linkages = linkage(arr, 'single')
    with plt.rc_context({'lines.linewidth': 1.25}):
        columns = cluster_names
        grid = sns.clustermap(pd.DataFrame(arr, [columns, ['']*len(columns)]),
                            method='single',
                            row_cluster=True, col_cluster=True,
                            row_linkage=linkages, col_linkage=linkages,
                            annot=isText, fmt='.2f',
                            cmap=sns.cubehelix_palette(start=9, rot=-0.2, dark=0.15, light=.85, as_cmap=True))

    _ = grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), size=9)
    _ = grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), size=9, rotation=0)
    plt.show()
    heatmap_file = 'CtA_8oxoG-cosmic_corr_heatmap.pdf'
#    plt.savefig(heatmap_file, dpi=450, bbox_inches='tight')

def plot_uhc_corr_dendrogram(cosmic_list, ref_sig, cluster_names):
    spectra = [list(ref_sig.values())]
    for sig in cosmic_list:
        spectra.append(list(cosmic_list[sig].values()))

    arr = [None]*len(spectra)
    for i in range(0,31):
        arr[i] = [sc.spearmanr(spectra[i],spectra[j])[0] for j in range(0,31)]

    # Use 'single' as UHC method
    corr = hac.linkage(arr, 'ward')
    plt.figure(figsize=(25,10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Index of Signatures')
    plt.ylabel('Distance')
    with plt.rc_context({'lines.linewidth':0.5}):
        hac.dendrogram(
            corr,
            leaf_rotation=90,
            leaf_font_size=10,
            link_color_func=lambda x:'k'
            )
    plt.show()
    den_file = 'ward_TKO-cosmic_corr_dendrogram.pdf'
#    plt.savefig(den_file, dpi=450)

# First column is context
# Second column is frequency
def get_con_freq(in_file):
    con_freq = {}
    f = open(in_file, 'r')
    for line in f:
        line = line.strip().split(',')
        con_freq[line[0]] = line[1]
    return con_freq

# First column is mutation, ex: 'C>A'
# Second is mouse origin
# Third is context, 50 [C/A] 50
# Output OrderedDict of mutational spectrum
def convert_tko_muts(in_file, con_freq):
    base_pair = {'A':'T','C':'G','G':'C','T':'A'}
    pyrimidines = ['C','T']

    spect_dict = OrderedDict()
    firstLine = True
    f = open(in_file, 'r')
    for line in f:
        mut = None
        con = None
        line = line.strip().split(',')
        if firstLine:
            firstLine = False
            continue
        if line[0][0] in pyrimidines:
            mut = line[0]
            con = line[2][49] + line[2][51] + line[2][55]
        else:
            mut = base_pair[line[0][0]]+">"+base_pair[line[0][2]]
            con = base_pair[line[2][55]] + base_pair[line[2][51]] + base_pair[line[2][49]]
        spect_dict[(mut,con)] = spect_dict.setdefault((mut,con),0) + 1

    # Find probabilities and normalize
    for mut, con in list(spect_dict.keys()):
        spect_dict[(mut,con)] = spect_dict[(mut,con)] / float(con_freq[con])
    total = sum(list(spect_dict.values()))
    for mut, con in list(spect_dict.keys()):
        spect_dict[(mut,con)] = spect_dict[(mut,con)] / total
    return refine(spect_dict)

def draw_spectrum(spec, title, fileName): # OrderedDict of spectrum/signature
    mut_types = ['C>A','C>G','C>T','T>A','T>C','T>G']
    spectrum_map(x=1, y=1,
                 heights=[[x*100 for x in list(spec.values())]],
                 xlabels=[[i[1] for i in list(spec.keys())]],
                 labels=mut_types,
                 titles=[title])
    plt.show()
#    plt.savefig(fileName)
