#!/usr/env/py27
# Visualize C>A portions of mutational spectra in Stratton format
# Runs in py36 too

from figures import spectrum_map
from collections import OrderedDict
import cospec as sa
import matplotlib.pyplot as plt
import scipy.stats as sc
import scipy.cluster.hierarchy as hac
from sklearn.metrics.pairwise import cosine_similarity

###----- PARAMETERS TO EDIT -----
cosmic_txt = "data/cosmic_signature_probabilities.txt"
norm_8oxoG_all = "data/8oxoG-mutations-all_normalized.txt"
context_freq_txt = "data/tko/tko-mice-context.txt"
mut_txt = "data/tko/tko-mice-mut.txt"

#sigs = ['signature_4', 'signature_18', 'signature_24', 'signature_29']
###----- END EDIT -----
mut_types = ['C>A','C>G','C>T','T>A','T>C','T>G']

#cosmic_sig = sa.refine(sa.take_CtA(sa.convert_cosmic_sigs(cosmic_txt)))
#muty_spec = sa.refine(sa.take_CtA(sa.convert_mutations(norm_8oxoG_all)))
muty_spec = sa.refine(sa.convert_mutations(norm_8oxoG_all))
cosmic_sigs = sa.refine(sa.convert_cosmic_sigs(cosmic_txt))
sigs = list(cosmic_sigs.keys())

compare_spec = OrderedDict()
compare_spec['mutY'] = muty_spec
for sig in sigs:
    compare_spec[sig] = cosmic_sigs[sig]

con_freq = sa.get_con_freq(context_freq_txt)
tko_spec = sa.refine(sa.convert_tko_muts(mut_txt, con_freq))

spectra = [list(muty_spec.values())] + [list(tko_spec.values())]
for sig in list(cosmic_sigs.keys()):
    spectra.append(list(cosmic_sigs[sig].values()))
cos_list = cosine_similarity(spectra)[0]
sig_names = ['mutY','TKO'] + [sig for sig in list(cosmic_sigs.keys())]
for i in range(len(cos_list)):
    print(sig_names[i] + ": " + str(cos_list[i]))

#sa.plot_uhc_dendrogram(cosmic_sigs, tko_spec, ['TKO']+list(cosmic_sigs.keys()))
#sa.plot_uhc_corr_dendrogram(cosmic_sigs, tko_spec, ['TKO']+list(cosmic_sigs.keys()))
#sa.plot_uhc_corr_heatmap(cosmic_sigs, muty_spec, ['mutY']+list(cosmic_sigs.keys()), False)
#sa.plot_uhc_corr_dendrogram(sa.take_CtA(cosmic_sigs), sa.take_CtA(muty_spec), ['mutY']+list(cosmic_sigs.keys()))
#sa.plot_uhc_corr_heatmap(sa.take_CtA(cosmic_sigs), sa.take_CtA(muty_spec), ['mutY']+list(cosmic_sigs.keys()), False)


'''
# Calculate Spearman rank-order correlation coefficient and p-value
# Correlations of -1 or +1 imply exact monotonic relationships.
# P-values are "not entirely reliable"
spectra = [list(compare_spec[sig].values()) for sig in list(compare_spec.keys())]
#print(sc.spearmanr(spectra)[0])
#print(sc.spearmanr(spectra)[0][0])
#print(sc.spearmanr(spectra[0],spectra[1]))
arr = [None]*len(spectra)
for i in range(0,31):
    arr[i] = [sc.spearmanr(spectra[i],spectra[j])[0] for j in range(0,31)]

#for i in range(0,31):
#    print(str(i) + ": " + str(arr[0][i]))

# Use 'single' as UHC method
corr = hac.linkage(arr, 'single')
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
#plt.show()
den_file = '8oxoG-cosmic_corr_dendrogram.pdf'
plt.savefig(den_file, dpi=450)
'''

'''
con_freq = sa.get_con_freq(context_freq_txt)
tko_spec = sa.refine(sa.convert_tko_muts(mut_txt, con_freq))
title = "Mutational Spectrum for 8oxoG TKO"
spectrum_map(x=1, y=1,
             heights=[[x*100 for x in tko_spec.values()]],
             xlabels=[[i[1] for i in tko_spec.keys()]],
             labels=mut_types,
             titles=[title])
plt.savefig('tko_spec_v03.pdf')
'''

'''
title = "Mutational Spectrum for mutY- KO"
spectrum_map(x=1, y=1,
             heights=[[x*100 for x in muty_spec.values()]],
             xlabels=[[i[1] for i in muty_spec.keys()]],
             labels=mut_types,
             titles=[title])
plt.savefig('mutY_spec_vXY.pdf')
'''

'''
for sig in compare_spec.keys():
    title = "Spectrum for " + str(sig)
    spectrum_map(x=1, y=1,
                 heights=[[x*100 for x in compare_spec[sig].values()]],
                 xlabels=[[i[1] for i in compare_spec[sig].keys()]],
                 labels=mut_types,
                 titles=[title])
    plt.savefig('mut-spec_'+str(sig)+'.pdf')
'''
