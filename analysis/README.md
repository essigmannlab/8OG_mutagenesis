## 8oxoG Analysis
### Version:
1.4

### Purpose:
- Import and convert COSMIC signatures for processing
- Import and convert 8-oxoG mutational spectrum for processing
- Match cosine similarities of COSMIC signatures to various spectra
- Perform unsupervised hierarchical clustering on signatures
- Create and save heatmap and dendrogram visuals
- Organize data to look at segments of mutations (e.g. C>A)

### Authors:
Essigmann Lab, MIT; Lina Kim, klkim [at] mit dot edu


## cospec.py
### Inputs:
Names of output file names

### Outputs:
Only when called by `visualize_CtA.py`

### Parameters:
graphics numbers

### Usage:
Called by `visualize_CtA.py`

### Dependencies:
Used on Python 3.6.0.

- collections
- scipy
- fastcluster
- numpy
- pandas
- seaborn
- matplotlib
- sklearn
- siglib.figures (Clint Valentine's library)

### Notes:
- additional functions included to strip COSMIC signatures to C>A
- use print functions for numbers on cosine similarities
- use print functions for cluster array
- TODO: fix heatmap/dendrogram functionality; only one can plot at once

## visualize_CtA.py
### Inputs:
- COSMIC mutational signatures in .txt
- 8-oxoG/TKO mutational spectrum in .txt
- (if TKO data included) TKO trinucleotide frequencies in .txt

### Outputs:
- Stratton spectra as called; in PDF
- cosine similarity heatmap
- clustering dendrogram

### Parameters:
- location of inputs
- output file names
- graphics numbers (?)

### Usage:
`python visualize_CtA.py`

### Dependencies:
Used on Python 3.6

- collections
- matplotlib
- scipy
- cospec
- siglib.figures
- all dependencies of `siglib.figures`
- all dependencies of `cospec`

### Notes:
- depends on `siglib.figures`, written by Clint Valentine; needed in directory or included in Python library - put in directory for portability

## Updates

### v1.0 (30 March 2017)
- first commit

### v1.1 (31 March 2017)
- `refine()` function added to equalize order of (mut,con) keys in COSMIC/8oxoG
- `main()` function added
- sklearn's `cosine_similarity` function used instead of scipy; also used only in aggregate cosine similarity, not one-on-one
- `take_CtA()` and `refine()` modified to allow dict or nested dict input

### v1.2 (03 April 2017)
- fix `refine()` function to actually make refinements, not just duplicate dicts

### v1.3 (10 April 2017)
- add `visualize_CtA.py` script to create Stratton spectra from data
- `visualize` depends on `siglib.figures`; added to directory

### v1.4 (02 May 2017)
- fixed matplotlib's issue with only a single heatmap or dendrogram plotting at once, made them separately called scripts
- renamed sig_8oxoG_analysis library `cospec.py`
- integrate `refine` method into spectra-creating methods (no need to be called explicitly)
