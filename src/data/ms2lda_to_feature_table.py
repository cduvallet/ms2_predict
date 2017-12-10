#!/usr/bin/env python
"""
Convert tidy ms2lda results into wide dataframe for use with classifiers.
"""

import pandas as pd
import json
import argparse

p = argparse.ArgumentParser()
p.add_argument('ms2lda_results', help='output from run_ms2lda.py')
p.add_argument('in_json', help='path to json that was used to make ms2lda '
    + 'results. The important thing is that the keys here match what is in '
    + 'the "spec" column of the ms2lda_results file, and that they have '
    + '"kingdom", "class", and "sub_class" entries.')
p.add_argument('out_table', help='path to output feature table')
args = p.parse_args()

# Read in json with metadata
with open(args.in_json, 'r') as f:
    spectra = json.load(f)

# Convert json to dataframe for later merging
spec_lst = []
for s in spectra:
    spec_lst.append(
        [s,
         spectra[s]['kingdom'],
         spectra[s]['class'],
         spectra[s]['sub_class']]
    )
spec_df = pd.DataFrame(spec_lst,
    columns=['spec', 'kingdom', 'class', 'sub_class'])

# Read in ms2lda results
ldares = pd.read_csv(args.ms2lda_results, sep='\t')
ldares = ldares.drop('motif_dup', axis=1)

# Convert each "motif" label to "motif-prob" and "motif-overlap"
resmelt = pd.melt(ldares, id_vars=['motif', 'annotation', 'spec'])
resmelt['feature'] = resmelt['motif'] + '-' + resmelt['variable']

# Pivot LDA results to wideform, features in columns and spectra in index
wideres = resmelt.pivot(index='spec', columns='feature', values='value')
wideres = wideres.fillna(0)

# Merge metadata and wide ms2lda results
wideres = pd.merge(spec_df, wideres, left_on='spec', right_index=True)

# Write to file
wideres.to_csv(args.out_table, sep='\t', index=False)
