#!/usr/bin/env python
"""
Convert feature table with merged spectra into a json that can be
parsed and run through ms2lda
"""

import pandas as pd
import json
import argparse

# If there are duplicate peaks, we'll just pick the one with the highest intensity
def remove_dup_mzs(df):
    """
    Remove duplicate mz's in df, keeping only the ones with the highest intensity.
    """
    return (df
        .sort_values(by='intensity', ascending=False)
        .drop_duplicates(subset=['inchi', 'mz'], keep='first'))

p = argparse.ArgumentParser()
p.add_argument('injson', help='input json file, with each individual spectrum '
    + ' as a separate entry.')
p.add_argument('outjson', help='file to write json with merged spectra to')
args = p.parse_args()

# Read in the unconcatenated data
with open(args.injson, 'r') as f:
    spectra = json.load(f)

# Make tidy dataframe with spectra-related metadata, mz, intensity
dflst = []

spec_keys = ['inchi', 'ionization', 'kingdom', 'sub_class',
             'class', 'parentmass']

for spec_id, spectrum in spectra.iteritems():
    # Get the spectrum-related metadata
    spec_metadata = [spectrum[k] for k in spec_keys]
    # Get the spectra number from the label
    spec_metadata += [spec_id.split('_')[1].split('-')[1]]
    # Add each mz, int pair as its own entry
    dflst += [spec_metadata + [p[0], p[1]] for p in spectrum['peaks']]
df = pd.DataFrame(dflst, columns=spec_keys + ['scan_id', 'mz', 'intensity'])

# Remove duplicate mz's
df = remove_dup_mzs(df)

# Now convert back to json
spec_dict = {}
for g, subdf in df.groupby('inchi'):
    spec_dict.update({g: {
        'parentmass': subdf['parentmass'].iloc[0],
        'kingdom': subdf['kingdom'].iloc[0],
        'class': subdf['class'].iloc[0],
        'sub_class': subdf['sub_class'].iloc[0],
        'peaks': [(i, j) for i, j in zip(subdf['mz'], subdf['intensity'])]
    }})

with open(args.outjson, 'w') as f:
    json.dump(spec_dict, f)
