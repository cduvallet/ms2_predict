#!/usr/bin/env python
"""
Make binned mz feature tables
"""
import pandas as pd
import json
import copy
import argparse
import os

def remove_dup_mzs(df):
    """
    Remove duplicate mz's in df, keeping only the ones with the highest intensity.
    """
    return (df
        .sort_values(by='intensity', ascending=False)
        .drop_duplicates(subset=['inchi', 'mz'], keep='first'))

def pivot_and_add_labels(df):
    """
    Pivot a dataframe on mz and inchi (values=intensity) and add
    taxonomy labels.
    """
    widedf = df.pivot(index='inchi', columns='mz', values='intensity')
    widedf = widedf.fillna(0.0)
    widedf = pd.merge(
        widedf,
        df[['inchi', 'kingdom', 'sub_class', 'class']].drop_duplicates(),
        left_index=True,
        right_on='inchi',
        how='left')
    return widedf

p = argparse.ArgumentParser()
p.add_argument('infile', help='input json file with all spectra')
p.add_argument('outdir', help='directory to save feat tables in')
args = p.parse_args()

with open(args.infile, 'r') as f:
    all_spectra = json.load(f)

# Make tidy dataframe with spectra-related metadata, mz, intensity
dflst = []

spec_keys = ['inchi', 'ionization', 'kingdom', 'sub_class',
             'class', 'parentmass']

for spec_id, spectrum in all_spectra.iteritems():
    # Get the spectrum-related metadata
    spec_metadata = [spectrum[k] for k in spec_keys]
    # Get the spectra number from the label
    spec_metadata += [spec_id.split('_')[1].split('-')[1]]
    # Add each mz, int pair as its own entry
    dflst += [spec_metadata + [p[0], p[1]] for p in spectrum['peaks']]
df = pd.DataFrame(dflst, columns=spec_keys + ['scan_id', 'mz', 'intensity'])

# Subset by ionization mode
pos = df.query('ionization == "Positive"')
neg = df.query('ionization == "Negative"')
# and df has all of the scans, including the ones labeled N/A

# Remove duplicate mz's and merge all scans per molecule
pos = remove_dup_mzs(pos)
neg = remove_dup_mzs(neg)
both = remove_dup_mzs(df)

# Convert to wideform and save to disk
widepos = pivot_and_add_labels(pos)
fname = os.path.join(args.outdir, 'raw_mz.positive.txt')
widepos.to_csv(fname, sep='\t', index=False)

wideneg = pivot_and_add_labels(neg)
fname = os.path.join(args.outdir, 'raw_mz.negative.txt')
wideneg.to_csv(fname, sep='\t', index=False)

wideboth = pivot_and_add_labels(both)
fname = os.path.join(args.outdir, 'raw_mz.all_scans.txt')
wideboth.to_csv(fname, sep='\t', index=False)

### mz's rounded to the nearest integer
def round_and_save(df, fname):
    df = copy.deepcopy(df)
    df['mz'] = df['mz'].astype(int)
    df = remove_dup_mzs(df)
    widedf = pivot_and_add_labels(df)
    widedf.to_csv(fname, sep='\t', index=False)

fname = os.path.join(args.outdir, 'mz_integer.positive.txt')
round_and_save(pos, fname)

fname = os.path.join(args.outdir, 'mz_integer.negative.txt')
round_and_save(neg, fname)

fname = os.path.join(args.outdir, 'mz_integer.all_scans.txt')
round_and_save(both, fname)
