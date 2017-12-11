#!/usr/bin/env
"""
This script pulls the molecules from our CSV and runs them through MS2LDA.
"""
import json
import requests
import argparse
import time
import pandas as pd

# User-defined modules
import os, sys
src_dir = os.path.normpath('src/util')
sys.path.insert(0, src_dir)
import util

p = argparse.ArgumentParser()
p.add_argument('infile', help='input json file, spectra as main keys')
p.add_argument('outfile', help='outfile')
args = p.parse_args()

# Read in the csv data
with open(args.infile, 'r') as f:
    spectra = json.load(f)

ms2lda_spectra = []
for spec in spectra.keys():
    ms2lda_spectra.append(
        (spec,
         spectra[spec]['parentmass'],
         [tuple(p) for p in spectra[spec]['peaks']])
    )

# Run ms2lda - post request first
ms2lda_args = {'spectra': json.dumps(ms2lda_spectra),
              'motifset': 'massbank_motifset'}
url = 'http://ms2lda.org/decomposition/api/batch_decompose/'

r = requests.post(url, ms2lda_args)
t0 = time.time()
print("Request posted at {}.".format(t0))

# Check for results
result_id = r.json()['result_id']
url2 = 'http://ms2lda.org/decomposition/api/batch_results/{}/'.format(result_id)
print(url2)
r2 = requests.get(url2)

# Keep checking if results are back before doing anything else
while 'status' in r2.json():
    time.sleep(30)
    r2 = requests.get(url2)
t1 = time.time()
print('MS2LDA took {:.2f} s for {} spectra'.format(
        t1 - t0, len(ms2lda_spectra)))

# Convert results into dataframe
res = r2.json()['decompositions']
reslst = []
for k in res:
    tmp = pd.DataFrame(res[k],
             columns=['motif', 'motif_dup', 'prob', 'overlap', 'annotation'])
    tmp['spec'] = k
    reslst.append(tmp)

ldadf = pd.concat(reslst)
ldadf.to_csv(args.outfile, sep='\t', encoding="utf-8", index=False)
