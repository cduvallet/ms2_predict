#!/usr/bin/env python
"""
This script unpacks the CSV that contains the HMDB parsed data and writes a JSON
file with only the molecules we want to keep.

It keeps spectra which have:
- associated parent mass
- at least npeaks peaks (in the actual MS2, not the peak_counter field)
- kingdom, sub_class, and class taxonomy labels
- ionization mode

It writes a json file with:
{
 spec_id:
    {
     inchi: str, parentmass: float, kingdom: str, sub_class: str,
     class: str, ionization: str, peaks: [(mz, int), (mz, int), ...]
     }
 ...
}

where spec_id is inchikey--MS2_i--ionization (where i is the index of the
MS2 spectrum, in the original unpacked CSV).
"""

import json
import argparse

# User-defined modules
import os, sys
src_dir = os.path.normpath('src/util')
sys.path.insert(0, src_dir)
import util

p = argparse.ArgumentParser()
p.add_argument('infile', help='input csv file (to unpack)')
p.add_argument('outfile', help='json file to write output to')
p.add_argument('--npeaks', help='min number of peaks in MS2 spectrum. '
    + '[default: %(default)s]', default=3, type=int)
args = p.parse_args()

# Read in the csv data
fname = args.infile
all_mtabs = util.unpackCSV(fname)

all_spectra = {}
for m in all_mtabs:
    mtab = all_mtabs[m]

    # Only look at metabolites with at least one MS2 spectrum
    n_ms2 = len(mtab.MS2)
    if n_ms2 == 0:
        continue

    # Keep only metabolites with a parent mass
    parentmass = mtab.monisotopic_molecular_weight
    if parentmass is None:
        continue
    else:
        parentmass = float(parentmass)

    # Keep only metabolites with taxonomy
    taxonomy_dict = mtab.taxonomy_dict
    if taxonomy_dict is not None:
        kingdom = taxonomy_dict['kingdom']
        sub_class = taxonomy_dict['sub_class']
        mclass = taxonomy_dict['class']
    else:
        continue
    # Make sure at least one of the taxonomies of interest is not empty str
    if not kingdom and not sub_class and not mclass:
        continue

    ## Now that we definitely want to keep this metabolite, go through and
    ## get each of its MS2 spectra
    #print(m)
    # For each MS2 spectrum associated with that metabolite,
    for i in range(n_ms2):
        # Get the number of peaks in the spectrum
        peaks = mtab.MS2[i].peaks
        if peaks is None:
            continue
        else:
            npeaks = len(peaks)
        ionization = mtab.MS2[i].ionization_mode
        if ionization is not None and npeaks > args.npeaks:
            # Give this spectra unique ID: inchi--MS2_i--ionization_mode
            spec_id = m + '_MS2-' + str(i) + '_' + ionization

            all_spectra[spec_id] = {
                'inchi': m,
                'parentmass': parentmass,
                'kingdom': kingdom,
                 'sub_class': sub_class,
                 'class': mclass,
                 'ionization': ionization,
                 'peaks': peaks
                 }

with open(args.outfile, 'w') as f:
    json.dump(all_spectra, f)
