# Makefile for 6.867 project
# Authors: Julian Alverio, Ashvin Bashyam, Claire Duvallet

all: data

## DATA
RAW = data/raw
CLEAN = data/clean

SRCDATA = src/data

## Define data files
# MS2 spectra files
raw_ms2 = $(RAW)/hmdb_spectra_xml.zip
concat_ms2_clean = $(CLEAN)/hmdb_spectra.concat.clean.xml
concat_ms2_tmp = $(RAW)/hmdb_spectra.concat.excess_decs.xml

# HMDB metabolite info
raw_hmdb = $(RAW)/hmdb_metabolites.zip
metabolites_clean = $(CLEAN)/hmdb_metabolites.xml

# CSV version of xml files
csv_data = $(CLEAN)/metabolites_and_spectra.csv

# JSON version of csv data, with curated metabolites
json_data = $(CLEAN)/clean_spectra.json

# Feature table with all scans (positive, negative, and n/a) merged (duplicate mz's removed, highest intensity peak retained)
feat_table = data/feature_tables/raw_mz.all_scans.txt

data: $(concat_ms2_clean) $(concat_hmdb_clean) $(csv_data) $(json_data)

################################
#                              #
#          MAKE DATA           #
#                              #
################################

# Download and unzip raw XML files from HMDB. The touch commands are to update
# the timestamp to the time of download
$(raw_ms2):
	wget -O $@ http://specdb.wishartlab.com/downloads/exports/spectra_xml/hmdb_spectra_xml.zip
	touch $@
	unzip $@ -d $(RAW)/ms2_xmls/

$(raw_hmdb):
	wget -O $@ http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
	touch $@

# Concatenate them using find (bc otherwise argument list is too long for shell)
# and remove the excess declarations from the concatenated file
$(concat_ms2_clean): $(raw_ms2) \
			$(SRCDATA)/eliminate_remove_excess_xml_declarations.py
	find $(RAW)/ms2_xmls/ -name '*.xml'  -exec cat {} + > $(concat_ms2_tmp)
	python $(SRCDATA)/eliminate_remove_excess_xml_declarations.py $(concat_ms2_tmp) $@

$(metabolites_clean): $(raw_hmdb)
	unzip $(raw_hmdb) -d $(CLEAN)
	touch $@

# Use xml_parser.py to combine the metabolite metadata with spectra info
$(csv_data): $(SRCDATA)/parse_xml_files.py $(concat_ms2_clean) $(metabolites_clean)
	python $< $(concat_ms2_clean) $(metabolites_clean) $@

# Convert csv into easy-to-read json with only metabolites of interest
$(json_data): $(SRCDATA)/clean_csv.py $(csv_data)
	python $< $(csv_data) $@ --npeaks 3

# Make similar json, but with the merged spectra (inchikey as the top keys)
$(merged_json): $(SRCDATA)/make_merged_json.py $(json_data)
	python $< $(json_data) $(merged_json)

# MS2LDA results on the non-collapsed spectra
$(ms2lda_data): $(SRCDATA)/run_ms2lda.py $(json_data)
	python $< $(json_data) $@

# MS2LDA results on collapsed spectra
$(collapsed_ms2lda_data): $(SRCDATA)/run_ms2lda.py $(merged_json)
	python $< $(merged_json) $@

### FEATURE TABLES

## Convert json to positive, negative, and all_scans feature tables
# python src/data/make_mz_feature_tables.py data/clean/clean_spectra.json data/feature_tables/

## MS2LDA feature tables
#python src/data/ms2lda_to_feature_table.py data/clean/ms2lda_results.txt data/clean/clean_spectra.json data/feature_tables/ms2lda_feature_table.txt

#python src/data/ms2lda_to_feature_table.py data/clean/ms2lda_results.merged_spectra.txt data/clean/clean_spectra.by_inchi.json data/feature_tables/ms2lda_feature_table.merged_spectra.txt
