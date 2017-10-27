#!/usr/bin/env python
"""
This script reads in the xml file containing all MS2 spectra and the xml
file containing HMDB metabolite metadata, combines them based on inchikeys,
and writes them to a csv.
"""
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import numpy as np
import csv
from unidecode import unidecode
import argparse

# Assuming that you're calling this script from the top directory in the repo,
# as the Makefile does, these statements add src/util to the path so we
# can import modules from files we've written and placed there
import os
import sys
src_dir = os.path.normpath(os.path.join(os.getcwd(), 'src/util'))
sys.path.insert(0, src_dir)
from MetabolomicsObjects import Metabolite, MS2, MSPeak

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('ms2_concat', help='path to file with all MS2s '
        + 'concatenated into one xml file.')
    p.add_argument('metabolites_info', help='path to file with all '
        + 'HMDB metabolites in one xml file.')
    p.add_argument('out', help='path to write output csv file to.')
    return p.parse_args()

def metabolitePreprocessing(xml_file, metabolite_feature_set):
    '''
    Reads in metabolite metadata xml file and generates Metabolite objects. All
        Metabolite objects have intentionally unpopulated MS2 fields.
    Populates metabolite attributes as one of:
        metabolite.attribute = value
        metabolite.attribute = [value1, value2, value3...]
        metabolite.attribute = {key1:item1, key2:item2, ...}

    Attributes populated:
      - accession, secondary accessions (list), name, chemical_formula,
        monisotopic_molecular_weight, iupac_name, traditional_iupac,
        cas_registry_number, smiles, inchi, inchikey
      - taxonomy, except substitutents, parents, and alternative parents
      - biofluid locations
      - IDs (features with _id in the name or "id" as the name)

    Args:
        xml_file: path to metabolite metadata xml_file
    Returns:
        Populated dictionary in the format of {inchikey:Metabolite}
    '''
    print 'Generating Metabolite XML Tree...'
    tree = ET.parse(xml_file)
    print 'Done. \nPreprocessing Metabolite Data...'
    root = tree.getroot()
    metabolite_dict = dict()
    for metabolite_tag in root:
        metabolite = Metabolite()
        id_dict = dict()
        #xml_object = highest level tag within <c-ms> tag
        for xml_object in metabolite_tag:
            clean_tag = xml_object.tag.replace('{http://www.hmdb.ca}', '')
            #filter out unwanted metadata
            if not (clean_tag in metabolite_feature_set or '_id' in clean_tag):
                continue
            #collect secondary accessions in list
            if clean_tag == 'secondary_accessions':
                secondary_accessions = []
                for accessions_object in xml_object:
                    secondary_accessions.append(accessions_object.text)
                metabolite.secondary_accessions = secondary_accessions
            #collect biofluid locations in list
            elif clean_tag == 'biofluid_locations':
                locations = []
                for location_object in xml_object:
                    locations.append(location_object.text)
                metabolite.biofluid_locations = locations

            #collect all ID's in a dictionary
            elif '_id' in clean_tag:
                id_dict[clean_tag] = xml_object.text

            #collect taxonomy metadata in dictionary
            elif clean_tag == 'taxonomy':
                taxonomy_dict = dict()
                for taxonomy_object in xml_object:
                    clean_taxonomy_tag = taxonomy_object.tag.replace('{http://www.hmdb.ca}', '')
                    if clean_taxonomy_tag != 'substituents' and 'parent' not in clean_taxonomy_tag:
                        taxonomy_dict[clean_taxonomy_tag] = taxonomy_object.text
                metabolite.taxonomy_dict = taxonomy_dict

            #collect miscellaneous, desireable attributes
            else:
                if 'inchi' in clean_tag and 'key' in clean_tag:
                    setattr(metabolite, 'inchikey', xml_object.text)
                    setattr(metabolite, 'inchi_key', xml_object.text)
                else:
                    setattr(metabolite, clean_tag, xml_object.text)

        metabolite.id_dict = id_dict
        metabolite_dict[metabolite.inchikey] = metabolite
    print 'Done.'
    return metabolite_dict

def MS2Preprocessing(xml_file, feature_set, metabolite_dict):
  '''
  Takes in desired features and metabolite dictionary (from
  metabolitePreprocessing), reads through MS2 metadata, populates MS2 objects,
  and populates the ms2 field of the Metabolite objects in the metabolite
  dictionary with these MS2 objects. Filters for MS2 data generated with
  ms-ms peaks.

  Attributes populated as one of:
    MS2.attribute = value
    MS2.attribute = {key1:item1, key2:item2, ...}
    MS2.ms_peak = [MSPeak1, MSPeak2, ...]

  Attributes populated:
    -  IDs (features with _id in the name or "id" as the name)
    -  References (anything under the references tag)
    -  ms-ms peaks (anything under the ms-ms-peaks tag)
    -  Other features from the feature_set input

  Args:
    xml_file: path to MS2 metadata xml file
    feature_set: hash set of desired features with which to populate MS2
      object fields
    metabolite_dict: dictionary formatted like: {inchikey: Metabolite}; the
      output of metabolitePreprocessing()
  Returns:
    Dictionary formatted like {inchikey: Metabolite} where the Metabolite
      has the MS2 field populated
    List of MS2 objects for which the input metabolite_dict did not have a
      matching inchikey key
  '''
  print 'Generating MS2 XML Tree...'
  tree = ET.parse(xml_file)
  print 'Done. \nPreprocessing MS2 Data...'
  root = tree.getroot()
  failure_list = []
  for metabolite_analysis in root:
    ms2_object = MS2()
    id_dict = dict()
    failure = False
    for feature in metabolite_analysis:
      #filter out non-ms-ms data
      if 'peaks' in feature.tag and feature.tag != 'ms-ms-peaks':
        failure = True
        break
      #filter out non-GC chromatagraphy data
      #if feature.tag == 'chromatography-type' and feature.text != 'GC':
      #  failure = True
      #  break

      #collect all id data in dictionary
      if '_id' in feature.tag or feature.tag == 'id' or '-id' in feature.tag:
        id_dict[feature.tag.replace('-', '_')] = feature.text

      #collect references data in dictionary
      elif feature.tag == 'references':
        references_dict = dict()
        for reference_tag in feature:
          for references_feature in reference_tag:
            if references_feature.text:
              references_dict[references_feature.tag.replace('-', '_')] = references_feature.text
            else:
              references_dict[references_feature.tag.replace('-', '_')] = ''
        ms2_object.references_dict = references_dict

      #collect peak data as list of MSPeak objects
      elif feature.tag == 'ms-ms-peaks':
        peak_objects = []
        for peak_object in feature:
          ms_peak = MSPeak()
          for peak_object_attribute in peak_object:
            if peak_object_attribute.text:
              setattr(ms_peak, peak_object_attribute.tag.replace('-', '_'), peak_object_attribute.text)
            else:
              setattr(ms_peak, peak_object_attribute.tag.replace('-', '_'), '')
          peak_objects.append(ms_peak)
        ms2_object.peaks = peak_objects

      #collect other data attributes of interest
      else:
        if feature.text:
          setattr(ms2_object, feature.tag.replace('-', '_'), feature.text)
        else:
          setattr(ms2_object, feature.tag.replace('-', '_'), '')
        if 'inchi' in feature.tag and 'key' in feature.tag:
          setattr(ms2_object, 'inchikey', feature.text)
          setattr(ms2_object, 'inchi_key', feature.text)

    ms2_object.id_dict = id_dict
    #check for failure i.e. this data should be ignored
    if failure: continue
    #add MS2 feature to metabolite dictionary
    if metabolite_dict.get(ms2_object.inchi_key):
      metabolite_dict[ms2_object.inchi_key].MS2.append(ms2_object)
    else:
      failure_list.append(ms2_object)

  print 'Done.'
  print 'Found', len(failure_list), 'MS2 objects without Metabolite pairs'
  return metabolite_dict, failure_list


def writeToCSV(matched_dict, file_path):
  '''
  Takes in output from MS2Preprocessing() and writes dictionary to CSV.
  Format is one of:
    attribute=value
    attribute=value1,value2,value3,...
    attribute=key1:value1,key2:value2,...
  Each CSV cell holds one Metabolite field or MS2 field, regardless of field type

  Args:
    matched_dict: dictionary in the form of: {inchikey: Metabolite}, where
      Metabolite object has a populated MS2 field. This is the output of
      MS2Preprocessing().
    file_path: path to which the csv file will be written
  '''
  print 'Writing to CSV...'
  csv_file = open(file_path, 'w+')
  writer = csv.writer(csv_file)
  #for each metabolite
  for metabolite in matched_dict.itervalues():
    attributes = dir(metabolite)
    attributes = [attribute for attribute in attributes if '__' not in attribute and 'MS2' not in attribute]
    writestring_list = [metabolite.inchikey]
    for attribute in attributes:
      value = getattr(metabolite, attribute)
      writestring = str(attribute) + '='
      #if empty field
      if not value:
        pass
      #if field=value format
      elif type(value) == str:
        writestring += value
      # If attribute contains a dictionary, write it as key:value;;key2:value2
      elif type(value) == dict:
        for item in value.iteritems():
          writestring += item[0] + ':'
          if item[1]:
            if type(item[1]) == unicode:
              writestring += unidecode(item[1]) + ';;'
            else:
              writestring += item[1] + ';;'
          else:
            writestring += ';;'
        writestring = writestring[:-1]
      #if field=[value, value, ...] format
      elif type(value) == list:
        for list_item in value:
          writestring += list_item + ','
        writestring = writestring[:-1]
      writestring_list.append(writestring)
    writer.writerow(writestring_list)

    ms2_object_list = metabolite.MS2
    # Write each MS2 object on a new row
    for ms2_object in ms2_object_list:
      attributes = [attribute for attribute in dir(ms2_object)
                    if '__' not in attribute]
      writestring_list = ['']
      for attribute in attributes:
        value = getattr(ms2_object, attribute)
        writestring = str(attribute) + '='
        #if empty field
        if not value:
          pass
        #if field=value format
        elif type(value) == str:
          writestring += value
        #if field={key:value, ...} format
        elif type(value) == dict:
          for item in value.iteritems():
            writestring += item[0] + ':'
            if item[1]:
              if type(item[1]) == unicode:
                writestring ++ unidecode(item[1]) + ','
              else:
                writestring += item[1] + ','
            else:
              writestring += ','
          writestring = writestring[:-1]
        #if field=[value, value, ...] format
        elif type(value) == list:
          if str(value[0]) != 'MSPeak':
            for list_item in value:
              writestring += list_item + ','
            writestring = writestring[:-1]
          else:
            ms2_attributes = [attribute for attribute
                              in dir(value[0]) if '__' not in attribute]
            for peak in value:
              for ms_attribute in ms2_attributes:
                try:
                  ms_value = getattr(peak, ms_attribute)
                  if not ms_value:
                    ms_value = ''
                except AttributeError:
                  ms_value = ''
                writestring += '!'
                writestring += str(ms_attribute) + '=' + str(ms_value)
            writestring += ','
            writestring = writestring[:-1]

        writestring_list.append(writestring)
      writer.writerow(writestring_list)
  csv_file.close()

if __name__ == '__main__':
  args = parse_args()
  #######    TO GENERATE CSV    #######
  ms2_feature_set =  set(
    ['inchi_key', 'frequency', 'instrument_type', 'id', 'energy_field',
     'base_peak', 'sample_mass_units', 'chromatography_type', 'searchable',
     'sample_mass', 'derivative_mw', 'retention_time', 'updated_at',
     'sample_assessment', 'derivative_formula', 'derivative_type',
     'database_id', 'ref_text', 'mass_charge', 'collision_energy_voltage',
     'sample_concentration', 'spectra_assessment', 'solvent', 'nucleus_y',
     'ionization_mode', 'collection_date', 'nucleus',
     'sample_temperature_units', 'sample_ph', 'spectra_id', 'c_ms_id',
     'sample_concentration_units', 'sample_source', 'nil_classes', 'nucleus_x',
     'database', 'notes', 'created_at', 'sample_temperature', 'ri_type',
     'pubmed_id', 'molecule_id', 'column_type', 'retention_index',
     'collision_energy_level', 'references', 'name', 'accession',
     'chemical_formula', 'monoisotopic_molecular_weight', 'iupac_name',
     'traditional_iupac', 'cas_registry', 'smiles', 'inchi', 'inchikey',
     'taxonomy', 'biofluid_locations', 'ids', 'peak_counter', 'mono_mass',
     'ms_ms_id', 'spectra_type'])

  # Set up some file paths
  ms2_xml_file = args.ms2_concat
  metabolite_xml_file = args.metabolites_info
  metabolite_feature_set = set(['accession', 'secondary_accessions', 'name',
                'chemical_formula',
                'monisotopic_molecular_weight', 'iupac_name',
                'traditional_iupac', 'cas_registry_number',
                'smiles', 'inchi', 'inchikey',
                'biofluid_locations', 'taxonomy'])

  metabolite_dict = metabolitePreprocessing(
    metabolite_xml_file, metabolite_feature_set)
  matched_dict, failure_list = MS2Preprocessing(
    ms2_xml_file, ms2_feature_set, metabolite_dict)

  csv_file_path = args.out
  writeToCSV(matched_dict, csv_file_path)
