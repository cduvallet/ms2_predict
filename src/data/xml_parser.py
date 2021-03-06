try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import numpy as np
import csv
import os

import sys
sys.path.append('/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
from unidecode import unidecode


class MS2:
    '''
    A container for all MS2 metadata. Final MS2 objects may have empty
    or additional fields, depending on any anomalies in the data.
    '''
    def __init__(self):
        self.peak_counter = None
        self.mono_mass = None
        self.spectra_type = None
        self.inchi_key = None
        self.frequency = None
        self.instrument_type = None
        self.energy_field = None
        self.base_peak = None
        self.sample_mass_units = None
        self.chromatography_type = None
        self.searchable = None
        self.sample_mass = None
        self.derivative_mw = None
        self.retention_time = None
        self.updated_at = None
        self.sample_assessment = None
        self.derivative_formula = None
        self.derivative_type = None
        self.database_id = None
        self.ref_text = None
        self.mass_charge = None
        self.collision_energy_voltage = None
        self.sample_concentration = None
        self.spectra_assessment = None
        self.solvent = None
        self.nucleus_y = None
        self.ionization_mode = None
        self.collection_date = None
        self.nucleus = None
        self.sample_temperature_units = None
        self.sample_ph = None
        self.spectra_id = None
        self.c_ms_id = None
        self.sample_concentration_units = None
        self.sample_source = None
        self.nil_classes = None
        self.nucleus_x = None
        self.database = None
        self.notes = None
        self.created_at = None
        self.sample_temperature = None
        self.ri_type = None
        self.pubmed_id = None
        self.molecule_id = None
        self.column_type = None
        self.retention_index = None
        self.collision_energy_level = None
        self.references = None
        self.name = None
        self.accession = None
        self.chemical_formula = None
        self.monoisotopic_molecular_weight = None
        self.iupac_name = None
        self.traditional_iupac = None
        self.cas_registry = None
        self.smiles = None
        self.inchi = None
        self.inchikey = None
        self.taxonomy = None
        self.biofluid_locations = None
        self.id_dict = None
        self.peaks = None

    def __str__(self):
        return 'MS2'


class MSPeak:
    '''
    Container for data about an ms-ms peak. Final MSPeak objects may have
    empty or additional fields, depending on any anomalies in the data.
    '''
    def __init__(self):
        self.id = None
        self.ms_ms_id = None
        self.intensity = None
        self.mass_charge = None

    def __str__(self):
        return 'MSPeak'


class Metabolite:
    '''
    Container for metadata about a metabolite. Final Metabolite objects may have
    empty or additional fields, depending on any anomalies in the data.
    '''
    def __init__(self):
        self.accession = None
        self.secondary_accessions = None
        self.name = None
        self.chemical_formula = None
        self.monisotopic_molecular_weight = None
        self.iupac_name = None
        self.traditional_iupac = None
        self.cas_registry = None
        self.smiles = None
        self.inchi = None
        self.inchikey = None
        self.biofluid_locations = None
        self.id_dict = None
        self.taxonomy_dict = None
        self.MS2 = []

    def __str__(self):
        return 'Metabolite'


def findAllTags(xml_file):
    '''
    Depth-first search implementation. Used to find all unique tags in the MS2
    metadata. Returns list of unique tags to allow for designation of desired
    metadata attributes.

    Args:
        xml_file: file path to xml_file of which to find attributes

    Returns:
        list of unique attributes in xml file
    '''
    print 'Generating XML Tree...'
    tree = ET.parse(xml_file)
    print 'Done. \nFinding All Tags...'
    root = tree.getroot()
    discovered_set = set()
    queue = [child for child in root]
    while queue:
        current_node = queue.pop()
        children = [child for child in current_node]
        queue.extend(children)
        for child in children:
            try:
                discovered_set.add(child.tag)
            except TypeError:
                import pdb; pdb.set_trace()
    return list(discovered_set)


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
  ms-ms peaks and GC chromotography.

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
      #if field={key:value, ...} format
      elif type(value) == dict:
        for item in value.iteritems():
          writestring += item[0] + ':'
          if item[1]:
            if type(item[1]) == unicode:
              writestring += unidecode(item[1]) + ','
            else:
              writestring += item[1] + ','
          else:
            writestring += ','
        writestring = writestring[:-1]
      #if field=[value, value, ...] format
      elif type(value) == list:
        for list_item in value:
          writestring += list_item + ','
        writestring = writestring[:-1]
      writestring_list.append(writestring)
    writer.writerow(writestring_list)

    ms2_object_list = metabolite.MS2
    #for MS2 object
    for ms2_object in ms2_object_list:
      attributes = dir(ms2_object)
      attributes = [attribute for attribute in attributes if '__' not in attribute]
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
            ms2_attributes = [attribute for attribute in dir(value[0]) if '__' not in attribute]
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

def unpackCSV(file_path):
  '''
  Reads from csv, generates {inchikey:Metabolite} dictionary w/populated MS2

  Args:
    file_path: path to csv file to read
  Returns:
    {inchikey:Metabolite} dictionary
  '''
  csv_file = open(file_path, 'r')
  reader = csv.reader(csv_file)
  row_count = 0
  metabolite_dict = dict()
  for row in reader:
    row_count += 1
    #Row is for a Metabolite object
    if row[0]:
      metabolite = Metabolite()
      ms2_object_lst = []
      for cell in row:
        #dictionary in cell
        if ':' in cell:
          attribute_dictionary = dict()
          attribute = cell.split('=')[0]
          dictionary = cell.split('=')[1]
          split_dictionary = dictionary.split(',')
          for item in split_dictionary:
            key = item.split(':')[0]
            if len(item.split(':')) == 2:
              value = item.split(':')[1]
            else:
              value = None
            attribute_dictionary[key] = value
          setattr(metabolite, attribute, attribute_dictionary)
        #list in cell
        elif ',' in cell:
          attribute = cell.split('=')[0]
          cell_list = cell.split('=')[1]
          split_list = cell_list.split(',')
          setattr(metabolite, attribute, split_list)
        #sole (attribute, value) pair in cell but not just cell w/ only inchikey
        else:
          if '=' not in cell:
            setattr(metabolite, 'inchikey', cell)
            setattr(metabolite, 'inchi_key', cell)
          else:
            attribute = cell.split('=')[0]
            value = cell.split('=')[1]
            if value:
              setattr(metabolite, attribute, value)
            else:
              setattr(metabolite, attribute, None)
    # If there is nothing in first cell, row is for an MS2 object
    else:
      ms2_object = MS2()
      for cell in row[1:]:
        #dictionary in cell:
        if ':' in cell:
          attribute_dictionary = dict()
          attribute = cell.split('=')[0]
          dictionary = cell.split('=')[1]
          split_dictionary = dictionary.split(',')
          for item in split_dictionary:
            key = item.split(':')[0]
            if len(item.split(':')) == 2:
              value = item.split(':')[1]
            else:
              value = None
            attribute_dictionary[key] = value
          setattr(ms2_object, attribute, attribute_dictionary)
        #ms_peak data in cell
        elif 'ms_peak' in cell and '!' in cell:
          peak_list = cell.split('=')[1].split(',')
          for peak in peak_list:
            peak_attributes = peak_list.split('!')
            peak_object = MSPeak()
            for peak_item in peak_items:
              peak_attribute = peak_items.split('=')[0]
              peak_value = peak_items.split('=')[1]
              setattr(peak_object, peak_attribute, peak_value)
        #sole attribute data in cell
        else:
          peak_attribute = cell.split('=')[0]
          peak_value = cell.split('=')[1]
          setattr(ms2_object, peak_attribute, peak_value)
      ms2_object_lst.append(ms2_object)
    metabolite.MS2 = ms2_object_lst
    metabolite_dict[metabolite.inchikey] = metabolite

  return metabolite_dict

def testMetabolitePreprocessing(testingCSV=False):
    '''
    Tests metabolitePreprocessing()

    Cases Tested:
        Multiple Metabolite objects
        Handling biofluid_locations
        Handling IDs
        Handling taxonomy
        Handling secondary_accessions
        Handling other features in feature set
        Filturing out irrelevant features not in feature set

    Returns:
        Boolean indicating equality
        metabolite_dict in {inchikey: Metabolite}, if testingCSV if true
    '''
    print 'Verifying Correctness of metabolitePreprocessing()...'
    metabolites_xml_path = os.getcwd() + '/tests/metabolites_test.xml'
    metabolite_feature_set = set(['inchikey','database_id','taxonomy','biofluid_locations', 'secondary_accessions', 'accession'])
    metabolite_dict = metabolitePreprocessing(metabolites_xml_path, metabolite_feature_set)

    correct_metabolite_dict = dict()
    correct_metabolite_1 = Metabolite()
    correct_metabolite_1.accession = 'testaccession'
    correct_metabolite_1.secondary_accessions = ['testsecondaryaccession1', 'testsecondaryaccession2']
    correct_metabolite_1.biofluid_locations = ['testlocation1', 'testlocation2']
    correct_metabolite_1.inchikey = 'testinchikey'
    correct_metabolite_1.inchi_key = 'testinchikey'
    taxonomy_dict = dict()
    taxonomy_dict['description'] = 'testdescription'
    taxonomy_dict['kingdom'] = 'testkingdom'
    taxonomy_dict['super_class'] = 'testsuperclass'
    taxonomy_dict['class'] = 'testclass'
    taxonomy_dict['sub_class'] = 'testsubclass'
    taxonomy_dict['molecular_framework'] = 'testmolecularframwork'
    correct_metabolite_1.taxonomy_dict = taxonomy_dict
    id_dict = dict()
    id_dict['foodb_id'] = 'testfoodbid'
    id_dict['knapsack_id'] = 'testknapsackid'
    correct_metabolite_1.id_dict = id_dict
    correct_metabolite_dict[correct_metabolite_1.inchikey] = correct_metabolite_1

    correct_metabolite_2 = Metabolite()
    correct_metabolite_2.accession = 'testaccession'
    correct_metabolite_2.inchikey = 'testinchikey2'
    correct_metabolite_2.inchi_key = 'testinchikey2'
    correct_metabolite_2.id_dict = dict()
    correct_metabolite_dict[correct_metabolite_2.inchikey] = correct_metabolite_2

    if testingCSV: return metabolite_dict

    result1 = compareMetabolites(correct_metabolite_dict['testinchikey'], metabolite_dict['testinchikey'])
    result2 = compareMetabolites(correct_metabolite_dict['testinchikey2'], metabolite_dict['testinchikey2'])
    return result1 and result2


def compareMetabolites(correct_metabolite, generated_metabolite):
    test_attributes = [attribute for attribute in dir(generated_metabolite) if '__' not in attribute]
    correct_attributes = [attribute for attribute in dir(correct_metabolite) if '__' not in attribute]
    success = True
    try:
        assert set(test_attributes) == set(correct_attributes)
    except AssertionError:
        success = False
        print 'FAILURE. Generated set of attributes %s did not match the expected set: %s' % (test_attributes, correct_attributes)
    for attribute in test_attributes:
        try:
            if attribute != 'MS2':
                if getattr(correct_metabolite, attribute) or getattr(generated_metabolite, attribute):
                    assert getattr(correct_metabolite, attribute) == getattr(generated_metabolite, attribute)
        except AssertionError:
            success = False
            print 'FAILURE. Generated %s Value of %s did not match the expected %s' % (attribute, getattr(generated_metabolite, attribute), getattr(correct_metabolite, attribute))

    if not correct_metabolite.MS2:
        try:
            assert not generated_metabolite.MS2
        except AssertionError:
            success = False
            print "Failure. Unexpected MS2 object found"
    else:
        if compareMS2(correct_metabolite.MS2, generated_metabolite.MS2):
            success = False
            print "Failure. MS2s were not equal"
    return success


def compareMS2(ms2_1, ms2_2):
    '''
    Evaluates the equality of two MS2 objects.

    Args:
        ms2_1: an MS2 object
        ms2_2: an MS2 object

    Returns:
        Boolean, indicating equality or lack thereof

    '''
    if type(ms2_1) != type(ms2_2):
        return False
    if (set(dir(ms2_1)) != set(dir(ms2_2))):
        return False
    try:
        attributes = [attribute for attribute in dir(ms2_1) if '__' not in attribute]
        for attribute in attributes:
            if getattr(ms2_1, attribute) or getattr(ms2_2, attribute):
                if getattr(ms2_1, attribute) != getattr(ms2_2, attribute):
                    return False
    except:
        return False
    return True


def testMS2Preprocessing():
    '''
    Tests MS2Preprocessing()

    Cases Tested:
        Filtering out non-ms-ms data
        Filtering out non-GC chromatagraphy data
        Handling of ID data
        Handling of references data
        Handling of ms-ms-peak data
        Handling of miscellaneous metadata

    Returns:
        Boolean, whether or not test cases passed
    '''
    ms2_feature_set = set(['inchi_key','frequency','instrument_type','id','energy_field','base_peak','sample_mass_units','chromatography_type','searchable','sample_mass','derivative_mw','retention_time','updated_at','sample_assessment','derivative_formula','derivative_type','database_id','ref_text','mass_charge','collision_energy_voltage','sample_concentration','spectra_assessment','solvent','nucleus_y','ionization_mode','collection_date','nucleus','sample_temperature_units','sample_ph','spectra_id','c_ms_id','sample_concentration_units','sample_source','nil_classes','nucleus_x','database','notes','created_at','sample_temperature','ri_type','pubmed_id','molecule_id','column_type','retention_index','collision_energy_level','references','name','accession','chemical_formula','monoisotopic_molecular_weight','iupac_name','traditional_iupac','cas_registry','smiles','inchi','inchikey','taxonomy','biofluid_locations','ids','peak_counter', 'mono_mass', 'ms_ms_id', 'spectra_type'])
    ms2_xml_path = os.getcwd() + '/tests/ms2_test.xml'

    metabolite_dict = dict()
    metabolite_dict['test'] = Metabolite()
    metabolite_dict['testinchikey'] = Metabolite()
    metabolite_dict['testinchikey2'] = Metabolite()

    matched_dict, failure_list = MS2Preprocessing(ms2_xml_path, ms2_feature_set, metabolite_dict)

    correct_matched_dict = dict()
    correct_ms2_1 = MS2()
    id_dict = dict()
    id_dict['database_id'] = 'testid'
    id_dict['id'] = '1'
    correct_ms2_1.id_dict = id_dict
    correct_ms2_1.energy_field = ''
    correct_ms2_1.inchi_key = 'testinchikey'
    correct_ms2_1.inchikey = 'testinchikey'
    references_dict = dict()
    references_dict['database'] = 'HMDB'
    references_dict['database_id'] = 'testid'
    references_dict['ref_text'] = ''
    correct_ms2_1.references_dict = references_dict
    peak1 = MSPeak()
    peak1.id = '1'
    peak1.mass_charge = '1.1'
    peak1.molecule_id = ''
    peak2 = MSPeak()
    peak2.id = '2'
    correct_ms2_1.peaks = [peak1, peak2]
    correct_matched_dict[correct_ms2_1.inchi_key] = Metabolite()
    correct_matched_dict[correct_ms2_1.inchi_key].MS2 = [correct_ms2_1]

    correct_ms2_2 = MS2()
    correct_ms2_2.id_dict = {'id': '2'}
    correct_ms2_2.inchi_key = 'testinchikey2'
    correct_ms2_2.inchikey = 'testinchikey2'
    correct_matched_dict[correct_ms2_2.inchi_key] = Metabolite()
    correct_matched_dict[correct_ms2_2.inchi_key].MS2 = [correct_ms2_2]

    correct_matched_dict['test'] = Metabolite()

    success = True
    try:
        assert set(matched_dict.keys()) == set(correct_matched_dict.keys())
    except AssertionError:
        success = False
        print 'Failure. Generated dictionary keys %s did not match the expected: %s' % (matched_matched_dict.keys(), correct_dict.keys())

    for inchikey in correct_matched_dict.keys():
        correct_ms2_list = correct_matched_dict[inchikey].MS2
        generated_ms2_list = matched_dict[inchikey].MS2
        try:
            assert len(correct_ms2_list) == len(generated_ms2_list)
        except AssertionError:
            success = False
            print 'Failure. Incorrect number of MS2 objects generated.'
        for ms2_index in xrange(len(correct_ms2_list)):
            correct_ms2 = correct_ms2_list[ms2_index]
            generated_ms2 = generated_ms2_list[ms2_index]
            correct_attributes = [attribute for attribute in dir(correct_ms2) if '__' not in attribute]
            generated_attributes = [attribute for attribute in dir(generated_ms2) if '__' not in attribute]
            try:
                assert set(correct_attributes) == set(generated_attributes)
            except AssertionError:
                success = False
                print 'Failure. Generated MS2 attributes %s did not match the expected: %s' % (sorted(generated_attributes), sorted(correct_attributes))

            for attribute in correct_attributes:
                correct_value = getattr(correct_ms2, attribute)
                generated_value = getattr(generated_ms2, attribute)
                if attribute == 'peaks':
                    try:
                        if not correct_value:
                            assert not generated_value
                        if correct_value:
                            assert generated_value
                            assert len(correct_value) == len(generated_value)
                    except AssertionError:
                        success = False
                        print 'Failure. Generated peaks length does not match expected peaks length'
                    success_array = []
                    if correct_value and generated_value:
                        for peak_index in xrange(len(correct_value)):
                            success_array.append(compareMS2(correct_value[peak_index], generated_value[peak_index]))
                            try:
                                assert all(success_array)
                            except AssertionError:
                                success = False
                                print 'Failure. Peaks lists did not match.'

                else:
                    try:
                        assert correct_value == generated_value
                    except AssertionError:
                        success = False
                        print 'Failure. Generated %s value of %s did not matched the expected: %s' % (attribute, generated_value, correct_value)
    return success

def testWriteToCSV(testingUnpack=False):
    '''
    Tests writeToCSV()

    Tests all cases, compares to answer key CSV
    '''
    ms2_object = MS2()
    id_dict = dict()
    id_dict['database_id'] = 'testid'
    id_dict['id'] = '1'
    ms2_object.id_dict = id_dict
    ms2_object.energy_field = ''
    ms2_object.inchi_key = 'testinchikey'
    ms2_object.inchikey = 'testinchikey'
    references_dict = dict()
    references_dict['database'] = 'HMDB'
    references_dict['database_id'] = 'testid'
    references_dict['ref_text'] = ''
    ms2_object.references_dict = references_dict
    peak1 = MSPeak()
    peak1.id = '1'
    peak1.mass_charge = '1.1'
    peak1.molecule_id = ''
    peak2 = MSPeak()
    peak2.id = '2'
    ms2_object.peaks = [peak1, peak2]

    metabolite_dict = testMetabolitePreprocessing(True)
    for key in metabolite_dict.iterkeys():
        metabolite_dict[key].MS2 = [ms2_object]

    if testingUnpack: return metabolite_dict

    file_path = os.getcwd() + '/tests/test.csv'
    writeToCSV(metabolite_dict, file_path)

    generated_csv = open(file_path, 'r')
    generated_reader = csv.reader(generated_csv)
    answer_key_file_path = os.getcwd() + '/tests/answer_key.csv'
    answer_key = open(answer_key_file_path, 'r')
    answer_reader = csv.reader(answer_key)

    answer_rows = []
    for row in answer_reader:
        answer_rows.append(row)
    answer_key.close()

    row_count = 0
    for row in generated_reader:
        if row != answer_rows[row_count]:
            return False
        row_count += 1
    generated_csv.close()
    return True


def testUnpackCSV():
    '''
    Tests unpackCSV()

    '''
    correct_dict = testWriteToCSV(True)
    testWriteToCSV(False)
    csv_path = os.getcwd() + '/tests/test.csv'
    generated_dict = unpackCSV(csv_path)

    success = True
    if set(correct_dict.keys()) != set(generated_dict.keys()):
        print 'key mismatch'
        return False

    for key in correct_dict.iterkeys():
        correct_metabolite = correct_dict[key]
        generated_metabolite = generated_dict[key]
        if not compareMetabolites(correct_metabolite, generated_metabolite):
            return False
    return True

def testAll():
  '''
  Tests metabolitePreprocessing(), MS2Preprocessing(), writeToCSV(), and unpackCSV()
  with minimal, yet exhaustive test cases.

  Raises:
    AssertionError if not all test cases pass

  '''
  success = []
  success.append(testMetabolitePreprocessing())
  success.append(testMS2Preprocessing())
  success.append(testWriteToCSV())
  success.append(testUnpackCSV())
  assert all(success)

import argparse
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('ms2_concat', help='path to file with all MS2s '
        + 'concatenated into one xml file.')
    p.add_argument('metabolites_concat', help='path to file with all '
        + 'HMDB metabolites in one xml file.')
    p.add_argument('out', help='path to write output csv file to.')
    return p.parse_args()

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
  metabolite_xml_file = args.metabolites_concat
  metabolite_feature_set = set(['accession', 'secondary_accessions', 'name',
                'chemical_formula',
                'monisotopic_molecular_weight', 'iupac_name',
                'traditional_iupac', 'cas_registry_number',
                'smiles', 'inchi', 'inchikey',
                'biofluid_locations', 'taxonomy'])

  metabolite_dict = metabolitePreprocessing(metabolite_xml_file, metabolite_feature_set)
  matched_dict, failure_list = MS2Preprocessing(ms2_xml_file, ms2_feature_set, metabolite_dict)

  csv_file_path = args.out
  writeToCSV(matched_dict, csv_file_path)

  #######    TO READ IN CSV    #######
  # file_path = os.getcwd() + '/aggregated_metabolites.csv'
  # m = unpackCSV(file_path)


'''
74234 Metabolites
865 with MS2 objects

Open Questions:
-To use HMDB key or InChi key?
-How many molecules in HMDB (with metadata)
-Stats on chromotography types (preferably LC=liquid chromotography)
Of those molecules that have LC-MS/MS, what's the distribution of chemical classes in there
- <ionization-mode> {positive, negative}
-For molecules that have liquid chromotagraphy, how many have positive, negative, both modes?
'''
