import xml.etree.ElementTree as ET
import numpy as np
import csv
import os
'''
A container for all MS2 metadata. Final MS2 objects may have empty or additional
fields, depending on any anomalies in the data.
'''
class MS2:
	def __init__(self):
		self.peak_counter = None
		self.chemical_shift_y = None
		self.mono_mass = None
		self.spectra_type = None
		self.inchi_key = None
		self.frequency = None
		self.instrument_type = None
		self.energy_field = None
		self.chemical_shift_reference = None
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
		self.chemical_shift_x = None
		self.spectra_assessment = None
		self.solvent = None
		self.chemical_shift = None
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


'''
Container for data about an ms-ms peak. Final MSPeak objects may have empty or 
additional fields, depending on any anomalies in the data.
'''
class MSPeak:
	def __init__(self):
		self.id = None
		self.ms_ms_id = None
		self.intensity = None
		self.mass_charge = None


'''
Container for metadata about a metabolite. Final Metabolite objects may have
empty or additional fields, depending on any anomalies in the data.
'''
class Metabolite:
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
		self.MS2 = None



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


def metabolitePreprocessing(xml_file):
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

	metabolite_feature_set = set(['accession', 'secondary_accessions', 'name', 
							   'chemical_formula', 
							   'monisotopic_molecular_weight', 'iupac_name', 
							   'traditional_iupac', 'cas_registry_number',
							   'smiles', 'inchi', 'inchikey',
							   'biofluid_locations', 'taxonomy'])

	print 'Generating Metabolite XML Tree...'
	tree = ET.parse(xml_file)
	print 'Done. \nPreprocessing Metabolite Data...'
	root = tree.getroot()
	metabolite_dict = dict()
	for metabolite_tag in root:
		metabolite = Metabolite()
		id_dict = dict()
		for xml_object in metabolite_tag:
			clean_tag = xml_object.tag.replace('{http://www.hmdb.ca}', '')
			if not (clean_tag in metabolite_feature_set or '_id' in clean_tag):
				continue
			if clean_tag == 'secondary_accessions':
				secondary_accessions = []
				for accessions_object in xml_object:
					secondary_accessions.append(accessions_object.text)
				metabolite.secondary_accessions = secondary_accessions
			elif clean_tag == 'biofluid_locations':
				locations = []
				for location_object in xml_object:
					locations.append(location_object.text)
				metabolite.biofluid_locations = locations

			elif '_id' in clean_tag:
				id_dict[clean_tag] = xml_object.text

			elif clean_tag == 'taxonomy':
				taxonomy_dict = dict()
				for taxonomy_object in xml_object:
					clean_taxonomy_tag = taxonomy_object.tag.replace('{http://www.hmdb.ca}', '')
					if clean_taxonomy_tag != 'substitutents' and 'parent' not in clean_taxonomy_tag:
						taxonomy_dict[clean_taxonomy_tag] = taxonomy_object.text
				metabolite.taxonomy_dict = taxonomy_dict

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
	metabolite_dict = dict()
	ms2_object = MS2()
	failure_list = []
	for metabolite in root:
		id_dict = dict()
		failure = False
		for feature in metabolite:
			if 'peaks' in feature.tag and feature.tag != 'ms-ms-peaks':
				failure = True
				break
			if feature.tag == 'chromatography-type' and feature.text != 'GC':
				failure = True
				break

			if '_id' in feature.tag or feature.tag == 'id':
				id_dict[feature.tag] = feature.text

			elif feature.tag == 'references':
				references_dict = dict()
				for references_feature in feature:
					references_dict[references_feature.tag] = references_feature.text

			elif feature.tag == 'ms-ms-peaks':
				peak_objects = []
				for peak_object in feature:
					ms_peak = MSPeak()
					for peak_object_attribute in peak_object:
						setattr(ms_peak, peak_object_attribute.tag, peak_object_attribute.text)
					peak_objects.append(ms_peak)
				ms2_object.peaks = peak_objects

			else:
				setattr(ms2_object, feature.tag.replace('-', '_'), feature.text)

		if failure: break
		if metabolite_dict.get(ms2_object.inchikey):
			metabolite_dict[ms2_object.inchi_key].MS2 = ms2_object
		else:
			failure_list.append(ms2_object)
		
	print 'Done.'
	return metabolite_dict, failure_list


def writeToCSV(matched_dict, file_path):
	'''
	Takes in output from MS2Preprocessing() and writes dictionary to CSV.
	Format is one of:
		attribute=value
		attribute=value1,value2,value3,...
		attribute=key1:value1,key2:value2,...

	Args: 
		matched_dict: dictionary in the form of: {inchikey: Metabolite}, where
			Metabolite object has a populated MS2 field. This is the output of 
			MS2Preprocessing().
		file_path: path to which the csv file will be written
	'''
	print "Writing to CSV..."
	csv_file = open(file_path, 'w+')
	writer = csv.writer(csv_file)
	for metabolite_item in matched_dict.iteritems():
		metabolite = metabolite_item[1]
		ms2_object = getattr(metabolite, 'MS2')
		for metabolite_object in [metabolite, ms2_object]:
			attributes = dir(metabolite_object)
			attributes.remove('__doc__')
			attributes.remove('__init__')
			attributes.remove('__module__')
			attributes.remove('MS2')
			if type(metabolite_object) == MS2:
				writestring_list = ['']
			else:
				writestring_list = [metabolite_item[0]]
			for attribute in attributes:
				value = getattr(metabolite_object, attribute)
				if type(value) == str:
					writestring = str(attribute) + '=' + str(value)
				elif type(value) == dict:
					writestring = str(attribute) + '='
					for item in attribute.iteritems():
						writestring += str(item[0]) + ':' + str(item[1]) + ','
					writestring = writestring[:-1]
				elif type(value) == list:
					writestring = str(attribute) + '='
					for list_item in value:
						writestring += str(list_item) + ','
					writestring = writestring[:-1]
				elif type(value) == MSPeak:
					writestring = 'MS2='
					ms2_attributes = dir(value)
					ms2_attributes.remove('__doc__')
					ms2_attributes.remove('__init__')
					ms2_attributes.remove('__module__')
					for ms_attribute in ms2_attributes:
						ms_value = getattr(value, ms_attribute)
						writestring += '!'
						writestring += str(ms_attribute) + '=' + str(ms_value)
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
		if row_count % 2 == 1:
			metabolite = Metabolite()
			for cell in row:
				#dictionary in cell
				if ':' in cell:
					attribute_dictionary = dict()
					attribute = cell.split('=')[0]
					dictionary = cell.split('=')[1]
					split_dictionary = dictionary.split(',')
					for item in split_dictionary:
						key = item.split(':')[0]
						value = item.split(':')[1]
						attribute_dictionary[key] = value
					setattr(metabolite, attribute, attribute_dictionary)
				#list in cell
				elif ',' in cell:
					attribute = cell.split('=')[0]
					cell_list = cell.split('=')[1]
					split_list = cell_list.split(',')
					setattr(metabolite, attribute, split_list)
				#sole attribute, value pair in cell
				else:
					attribute = cell.split('=')[0]
					value = cell.split('=')[1]
					setattr(metabolite, attribute, value)
		#Row is for an MS2 object
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
						value = item.split(':')[1]
						attribute_dictionary[key] = value
					setattr(ms2_object, attribute, attribute_dictionary)
				#ms_peak data in cell
				elif 'ms_peak' in cell:
					peak_list = cell.split('=')[1]
					peak_attributes = peak_list.slit('!')
					peak_object = MSPeak()
					for peak_item in peak_items:
						peak_attribute = peak_items.split('=')[0]
						peak_value = peak_items.split('=')[1]
						setattr(peak_object, peak_attribute, peak_value)
				#sole attribute data in cell
				else:
					peak_attribute = cell.split('=')[0]
					peak_value = cell.split('=')[1]
					setattr(peak_object, peak_attribute, peak_value)
			metabolite.MS2 = ms2_object
			metabolite_dict[metabolite.inchikey] = metabolite



if __name__ == '__main__':
	feature_set = set(['inchi_key','frequency','instrument_type','id','energy_field','chemical_shift_reference','base_peak','sample_mass_units','chromatography_type','searchable','sample_mass','derivative_mw','retention_time','updated_at','sample_assessment','derivative_formula','derivative_type','database_id','ref_text','mass_charge','collision_energy_voltage','sample_concentration','chemical_shift_x','spectra_assessment','solvent','chemical_shift','nucleus_y','ionization_mode','collection_date','nucleus','sample_temperature_units','sample_ph','spectra_id','c_ms_id','sample_concentration_units','sample_source','nil_classes','nucleus_x','database','notes','created_at','sample_temperature','ri_type','pubmed_id','molecule_id','column_type','retention_index','collision_energy_level','references','name','accession','chemical_formula','monoisotopic_molecular_weight','iupac_name','traditional_iupac','cas_registry','smiles','inchi','inchikey','taxonomy','biofluid_locations','ids','peak_counter', 'chemical_shift_y', 'mono_mass', 'ms_ms_id', 'spectra_type'])

	ms2_xml_file = '/Users/julianalverio/Documents/CURRENT/6.867/metabolomics/hmdb_spectra_xml/combined_file_clean.xml'
	metabolite_xml_file = '/Users/julianalverio/Documents/CURRENT/6.867/metabolomics/hmdb_metabolites.xml'

	metabolite_dict = metabolitePreprocessing(metabolite_xml_file)

	matched_dict, failure_list = MS2Preprocessing(ms2_xml_file, feature_set, metabolite_dict)


	file_path = os.cwd() + '/aggregated_metabolites.csv'
	writeToCSV(matched_dict, file_path)


'''
Types of peaks:
'c-ms-peak', 'nmr-one-d-peak', 'c-ms-peaks', 'nmr-one-d-peaks', 'ms-ms-peaks', 'nmr-two-d-peaks', 'ei-ms-peak', 'nmr-two-d-peak'

Open Questions:
-To use HMDB key or InChi key?
-How many molecules in HMDB (with metadata)
-Stats on chromotography types (preferably LC=liquid chromotography)
Of those molecules that have LC-MS/MS, what's the distribution of chemical classes in there
- <ionization-mode> {positive, negative} 
-For molecules that have liquid chromotagraphy, how many have positive, negative, both modes?
'''

