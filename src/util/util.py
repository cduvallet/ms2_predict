#!/usr/bin/env
"""
This file contains useful functions used multiple times throughout this project.
"""
from MetabolomicsObjects import Metabolite, MS2, MSPeak
import csv

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
        # Dictionary in cell should have been written as key:value;;key2:value2
        attribute = cell.split('=')[0]
        if ':' in cell:
          attribute_dictionary = dict()
          dictionary = cell.split('=', 1)[1]
          split_dictionary = dictionary.split(';;')
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
          cell_list = cell.split('=', 1)[1]
          split_list = cell_list.split(',')
          setattr(metabolite, attribute, split_list)
        #sole (attribute, value) pair in cell but not just cell w/ only inchikey
        else:
          if '=' not in cell:
            setattr(metabolite, 'inchikey', cell)
            setattr(metabolite, 'inchi_key', cell)
          else:
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
        attribute = cell.split('=')[0]
        if ':' in cell:
          attribute_dictionary = dict()
          dictionary = cell.split('=', 1)[1]
          split_dictionary = dictionary.split(';;')
          for item in split_dictionary:
            key = item.split(':')[0]
            if len(item.split(':')) == 2:
              value = item.split(':')[1]
            else:
              value = None
            attribute_dictionary[key] = value
          setattr(ms2_object, attribute, attribute_dictionary)
        # #ms_peak data in cell
        # elif 'ms_peak' in cell and '!' in cell:
        #   peak_list = cell.split('=')[1].split(',')
        #   for peak in peak_list:
        #     peak_attributes = peak_list.split('!')
        #     peak_object = MSPeak()
        #     for peak_item in peak_items:
        #       peak_attribute = peak_items.split('=')[0]
        #       peak_value = peak_items.split('=')[1]
        #       setattr(peak_object, peak_attribute, peak_value)

        # The "peaks" attributes are separated by exclamation points
        elif '!' in cell and attribute == "peaks":
            mspeaks_str = cell.split('=', 1)[1]
            # This creates a list with ['peaks=', 'id=123',
            # 'intensity=123', ...]
            mspeaks = mspeaks_str.split('!')
            # Collect all values corresponding to 'mass_charge=' key
            mzs = [i.split('=')[1] for i in mspeaks
                   if i.startswith('mass_charge')]
            intensities = [i.split('=')[1] for i in mspeaks
                           if i.startswith('intensity')]
            # Store MS2 peaks as a list of tuples
            mz_i = [(float(i), float(j)) for i, j in zip(mzs, intensities)]
            ms2_object.peaks = mz_i

        #sole attribute data in cell
        else:
          peak_attribute = cell.split('=')[0]
          peak_value = cell.split('=')[1]
          setattr(ms2_object, peak_attribute, peak_value)
      ms2_object_lst.append(ms2_object)
    metabolite.MS2 = ms2_object_lst
    metabolite_dict[metabolite.inchikey] = metabolite

  return metabolite_dict
