'''

THIS FILE IS NOT YET DONE


'''


import re
import sys
import os
import time
import csv
sys.path.append('/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
from bs4 import BeautifulSoup
import requests
from unidecode import unidecode

sys.path.append(os.getcwd() + '/../util/')
import util


class MBMS2(object):

  def __init__(self):
    self.name = []
    self.formula = None
    self.iupac = None
    self.smiles = None
    self.inchikey = None
    self.MS2 = []



'''
Get the url to each individual text file of MS2 data
'''
def getLinks():
  url = 'http://www.massbank.jp/SVN/OpenData/record/'

  main_page = requests.get(url)
  soup = BeautifulSoup(main_page.content, 'html.parser')
  prefix = 'http://www.massbank.jp/SVN/OpenData/record/'
  all_links = [prefix+unidecode(link.get('href')) for link in soup.find_all('a')[1:-1]]

  link_queue = list(all_links)
  discovered_links = set()
  while link_queue:
    current_link = link_queue.pop(0)
    page = requests.get(current_link)
    if page.status_code != 200:
      print 'Bad Status Code'
    soup = BeautifulSoup(page.content, 'html.parser')
    child_links = []
    for link in soup.find_all('a'):
      clean_link = current_link + link.get('href')
      if '.txt' in clean_link:
        child_links.append(clean_link)
    for child in child_links:
      discovered_links.add(unidecode(child))
  return list(discovered_links)


'''
Takes in list of match tuples with match attributes and generates a MBMS2 object
'''
def setAttributes(match_tuples, groups):
  container = MBMS2()
  for match_tuple in match_tuples:
    for index, item in match_tuple:
      if item:
        tuple_index = index
    group = groups[tuple_index]
    captured_value = match_tuple[tuple_index]
    attr_value = getattr(container, group.lower())
    if attr_value:
      attr_value.append(captured_value)
    else:
      setattr(container, group, captured_value)
  return container



def generateMBMS2(data_links):
  MBMS2_list = []
  name_dict = dict()
  inchikey_dict = dict()
  formula_dict = dict()
  iupac_dict = dict()
  smiles_dict = dict()
  name_pattern = 'CH\$NAME:\s(?P<name>.*)\n'
  formula_pattern = 'CH\$FORMULA:\s(?P<formula>.*)\n'
  smiles_pattern = 'CH\$SMILES:\s(?P<smiles>.*)\n'
  inchikey_pattern = 'INCHIKEY\s(?P<inchikey>.*)\n'
  iupac_pattern = 'CH\$IUPAC:\s(?P<iupac>.*)\n'
  match_tuple = (name_pattern, formula_pattern, smiles_pattern, inchikey_pattern, iupac_pattern)
  regex_string = '%s|%s|%s|%s|%s' % match_tuple
  regex = re.compile(regex_string)
  groups = ('name', 'pattern', 'smiles', 'inchikey', 'iupac')

  MS2_regex_string = '\s\s(?P<mz_intensity>\d*\.?\d*\s\d*\.?\d*)\s(?:\d*\.?\d*)\n'
  MS2_regex = re.compile(MS2_regex_string)
  
  for link in data_links:
    content = requests.get(link).content
    match = regex.findall(content)
    container = setAttributes(match, groups)
    ms2_match = MS2.findall(content)
    ms2_list = []
    for match in ms2_match:
      ms2_list.append(tuple(match.split(' ')))
    container.MS2 = ms2_list
    MBMS2_list.append(container)
    name_dict[container.name] = container
    formula_dict[container.formula] = container
    iupac_dict[container.iupac] = container
    smiles_dict[container.smiles] = container
    inchikey_dict[container.inchikey] = container
  return MBMS2_list, name_dict, formula_dict, iupac_dict, smiles_dict, inchikey_dict




csv_path = os.getcwd() + '/'
metabolite_dict = util.unpackCSV(csv_path)
for metabolite in metabolite_dict.itervalues():
  if metabolite.


# def getInchis(discovered_links):
#   inchi_list = []
#   # print len(discovered_links), 'total discovered links'
#   for index, link in enumerate(discovered_links[4000:]):
#     try:
#       print float(index)/len(discovered_links[4000:]), 'percent done'
#       content = requests.get(link).content
#       regex = r'INCHIKEY \s (.*) \n'
#       match = re.search(regex, content)
#       if match:
#         inchi_list.append(match.group(0))
#     except requests.exceptions.ConnectionError:
#       print 'hit error on index %d, sleeping' % index
#       time.sleep(3)
#       inchi_list.extend(getInchis(discovered_links[index:]))
#       return inchi_list
#   return inchi_list

# inchi_set = set(getInchis(discovered_links))

# print 'finished getting inchis from MassBank, now getting from HMDB'


# metabolite_inchis = set()
# file_path = os.getcwd() + '/../../../metabolomics/aggregated_metabolites.csv'
# csv_file = open(file_path, 'r')
# reader = csv.reader(csv_file)
# for row in reader:
#   if row[0]:
#     metabolite_inchis.add(row[0])

# print 'starting here'
# print metabolite_inchis.issuperset(inchi_set)
# print len(inchi_set)
# print len(metabolite_inchis - inchi_set)
# print len(inchi_set - metabolite_inchis)

# print '\a'
# print '\a'
# print '\a'

if __name__ == '__main__':
  pass
  # data_links = getLinks()






