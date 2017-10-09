import os

path_to_xml = os.getcwd() + '/hmdb_spectra_xml/'
filenames = os.listdir(path_to_xml)
combined_filename = path_to_xml + 'combined_file.xml'
with open(combined_filename, 'w') as outfile:
    for input_filename in filenames:
        with open(path_to_xml + input_filename) as infile:
            for line in infile:
                outfile.write(line)
