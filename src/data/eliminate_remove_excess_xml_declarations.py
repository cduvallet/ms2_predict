import argparse
import os
import re

if __name__ == '__main__':
	inputxml = '/Users/julianalverio/Documents/CURRENT/6.867/Metabolomics\ Project/hmdb_spectra_xml/combined_file.xml'
	outputxml = '/Users/julianalverio/Documents/CURRENT/6.867/Metabolomics\ Project/hmdb_spectra_xml/combined_file_clean.xml'

	xml_top_tag = 'database'

	with open(inputxml, 'r+') as f, open(outputxml, 'w') as output:
		#Add a new root to the xml - you can only have one root
		for line_num, line in enumerate(f):
			#If the line contains an xml header, don't write it to
			#the new file, except the first tag
			if re.search('<\?xml version', line) and line_num == 0:
				#write out the xml declaration and the first tag
				output.write(line)
				output.write('<'+xml_top_tag+'>\n')
			elif re.search('<\?xml version', line):
				pass
			else:
				#keep the indentation the same as the original file
				output.write('  '+line)
				#write out something for user to see it's working.
				if line_num % 1000000 == 0:
					print "Writing line %s" % line_num
		#Close the root
		output.write('</'+xml_top_tag+'>')