#!/usr/bin/env
"""
These functions were used during development, but are not incorporated or
used by any of the final code.
"""
import sys
sys.path.insert(0, '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

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
