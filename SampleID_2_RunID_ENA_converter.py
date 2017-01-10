#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
SampleID_2_RunID_ENA_converter.py - Convert ENA/SRA Sample accession ID or
Secondary sample accession ID to ENA/SRA Run accession ID (via ENA repository)
<https://github.com/B-UMMI/getSeqENA/>

Copyright (C) 2016 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: September 25, 2016

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import urllib2
import xml.etree.ElementTree as ET
import argparse
import pickle
import multiprocessing
import sys
import os


version = '0.2'


def check_create_directory(directory):
	if not os.path.isdir(directory):
		os.makedirs(directory)


def saveVariableToPickle(variableToStore, outdir, prefix):
	pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
	with open(pickleFile, 'wb') as writer:
		pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
	with open(pickleFile, 'rb') as reader:
		variable = pickle.load(reader)
	return variable


def sampleID_2_RunID(sampleID):
	sample_info = {'sample_primary_ID': None, 'sample_secondary_ID': None, 'ena_run': None, 'ena_study': None, 'center_name': None, 'attributes': {}}

	url = "http://www.ebi.ac.uk/ena/data/view/" + sampleID + "&display=xml"

	try:
		content = urllib2.urlopen(url)
		xml = content.read()
		tree = ET.fromstring(xml)
	except:
		print 'It was not possible to connect ENA for ' + sampleID + ' sample ID'
	else:
		for child_1 in tree:
			if child_1.tag == 'SAMPLE':
					sample_info['center_name'] = child_1.attrib['center_name']

			for child_2 in child_1:
				if child_2.tag == 'IDENTIFIERS':
					for child_3 in child_2:
						if child_3.tag == 'PRIMARY_ID':
							sample_info['sample_primary_ID'] = child_3.text
						elif child_3.tag == 'EXTERNAL_ID':
							sample_info['sample_secondary_ID'] = child_3.text
				elif child_2.tag == 'SAMPLE_LINKS':
					for child_3 in child_2:
						for child_4 in child_3:
							tag_DB = False
							text_ENA_RUN = False
							text_ENA_STUDY = False
							for child_5 in child_4:
								if child_5.tag == 'DB':
									tag_DB = True

								if child_5.text == 'ENA-RUN':
									text_ENA_RUN = True
								elif child_5.text == 'ENA-STUDY':
									text_ENA_STUDY = True

								if tag_DB and text_ENA_RUN and child_5.tag == 'ID':
									if sample_info['ena_run'] is None:
										sample_info['ena_run'] = child_5.text
										tag_DB = False
										text_ENA_RUN = False
									else:
										sys.exit('SampleID with more than one RunID!')
								elif tag_DB and text_ENA_STUDY and child_5.tag == 'ID':
									sample_info['ena_study'] = child_5.text
									tag_DB = False
									text_ENA_STUDY = False
				elif child_2.tag == 'SAMPLE_ATTRIBUTES':
					for child_3 in child_2:
						tag_text = None
						for child_4 in child_3:
							if child_4.tag == 'TAG':
								tag_text = child_4.text

							if tag_text is not None and child_4.tag == 'VALUE':
								sample_info['attributes'][tag_text.replace(' ', '_')] = child_4.text
								tag_text = None
		if sample_info['ena_run'] is None:
			print 'It was not possible to retrieve ENA for ' + sampleID + ' sample ID'
		else:
			print sampleID + ' - DONE'
	return sample_info


def get_sample_info(sampleID, outdir):
	sample_info = sampleID_2_RunID(sampleID)
	saveVariableToPickle(sample_info, outdir, str(sampleID + '_sample_info'))


def gather_all_samples_info(outdir):
	info = {}

	attributes = set([])

	with open(os.path.join(outdir, 'sampleID_with_problems.txt'), 'wt') as writer:
		files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
		for file_found in files:
			if file_found.endswith('_sample_info.pkl'):
				file_path = os.path.join(outdir, file_found)

				sample = file_found.split('_', 1)[0]
				sample_info = extractVariableFromPickle(file_path)

				if sample_info['ena_run'] is not None:
					info[sample] = sample_info
					attributes = attributes.union(set(sample_info['attributes'].keys()))
				else:
					writer.write(sample + '\n')

				os.remove(file_path)

	attributes = sorted(list(attributes))

	return info, attributes


def check_attributes_present(list_all_attributes, dict_sample_attributes):
	dict_all_attributes = {}
	for attribute in list_all_attributes:
		dict_all_attributes[attribute] = 'NA'

	for attribute in dict_sample_attributes:
		dict_all_attributes[attribute] = dict_sample_attributes[attribute]

	return dict_all_attributes


def run_SampleID_2_RunID_ENA_converter(args):
	threads = args.threads

	outdir = os.path.abspath(args.outdir)
	check_create_directory(outdir)

	inputSampleIDlist = args.inputSampleIDlist.name

	with open(inputSampleIDlist, 'rtU') as reader:
		inputSampleIDlist = []
		for line in reader:
			if len(line) > 0:
				inputSampleIDlist.append(line.splitlines()[0])

	pool = multiprocessing.Pool(processes=threads)
	for sample in inputSampleIDlist:
		pool.apply_async(get_sample_info, args=(sample, outdir,))
	pool.close()
	pool.join()

	samples_info, samples_attributes = gather_all_samples_info(outdir)

	counter = 0
	with open(os.path.join(outdir, 'sampleID_to_runID.tab'), 'wt') as writer:
		partial_header = ['sample_primary_ID', 'sample_secondary_ID', 'ena_run', 'ena_study', 'center_name']
		writer.write('#' + '\t'.join(partial_header) + '\t' + '\t'.join(samples_attributes) + '\n')
		for sample in samples_info:
			sample_values = []
			for field in partial_header:
				sample_values.append(str(samples_info[sample][field]))

			dict_sample_all_attributes = check_attributes_present(samples_attributes, samples_info[sample]['attributes'])
			for attribute in samples_attributes:
				sample_values.append(str(dict_sample_all_attributes[attribute]))

			writer.write('\t'.join(sample_values) + '\n')

			counter += 1

	if counter == 0:
		sys.exit('No SampleIDs  were RunID converted!')

	print '\n' + 'SampleID_2_RunID_ENA_converter.py is FINISHED'
	print '\n' + 'Many information was retrieved for ' + str(counter) + ' SampleIDs'
	print '\n' + 'To get a list with only RunIDs try run the following command:'
	print "cut -f 3 -d '\\t' " + os.path.join(outdir, 'sampleID_to_runID.tab') + " > " + os.path.join(outdir, 'sampleID_to_runID.only_runID.txt')


def main():
	parser = argparse.ArgumentParser(prog='python SampleID_2_RunID_ENA_converter.py', description="Convert ENA/SRA Sample accession ID or Secondary sample accession ID to ENA/SRA Run accession ID (via ENA repository)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser_required = parser.add_argument_group('Required options')
	parser_required.add_argument('-i', '--inputSampleIDlist', type=argparse.FileType('r'), metavar='/path/to/file/with/sampleID/list.txt', help='Path to file containing a list (one per line) of Sample accession ID or Secondary sample accession ID to be converted to Run accession IDs', required=True)

	parser_optional = parser.add_argument_group('Facultative options')
	parser_optional.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/', help='Path for output directory', required=False, default='./')
	parser_optional.add_argument('-j', '--threads', metavar=('N'), type=int, help='Number of threads to be used', required=False, default=1)

	parser.set_defaults(func=run_SampleID_2_RunID_ENA_converter)

	args = parser.parse_args()

	args.func(args)


if __name__ == "__main__":
	main()
