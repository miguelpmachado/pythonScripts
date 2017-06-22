#!/usr/bin/env python

# -*- coding: utf-8 -*-


"""
getCompleteGenomes.py - Get Bacterial Complete Genomes from NCBI
INNUca.py - INNUENDO quality control of reads, de novo assembly and contigs quality assessment, and possible contamination search
<https://github.com/miguelpmachado/pythonScripts>

Copyright (C) 2016 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: August 05, 2016

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


import multiprocessing
import argparse
import shlex
import subprocess
from threading import Timer
import time
import os.path
import pickle
import sys
import urllib

version = '0.1'


def runCommandPopenCommunicate(command, shell_True, timeout_sec_None):
	run_successfully = False
	if isinstance(command, basestring):
		command = shlex.split(command)
	else:
		command = shlex.split(' '.join(command))

	print 'Running: ' + ' '.join(command)
	if shell_True:
		command = ' '.join(command)
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	else:
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if timeout_sec_None is None:
		stdout, stderr = proc.communicate()
	else:
		timer = Timer(timeout_sec_None, proc.kill)
		timer.start()
		stdout, stderr = proc.communicate()
		timer.cancel()

	if proc.returncode == 0:
		run_successfully = True
	else:
		print 'STDOUT'
		print stdout.decode("utf-8")
		print 'STDERR'
		print stderr.decode("utf-8")
	return run_successfully, stdout, stderr


def check_create_directory(directory):
	if not os.path.isdir(directory):
		os.makedirs(directory)


def retreiveSpecies(ncbi_genome_summary, genus):
	list_species = []
	with open(ncbi_genome_summary, 'rtU') as reader:
		blank_line = True
		for line in reader:
			line = line.splitlines()[0]
			if len(line) == 0:
				blank_line = True
			else:
				if blank_line:
					blank_line = False
					number = None
					species = None
					try:
						number, species = line.split(' ', 1)
					except:
						continue
					try:
						number = int(number.split('.')[0])
					except:
						continue
					if number > 0:
						species = species.split(' ')
						if species[0] == genus:
							if len(species) >= 3:
								if species[1] != 'phage' and species[2] != 'phage':
									list_species.append(species)
							elif len(species) >= 2:
								if species[1] != 'phage':
									list_species.append(species)
	return list_species


def getListGenomesSpecies(species_list, outdir):
	url = ['http://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=download&orgn=', '', '[orgn]&status=50&report=proks&group=--%20All%20Prokaryotes%20--&subgroup=--%20All%20Prokaryotes%20--&format=']
	url[1] = '%20'.join(species_list)
	url = ''.join(url)
	command = ['wget', '-O', os.path.join(outdir, str('_'.join(species_list) + '.NCBI_genomes_proks.completeGenomes.' + time.strftime("%Y%m%d-%H%M%S") + '.tab')), str('"' + url + '"')]
	run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
	variableToPickle = {'_'.join(species_list): run_successfully}
	saveVariableToPickle(variableToPickle, outdir, str('_'.join(species_list) + '.NCBI_genomes_proks.completeGenomes'))


# Rename sequences
def renameSequences(inputFasta, outputFasta):
	with open(outputFasta, 'wt') as writer:
		with open(inputFasta, 'rtU') as reader:
			for line in reader:
				if len(line) > 0:
					if line.startswith('>'):
						accession = line[1:].split(' ')[0]
						gi = convert_accession_2_gi(accession)
						line = '>' + str('gi|' + gi + '|') + ' ' + line[1:]
						writer.write(line)
					else:
						writer.write(line)


def convert_accession_2_gi(accession_number):
	url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + accession_number + '&rettype=gi'
	gi = urllib.urlopen(url).read().splitlines()[0]
	if len(gi) == 0:
		gi = None
	return gi


def getGenomes(file_list_complete_genomes, outdir):
	with open(file_list_complete_genomes, 'rtU') as reader:
		for line in reader:
			line = line.splitlines()[0]
			if len(line) > 0:
				if not line.startswith('#'):
					line = line.split('\t')
					ftp = line[19]
					sample = ftp.rsplit('/', 1)[1]

					url = ftp + '/' + sample + '_genomic.fna.gz'
					command = ['wget', '-O', os.path.join(outdir, str(sample + '_genomic.fna.gz')), str('"' + url + '"')]
					run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
					variableToPickle = {sample: run_successfully}
					saveVariableToPickle(variableToPickle, outdir, str(sample + '.fna'))
					if run_successfully:
						command_gz = ['gunzip', '--keep', os.path.join(outdir, str(sample + '_genomic.fna.gz'))]
						run_successfully, stdout, stderr = runCommandPopenCommunicate(command_gz, False, None)
						if run_successfully:
							renameSequences(os.path.join(outdir, str(sample + '_genomic.fna')), os.path.join(outdir, str(sample + '_genomic.fna.renamed.fasta')))

					command[2] = os.path.join(outdir, str(sample + '_genomic.gbff.gz'))
					url = ftp + '/' + sample + '_genomic.gbff.gz'
					command[3] = '"' + url + '"'
					run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
					variableToPickle = {sample: run_successfully}
					saveVariableToPickle(variableToPickle, outdir, str(sample + '.gbff'))


def saveVariableToPickle(variableToStore, outdir, prefix):
	pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
	with open(pickleFile, 'wb') as writer:
		pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
	with open(pickleFile, 'rb') as reader:
		variable = pickle.load(reader)
	return variable


def runTime(start_time):
	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken, 3600)
	minutes, seconds = divmod(rest, 60)
	print 'Runtime :' + str(hours) + 'h:' + str(minutes) + 'm:' + str(round(seconds, 2)) + 's'
	return time_taken


def runGetCompleteGenomes(args):
	general_start_time = time.time()

	threads = args.threads[0]
	input_ncbi_genome_summary = os.path.abspath(args.input_ncbi_genome_summary[0].name)
	genus = args.genus[0]
	outdir = os.path.abspath(args.outdir[0])
	check_create_directory(outdir)

	list_species_inListFormat = retreiveSpecies(input_ncbi_genome_summary, genus)

	folder_files_list_genomes = os.path.join(outdir, 'complete_genomes_files_list', '')

	pool = multiprocessing.Pool(processes=threads)
	folder_genomes = os.path.join(outdir, 'complete_genomes_files', '')
	check_create_directory(folder_genomes)
	files = [f for f in os.listdir(folder_files_list_genomes) if not f.startswith('.') and os.path.isfile(os.path.join(folder_files_list_genomes, f))]
	for file_found in files:
		file_found = os.path.join(folder_files_list_genomes, file_found)
		pool.apply_async(getGenomes, args=(file_found, folder_genomes,))
	pool.close()
	pool.join()
	files = [f for f in os.listdir(folder_genomes) if not f.startswith('.') and os.path.isfile(os.path.join(folder_genomes, f))]
	with open(os.path.join(outdir, 'bad.complete_genomes_files.txt'), 'wt') as writer:
		for file_found in files:
			if file_found.endswith('.pkl'):
				file_found = os.path.join(folder_genomes, file_found)
				file_run_successfully = extractVariableFromPickle(file_found)
				for i in file_run_successfully:
					if not file_run_successfully[i]:
						writer.write(i + '\n')
						writer.flush()
				os.remove(file_found)

	print ''
	runTime(general_start_time)

	if os.path.getsize(os.path.join(outdir, 'bad.complete_genomes_files_list.txt')) > 0 or os.path.getsize(os.path.join(outdir, 'bad.complete_genomes_files.txt')) > 0:
		print ''
		sys.exit('Is was not possible to download some files!')


def main():

	parser = argparse.ArgumentParser(prog='python getCompleteGenomes.py', description="Get Bacterial Complete Genomes from NCBI", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser_required = parser.add_argument_group('Required options')
	parser_required.add_argument('-i', '--input_ncbi_genome_summary', nargs=1, type=argparse.FileType('r'), metavar='/path/to/file/with/ncbi/genome/summary.txt', help='Path to text file containing the NCBI genome summary from a genus', required=True)
	parser_required.add_argument('-g', '--genus', nargs=1, type=str, metavar='Streptococcus', help='The genus name to look for', required=True)
	parser_optional = parser.add_argument_group('Facultative options')
	parser_optional.add_argument('-o', '--outdir', nargs=1, type=str, metavar='/path/to/output/directory/', help='Path to where to store the outputs', required=False, default=['.'])
	parser_optional.add_argument('-j', '--threads', nargs=1, metavar=('N'), type=int, help='Number of threads to be used', required=False, default=[1])

	parser.set_defaults(func=runGetCompleteGenomes)

	args = parser.parse_args()

	args.func(args)



if __name__ == "__main__":
	main()
