import multiprocessing
import argparse
import os
import pickle

version = '0.2'


def saveVariableToPickle(variableToStore, outdir, prefix):
	pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
	with open(pickleFile, 'wb') as writer:
		pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
	with open(pickleFile, 'rb') as reader:
		variable = pickle.load(reader)
	return variable


# Check whether a fastq file have all the required fields
def checkFastqFile(fastq, outdir):
	number_reads_components = [0, 0, 0, 0]
	with open(fastq, 'rtU') as reader:
		plus_line = True
		quality_line = True
		length_sequence = 0
		for line in reader:
			if len(line) > 0:
				if line.startswith('@') and plus_line and quality_line:
					number_reads_components[0] += 1
					plus_line = False
					quality_line = False
					length_sequence = 0
				elif line.startswith('+') and not plus_line:
					number_reads_components[2] += 1
					plus_line = True
				elif plus_line and not quality_line:
					number_reads_components[3] += 1
					quality_line = True
					line = line.splitlines()[0]
					if len(line) != length_sequence:
						print 'Sequence length and quality length are not equal!'
						break
				else:
					number_reads_components[1] += 1
					line = line.splitlines()[0]
					length_sequence = len(line)
	print fastq + ' -> ' + str(number_reads_components)
	saveVariableToPickle(outdir, [os.path.basename(fastq), number_reads_components[0], number_reads_components[0] == number_reads_components[1] == number_reads_components[2] == number_reads_components[3]], os.path.basename(fastq))


def runCheckFastq(args):
	threads = args.threads
	outdir = os.path.abspath(args.outdir)
	inputFastqFiles = args.inputFastqFiles
	for i in range(0, len(inputFastqFiles)):
		inputFastqFiles[i] = inputFastqFiles[i].name
	pool = multiprocessing.Pool(processes=threads)
	for fastq in inputFastqFiles:
		pool.apply_async(checkFastqFile, args=(fastq, outdir,))
	pool.close()
	pool.join()

	with open(os.path.join(outdir, 'report.number_reads.tab'), 'wt') as writer:
		writer.write('#file' + '\t' + 'numberReads' + '\t' + 'fastq_well_formatted' + '\n')
		files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
		for file_found in files:
			if file_found.endswith('.pkl'):
				file_path = os.path.join(outdir, file_found)
				sampleINFO = extractVariableFromPickle(file_path)
				writer.write('\t'.join(sampleINFO) + '\n')
				os.remove(file_path)


def main():

	parser = argparse.ArgumentParser(prog='python checkFastqFiles.py', description="Check whether fastq files are well formatted", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser_required = parser.add_argument_group('Required options')
	parser_required.add_argument('-i', '--inputFastqFiles', nargs='+', type=argparse.FileType('r'), metavar='/path/to/reference/genome/fastq/file.fa', help='Path to FASTQ file or files (if more than one are provided, they must be separated by space)', required=True)
	parser_optional = parser.add_argument_group('Facultative options')
	parser_optional.add_argument('-o', '--outdir', type=str, metavar='/output/directory/', help='Path for output directory', required=False, default='.')
	parser_optional.add_argument('-j', '--threads', metavar=('N'), type=int, help='Number of threads to be used', required=False, default=1)

	parser.set_defaults(func=runCheckFastq)

	args = parser.parse_args()

	args.func(args)


if __name__ == "__main__":
	main()
