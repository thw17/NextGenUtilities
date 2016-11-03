import argparse
import subprocess
import sys


def main():
	""" Main Function """
	args = parse_args()
	print args

	if args.wrap_length is not None:
		wrap = int(args.wrap_length)

	# Check bioawk install
	bioawk = args.bioawk
	a = subprocess.call(
		"{} -c help".format(bioawk), shell=True)
	# The -c help command returns 1 as an exit code
	if a != 1:
		print "Error. Check bioawk installation and path"
		sys.exit(1)

	# Create temp fasta for large enough contigs/scaffolds
	a = subprocess.call(
		"""{} -c fastx '{{ if (length($seq) > {} ) print ">"$name; print $seq}}' {} > tmp_large_enough.fa""".format(
			bioawk, int(args.size) - 1, args.fasta), shell=True)
	# Create a file of just sequence for too short scaffolds/contigs
	a = subprocess.call(
		"""{} -c fastx '{{ if (length($seq) < {} ) print $seq}}' {} > tmp_toosmall_seq.txt""".format(
			bioawk, args.size, args.fasta), shell=True)
	# Create a file of just sequence names for too short scaffolds/contigs
	a = subprocess.call(
		"""{} -c fastx '{{ if (length($seq) < {} ) print $name}}' {} > tmp_toosmall_name.txt""".format(
			bioawk, args.size, args.fasta), shell=True)
	# Concatenate too short sequences
	with open("tmp_toosmall_seq.txt", "r") as f:
		sequence = ""
		lengths = []
		for line in f:
			sequence += line
			lengths.append(len(line))
	# Collect sequence names
	with open("tmp_toosmall_name.txt", "r") as f:
		names = []
		for line in f:
			names.append(line)
	# Write fasta of supercontig
	length = len(sequence)
	with open("tmp_supercontig.fa", "w") as f:
		f.write(">{}\n".format(args.supercontig_name))
		if wrap is not None:
			for i in range(0, length, wrap):
				f.write(sequence[i: i + wrap] + "\n")
		else:
			f.write(sequence + "\n")
	# Write bed containing coordinates of contigs in the supercontig
	with open(args.output_bed, "w") as f:
		start = 0
		stop = 0
		for idx, i in lengths:
			stop += i
			f.write("{}\t{}\t{}\t{}".format(args.supercontig_name, start, stop, names[idx]))
			start = stop
	# Concatenate fasta of sequences large enough with supercontig fasta
	a = subprocess.call(
		"cat tmp_large_enough.fa tmp_supercontig.fa > {}".format(
			args.output_fasta), shell=True)

def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		description="Concatenates all contigs smaller than a given length")

	parser.add_argument(
		"--fasta", required=True,
		help="Fasta reference file to be processed")

	parser.add_argument(
		"--output_fasta", required=True,
		help="Output file to write new reference")

	parser.add_argument(
		"--output_bed", required=True,
		help="Output bed file for annotations of new fasta")

	parser.add_argument(
		"--size", type=int, default=1000,
		help="Minimum length for contig/scaffold to be included. Default = 1000")

	parser.add_argument(
		"--supercontig_name", default="chrUn",
		help="The name of the supercontig created from the filtered contigs. "
		"Default is chrUn")

	parser.add_argument(
		"--wrap_length", default=50,
		help="Provide either intger length of each line to use "
		"when writing sequence data or None to prevent wrapping and write "
		"sequence on a single line. Default = 50.")

	parser.add_argument(
		"--bioawk", default="bioawk",
		help="Path to bioawk.  Default is 'bioawk'")

	args = parser.parse_args()

	return args


if __name__ == '__main__':
	main()
