from __future__ import print_function
import argparse
from itertools import groupby
import subprocess
import sys


def main():
	""" Main Function """
	args = parse_args()

	if args.wrap_length is not None:
		wrap = int(args.wrap_length)

	# Check bioawk install
	bioawk = args.bioawk
	a = subprocess.call(
		"{} -c help".format(bioawk), shell=True)
	# The -c help command returns 1 as an exit code
	if a != 1:
		print("Error. Check bioawk installation and path")
		sys.exit(1)

	if args.mode == "LENGTH":
		# Create temp fasta for large enough contigs/scaffolds
		a = subprocess.call(
			"""{} -c fastx '{{ if (length($seq) > {} ) print ">"$name"\n"print $seq}}' {} > tmp_pass.fa""".format(
				bioawk, int(args.size) - 1, args.fasta), shell=True)

		# Create a file of just sequence for too short scaffolds/contigs
		a = subprocess.call(
			"""{} -c fastx '{{ if (length($seq) < {} ) print $seq}}' {} > tmp_toconcat_seq.txt""".format(
				bioawk, args.size, args.fasta), shell=True)

		# Create a file of just sequence names for too short scaffolds/contigs
		a = subprocess.call(
			"""{} -c fastx '{{ if (length($seq) < {} ) print $name}}' {} > tmp_toconcat_name.txt""".format(
				bioawk, args.size, args.fasta), shell=True)

	elif args.mode == "ID":
		# Code to read id text file
		with open(args.ids, "r") as f:
			id_list = [x.strip() for x in f]
			while "" in id_list:
				id_list.remove("")
		# base code on this kind of a line (this example prints yes if a sequence name is either chr1 or chr4):
		# bioawk -c fastx '{split("chr1 chr4",a, " "); for (i in a) value[a[i]]; if ($name in value == 1) print "Yes" }' test.fa

		# Create temp fasta for large enough contigs/scaffolds
		a = subprocess.call(
			"""{} -c fastx '{{split("{}", a, " "); for (i in a) value[a[i]]; if ($name in value == 1) print ">"$name{}$seq}}' {} > tmp_pass.fa""".format(
				bioawk, " ".join(id_list), repr('"\n"'), args.fasta), shell=True)
		print(
			"""{} -c fastx '{{split("{}", a, " "); for (i in a) value[a[i]]; if ($name in value == 1) print ">"$name{}print $seq}}' {} > tmp_pass.fa""".format(
			bioawk, " ".join(id_list), repr('"\n"'), args.fasta))

		# Create a file of just sequence for too short scaffolds/contigs
		a = subprocess.call(
			"""{} -c fastx '{{split("{}", a, " "); for (i in a) value[a[i]]; if ($name in value == 0) print $seq}}' {} > tmp_toconcat_seq.txt""".format(
				bioawk, " ".join(id_list), args.fasta), shell=True)

		# Create a file of just sequence names for too short scaffolds/contigs
		a = subprocess.call(
			"""{} -c fastx '{{split("{}", a, " "); for (i in a) value[a[i]]; if ($name in value == 0) print $name}}' {} > tmp_toconcat_name.txt""".format(
				bioawk, " ".join(id_list), args.fasta), shell=True)

	else:
		print("Please set --mode to either ID or LENGTH")
		sys.exit(1)

	# Concatenate too short sequences
	with open("tmp_toconcat_seq.txt", "r") as f:
		sequence = ""
		lengths = []
		for line in f:
			sequence += line.rstrip()
			lengths.append(len(line))

	# Collect sequence names
	with open("tmp_toconcat_name.txt", "r") as f:
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
		for idx, i in enumerate(lengths):
			stop += i
			f.write("{}\t{}\t{}\t{}".format(
				args.supercontig_name, start, stop, names[idx]))
			start = stop

	# Ensure large enough fasta sequences are wrapped correctly (as bioawk
	# doesn't wrap) and in a way consistent with the supercontig record
	with open("tmp_pass.fa", "r") as f:
		with open("tmp_pass.WRAPPED.fa", "w") as o:
			# the following text for parsing fasta based on Brent Pedersen's
			# response here:  https://www.biostars.org/p/710/
			faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
			for header in faiter:
				# drop the ">"
				header = header.next()[1:].strip()
				# join all sequence lines to one.
				seq = "".join(s.strip() for s in faiter.next())

				# Write output to file using same wrapping as supercontig above
				o.write(">{}\n".format(header))
				length = len(seq)
				if wrap is not None:
					for i in range(0, length, wrap):
						o.write(seq[i: i + wrap] + "\n")
				else:
					o.write(sequence + "\n")

	# Concatenate fasta of sequences large enough with supercontig fasta
	with open(args.output_fasta, "w") as f:
		a = subprocess.call(
			["cat", "tmp_pass.WRAPPED.fa", "tmp_supercontig.fa"],
			stdout=f)


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		description="Concatenates all contigs smaller than a given length")

	parser.add_argument(
		"--fasta", required=True,
		help="Fasta reference file to be processed")

	parser.add_argument(
		"--mode", required=True, choices=["LENGTH", "ID"],
		help="Mode to determine which contigs to leave intact, and which to "
		"combine. Two options: LENGTH and ID.  LENGTH will filter by contig "
		"length (via --size flag). ID will leave intact contigs with IDs "
		"included in a .txt file (via --ids flag)")

	parser.add_argument(
		"--output_fasta", required=True,
		help="Output file to write new reference")

	parser.add_argument(
		"--output_bed", required=True,
		help="Output bed file for annotations of new fasta")

	parser.add_argument(
		"--ids", default=None,
		help="Text file listing contig ids to leave as-is in a single column, one "
		"per row.  All other contig ids will be combined into the supercontig.")

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
