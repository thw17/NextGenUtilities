import argparse
from itertools import groupby
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
		for idx, i in enumerate(lengths):
			stop += i
			f.write("{}\t{}\t{}\t{}".format(
				args.supercontig_name, start, stop, names[idx]))
			start = stop

	# Ensure large enough fasta sequences are wrapped correctly (as bioawk
	# doesn't wrap) and in a way consistent with the supercontig record
	with open("tmp_large_enough.fa", "r") as f:
		with open("tmp_large_enough.WRAPPED.fa", "w") as o:
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
			["cat", "tmp_large_enough.WRAPPED.fa", "tmp_supercontig.fa"],
			stdout=f)

	# Remove blank lines (in place) from final output file
	a = subprocess.call(
		["sed", "-i '/^$/d'", args.output_fasta])


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

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence

	From Brent Pedersen: https://www.biostars.org/p/710/
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

if __name__ == '__main__':
	main()
