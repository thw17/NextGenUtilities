from __future__ import print_function
import argparse
import csv
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

	# Read NCBI contamination file
	# Assumes four columns looking something like:
	# scaffold1000	594793	294651..294688	adaptor:NGB00843.1
	# where the third column (e.g., 294651..294688) contains the one-based index
	# start and end position of the problem sequence to be removed
	cont_coords = {}
	with open(args.ncbi_tab, "r") as f:
		for line in csv.reader(f, delimiter="\t"):
			line2 = line.rstrip()
			if line2 != "":
				if "," in line2[2]:
					
				else:
					coords = line2[2].split("..")
					start = coords[0] - 1
					stop = coords[1]
					cont_coords[coords[0]] = [start, stop]

	# Read through fasta, find scaffolds wild





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
		"--ncbi_tab", required=True,
		help="Tab-delimited output file from NCBI containing contamination "
		"results.  Assumes scaffold name is first (Python = 0th) column and "
		"one-based coordinates are in the third (Python = 2nd) column.  Also "
		"assumes coordinates are separated by the delimiter: ..")

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
