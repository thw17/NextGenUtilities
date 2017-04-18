"""
Remove_and_split_contamination_NCBI.py

Author: Tim Webster, Arizona State University
Date: April 11, 2017

Currently requires Python 2.7

This script takes as input a fasta file and a contamination file from NCBI. The
latter contains coordinates (1-based indexing) of contamination in the former.
This script will run through the fasta file and if a sequence does not contain
any contamination, then the script will simply print the sequence with a fasta
id. However, if a sequence does contain contamination, the script will remove
the contaminated sequence and split the fasta sequence at those locations.  For
example, a sequence with contamination from position 27-44 would be split into
two sequences: beginning to 27, and 44 until the end.

The script also has the option to strip repeated elements at the beginning of
fasta ids (e.g., strip "DNA|" from ">DNA|scaffold1", ">DNA|scaffold2", and so
on).  It also allows wrapping of fasta sequences at arbitrary lengths.
"""

from __future__ import print_function
import argparse
import csv
from itertools import groupby
import re


def main():
	""" Main Function """
	args = parse_args()

	if args.wrap_length != "None":
		wrap = int(args.wrap_length)
	else:
		wrap = None

	# Read NCBI contamination file
	# Assumes four columns looking something like:
	# scaffold1000	594793	294651..294688	adaptor:NGB00843.1
	# where the third column (e.g., 294651..294688) contains the one-based index
	# start and end position of the problem sequence to be removed, and commas
	# without spaces separate pairs of positions if more than one exist
	# (e.g., 294651..294688,300002..300036)
	cont_coords = {}
	with open(args.ncbi_tab, "r") as f:
		for line in csv.reader(f, delimiter="\t"):
			line2 = [x.rstrip() for x in line]
			if line2 != "":
				regex_pattern = '|'.join(map(re.escape, args.delimiters))
				coords = re.split(regex_pattern, line2[2])
				coords = [int(x) for x in coords]
				cont_coords[line2[0]] = coords

	# Process fasta
	f_strip = len(args.fasta_id_junk) + 1
	with open(args.fasta, "r") as f:
		with open(args.output_fasta, "w") as o:
			# the following text for parsing fasta based on Brent Pedersen's
			# response here:  https://www.biostars.org/p/710/
			faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
			for header in faiter:
				# drop the ">" and other junk
				header = header.next()[f_strip:].strip()
				header_split = header.split(" ")
				header_id = header_split[0]
				# join all sequence lines to one.
				seq = "".join(s.strip() for s in faiter.next())

				if header_id in cont_coords:
					suffix = 0
					temp_coord_list = cont_coords[header_id]
					# if len(temp_coord_list) > 2:
					for idx, co in enumerate(temp_coord_list[:-1]):
						if idx == 0:
							suffix += 1
							o.write(">{}\n".format(header_id + "_{}".format(suffix)))
							seq_1 = seq[0:co]
							length = len(seq_1)
							if wrap is not None:
								for i in range(0, length, wrap):
									o.write(seq_1[i: i + wrap] + "\n")
							else:
								o.write(seq_1 + "\n")
						elif idx % 2 != 0:
							suffix += 1
							o.write(">{}\n".format(header_id + "_{}".format(suffix)))
							seq_1 = seq[co: temp_coord_list[idx + 1]]
							length = len(seq_1)
							if wrap is not None:
								for i in range(0, length, wrap):
									o.write(seq_1[i: i + wrap] + "\n")
							else:
								o.write(seq_1 + "\n")
					# Note this else is associated with the for above NOT the if/elif
					else:
						suffix += 1
						o.write(">{}\n".format(header_id + "_{}".format(suffix)))
						seq_1 = seq[temp_coord_list[-1]:]
						length = len(seq_1)
						if wrap is not None:
							for i in range(0, length, wrap):
								o.write(seq_1[i: i + wrap] + "\n")
						else:
							o.write(seq_1 + "\n")

				else:
					o.write(">{}\n".format(header_id))
					length = len(seq)
					if wrap is not None:
						for i in range(0, length, wrap):
							o.write(seq[i: i + wrap] + "\n")
					else:
						o.write(seq + "\n")


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		description="This program takes information about fasta contamination "
		"from NCBI output and splits problem contigs in a fasta file, excluding "
		"the contaminated regions.  Unaffected scaffolds/contigs are left "
		"unchanged.")

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
		"--delimiters", nargs="*", default=["..", ","],
		help="The delimiters in the coordinate column of the NCBI contamination "
		"file.  Default are '..' and ','")

	parser.add_argument(
		"--fasta_id_junk", type=str, default="",
		help="Any junk sequence before the scaffold id on the fasta id lines, "
		"not including the '>'.  For example, if your fasta ids all look like "
		"'>lcl|scaffold1 Organism_name', where 'scaffold1' is the scaffold id, "
		"then you would include '--fasta_id_junk lcl|'.  This flag works by string "
		"length, not pattern matching.")

	args = parser.parse_args()

	return args

if __name__ == '__main__':
	main()
