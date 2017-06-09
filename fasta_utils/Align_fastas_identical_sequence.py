"""
Align_fastas_identical_sequence.py

Author: Tim Webster, Arizona State University
Date: June 8, 2017

Currently requires Python 2.7

This script will only work for cases in which the fasta sequences of one file
are *completely and exactly contained* in another.  For example, our need for
this came after we removed contamination, trimmed Ns, and split scaffolds from
an original fasta and we needed to convert coordinates between assemblies.  If
there are potentially gaps, substitutions, etc., use a different program
(e.g., lastz).

This script takes as input two fasta files, one "primary organism"
and one "aligning organism" (see:
https://genome.ucsc.edu/goldenPath/help/axt.html). In our use, the new assembly
was the aligning organism and the original assembly was the primary organism.

Note that only the only output format is axt.  For the blastz score column,
this script will calculate a score based on (Charomonte et al., 2002),
assuming 100% identity.  Because Charmonte et al. scale scores within an
alignment so that the highest scoring pair is 100, this effectively means we
multiply the length of the alignment by 100 for a score.
Alternatively, there is an option to add an identical string as the score
to every alignment (e.g., 3500).

"""
import argparse
from itertools import groupby


def main():
	""" Main function """
	args = parse_args()

	a_strip = len(args.aligning_fasta_id_junk) + 1
	p_strip = len(args.primary_fasta_id_junk) + 1

	alignment_number = 0
	with open(args.output_prefix + ".axt", "w") as out:
		with open(args.aligning_fasta, "r") as ao:
			# note: the following code for parsing fastas is largely based on
			# Brent Pedersen's response here:  https://www.biostars.org/p/710/
			ao_faiter = (x[1] for x in groupby(ao, lambda line: line[0] == ">"))
			for ao_header in ao_faiter:
				# drop the ">" and other junk
				ao_header_id = ao_header.next()[a_strip:].strip()
				if args.no_splitting_aligning_id is False:
					ao_header_id_split = ao_header_id.split()
					ao_header_id = ao_header_id_split[0]
				# join all sequence lines to one.
				ao_seq = "".join(s.strip() for s in ao_faiter.next())
				ao_seq_length = len(ao_seq)

				# run through primary fasta to find match
				with open(args.primary_fasta, "r") as po:
					po_faiter = (x[1] for x in groupby(po, lambda line: line[0] == ">"))
					for po_header in po_faiter:
						# drop the ">" and other junk
						po_header_id = po_header.next()[p_strip:].strip()
						if args.no_splitting_primary_id is False:
							po_header_id_split = po_header_id.split()
							po_header_id = po_header_id_split[0]
						# join all sequence lines to one.
						po_seq = "".join(s.strip() for s in po_faiter.next())
						po_seq_length = len(po_seq)

						if ao_seq_length <= po_seq_length:
							idx = po_seq.find(ao_seq)
							# -1 indicates string not found
							if idx != -1:
								aligned_po = po_seq[idx:idx + ao_seq_length]
								if args.fixed_score == "script_calculates":
									score = ao_seq_length * 100
								else:
									score = args.fixed_score
								out_header_list = [
									alignment_number,
									po_header_id,
									idx + 1,
									idx + ao_seq_length,
									ao_header_id,
									1,
									ao_seq_length,
									"+",
									score]
								out_header_list = [str(x) for x in out_header_list]
								out_header = " ".join(out_header_list)
								out.write("{}\n".format(out_header))
								out.write("{}\n".format(aligned_po))
								out.write("{}\n".format(ao_seq))
								out.write("\n")
								alignment_number += 1


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		description="This program outputs axt alignments for cases comparing "
		"fastas for which the sequences of one file are completely and exactly "
		" contained in another.  For example, our need for this came after "
		"we removed contamination, trimmed Ns, and split scaffolds from an original "
		"fasta and we needed to convert coordinates between assemblies.")

	parser.add_argument(
		"--aligning_fasta", required=True,
		help="Fasta file for which sequences are smaller pieces of "
		"primary_fasta")

	parser.add_argument(
		"--primary_fasta", required=True,
		help="Fasta file for which sequences completely contain sequences "
		"from aligning_fasta")

	parser.add_argument(
		"--output_prefix", required=True,
		help="The file prefix (including path) for output file.  .axt will "
		"automatically be added.  For example, '--output_prefix /home/out' will "
		"create the output file: /home/out.axt")

	parser.add_argument(
		"--fixed_score", default="script_calculates", type=str,
		help="Default is script_calculates.  If script_calculates, then score "
		"will be calulated using BLASTZ method. Else, the provided string will "
		"be used as the score for all alignments")

	parser.add_argument(
		"--primary_fasta_id_junk", type=str, default="",
		help="Any junk sequence before the scaffold id on the fasta id lines in "
		"the primary_fasta, not including the '>'.  For example, "
		"if your fasta ids all look like "
		"'>lcl|scaffold1 Organism_name', where 'scaffold1' is the scaffold id, "
		"then you would include '--fasta_id_junk lcl|'.  This flag works by string "
		"length, not pattern matching.")

	parser.add_argument(
		"--aligning_fasta_id_junk", type=str, default="",
		help="Any junk sequence before the scaffold id on the fasta id lines in "
		"the aligning_fasta, not including the '>'.  For example, "
		"if your fasta ids all look like "
		"'>lcl|scaffold1 Organism_name', where 'scaffold1' is the scaffold id, "
		"then you would include '--fasta_id_junk lcl|'.  This flag works by string "
		"length, not pattern matching")

	parser.add_argument(
		"--no_splitting_primary_id", action='store_true', default=False,
		help="If flag provided, will not split fasta sequence ids in primary "
		"fasta. Otherwise, ids split and 0th index item used as ID.")

	parser.add_argument(
		"--no_splitting_aligning_id", action='store_true', default=False,
		help="If flag provided, will not split fasta sequence ids in aligning "
		"fasta. Otherwise, ids split and 0th index item used as ID.")

	args = parser.parse_args()

	return args

if __name__ == "__main__":
	main()
