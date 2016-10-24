import argparse
from Bio import SeqIO
from Bio.seq import Seq


def main():
	""" Main function """

	args = parse_args()

	output_bed_collector = []
	chrUn = Seq("", generic_dna)
	start = 0
	stop = 0
	for seq_record in SeqIO.parse(args.fasta, "fasta"):
		length = len(seq_record)
		if length < args.size:
			chrUn += seq_record.seq
			temp_line = 


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
		help="Minimum length for contig/scaffold to be included")

	parser.add_argument(
		"--supercontig_name", default="chrUn",
		help="The name of the supercontig created from the filtered contigs")

	args = parser.parse_args()

	return args
