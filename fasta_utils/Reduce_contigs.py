import argparse
import csv
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord


def main():
	""" Main function """

	args = parse_args()

	output_bed_collector = []
	chrUn = Seq("", generic_dna)
	start = 0
	stop = 0
	with open(args.output_fasta, "w") as outfasta:
		fasta_out = FastaIO.FastaWriter(outfasta)
		for seq_record in SeqIO.parse(args.fasta, "fasta"):
			length = len(seq_record)
			if length < args.size:
				chrUn += seq_record.seq
				stop += length
				output_bed_collector.append(
					[args.supercontig_name, start, stop, seq_record.id])
				start = stop
			else:
				output_bed_collector.append(
					[seq_record.id, 0, len(seq_record), seq_record.id])
				fasta_out.write_record(seq_record)
		chrUn_rec = SeqRecord(chrUn)
		chrUn_rec.id = args.supercontig_name
		fasta_out.write_record(chrUn_rec)

	# Write output
	with open(args.output_bed, "w") as f:
		w = csv.writer(f, dialect="excel-tab")
		w.writerows(output_bed_collector)


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

if __name__ == '__main__':
	main()
