import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
	""" Main function """

	args = parse_args()

	# Set up variables for loops
	counter = 0
	chrUn = Seq("", generic_dna)
	start = 0
	stop = 0
	if args.wrap_length is not None:
		wrap = int(args.wrap_length)
	with open(args.output_fasta, "w") as outfasta, open(args.output_bed, "w") as outbed:
		for seq_record in SeqIO.parse(args.fasta, "fasta"):
			try:
				length = len(seq_record)
			except TypeError:
				continue
			id = seq_record.id
			id_cleaned = id.replace("\n", " ").replace("\r", " ").replace("  ", " ")
			if length < 1:
				continue
			if length < args.size:
				chrUn += seq_record
				stop += length
				outbed.write("{}\n".format(
					"\t".join(
						[str(x) for x in [args.supercontig_name, start, stop, id_cleaned]])))
				start = stop
			else:
				print "\t".join([str(x) for x in [id_cleaned, 0, len(seq_record), id_cleaned]])
				outbed.write("{}\n".format(
					"\t".join(
						[str(x) for x in [id_cleaned, 0, len(seq_record), id_cleaned]])))
				outfasta.write(">{}\n".format(id_cleaned))
				if wrap is not None:
					for i in range(0, length, wrap):
						outfasta.write(str(seq_record.seq[i: i + wrap]) + "\n")
				else:
					outfasta.write(seq_record.seq + "\n")
			counter += 1
			if counter % 100 == 0:
				print "{} records processed".format(counter)
		chrUn_rec = SeqRecord(chrUn)
		chrUn_rec.id = args.supercontig_name
		length = len(SeqRecord(chrUn))
		outfasta.write(">{}\n".format(chrUn_rec.id))
		if wrap is not None:
			for i in range(0, length, wrap):
				outfasta.write(chrUn_rec[i: i + wrap].seq + "\n")
			else:
				outfasta.write(chrUn_rec.seq + "\n")


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

	parser.add_argument(
		"--wrap_length", default=50,
		help="Provide either intger length of each line to use "
		"when writing sequence data or None to prevent wrapping and write "
		"sequence on a single line")

	args = parser.parse_args()

	return args

if __name__ == '__main__':
	main()
