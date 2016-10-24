from __future__ import print_function
import argparse
import gzip
import re


def main():
	""" Main Program """

	args = parse_args()
	input_fastq = args.fastq
	logfile = args.logfile
	tiles = args.tiles

	try:
		with open(tiles, "r") as f:
			lines = [line.split() for line in f]
			filters = [x.strip() for k in lines for x in k]
			while "" in filters:
				filters.remove("")
	except IOError:
		filters = tiles

	counter = 0
	filtered = 0
	if input_fastq[-2:] == "gz":
		with gzip.open(input_fastq, "rb") as f:
			while True:
				id = f.readline()
				seq = f.readline()
				third = f.readline()
				qual = f.readline()
				if len(id) > 0:
					if tile_filter(id, filters) is False:
						print(id.rstrip())
						print(seq.rstrip())
						print(third.rstrip())
						print(qual.rstrip())
						counter += 1
					else:
						counter += 1
						filtered += 1
				else:
					break
	else:
		with open(input_fastq, "r") as f:
			while True:
				id = f.readline()
				seq = f.readline()
				third = f.readline()
				qual = f.readline()
				if len(id) > 0:
					if tile_filter(id, filters) is False:
						print(id.rstrip())
						print(seq.rstrip())
						print(third.rstrip())
						print(qual.rstrip())
						counter += 1
					else:
						counter += 1
						filtered += 1
				else:
					break
	with open(logfile, "w") as lf:
		lf.write("Processed fastq file: %s\n" % input_fastq)
		lf.write("%d total reads\n" % counter)
		lf.write("%d filtered reads\n" % filtered)
		lf.write("%d reads remaining\n" % (counter - filtered))


def parse_args():
	parser = argparse.ArgumentParser(
		description="Filters fastq file to remove reads from a list of tiles")

	parser.add_argument(
		"--fastq", default="-",
		help="Fastq file to be filtered")

	parser.add_argument(
		"--logfile", default="fastq_log.txt",
		help="File to write output log.  Will overwrite if already exists.")

	parser.add_argument(
		"--tiles", required=True,
		help="List of tiles to be removed. Can be space-separated on the "
		"command line OR a text file listing tiles space-separated on a "
		"single line or in a single horizontal column with one tile per line.")

	args = parser.parse_args()

	return args


def tile_filter(sequence_id, filter_list):
	""" Takes a sequence identifier - first line of the four lines for each
		sequence in a FASTQ file (begins with a @) - and a list of tiles to be
		removed, and returns and returns True if the tile is in the list
		and False if the tile is not in the list.  Assumes Casava 1.8 and later
		format, where the tile is the fifth item (the 6th index)
		in the identifier.
		"""

	parsed = re.split('[ :]+', sequence_id)
	if parsed[4] in filter_list:
		return True
	else:
		return False


if __name__ == "__main__":
	main()
