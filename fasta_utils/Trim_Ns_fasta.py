"""
Trim_Ns_Fasta.py

Author: Tim Webster, Arizona State University
Date: April 17, 2017

Currently requires Python 2.7

This script takes as input a fasta file and trims Ns from the beginning and
end.  It also allows for abitrarily sized "buffers" to ignore in this
procedure.  For example, if you set a buffer to 2 and a fasta sequence started
with "AGNNNNNNNNNTCCG", the "AG" and Ns would all be trimmed.

It also allows for wrapping of fasta sequences at arbitrary lengths.

"""

from __future__ import print_function
import argparse
import csv
from itertools import groupby
import re


def main():
	""" Main Function """
	args = parse_args()

	if args.wrap_length != 'None':
		wrap = int(args.wrap_length)
	else:
		wrap = None

	if args.hard_buffer > 0:
		no_buffer = False
		lead_soft_buffer = 0
		tail_soft_buffer = 0
		lead_hard_buffer = args.hard_buffer
		tail_hard_buffer = args.hard_buffer

	elif args.soft_buffer > 0:
		no_buffer = False
		lead_hard_buffer = 0
		tail_hard_buffer = 0
		lead_soft_buffer = args.soft_buffer
		tail_soft_buffer = args.soft_buffer

	else:
		lead_hard_buffer = args.lead_hard_buffer
		tail_hard_buffer = args.tail_hard_buffer
		lead_soft_buffer = args.lead_soft_buffer
		tail_soft_buffer = args.tail_soft_buffer

	if all(
		[lead_hard_buffer, lead_soft_buffer, tail_hard_buffer,
			tail_soft_buffer]) == 0:
				no_buffer = True
	else:
		no_buffer = False

	if args.min_n > 1:
		lead_min_n = args.min_n
		tail_min_n = args.min_n
	else:
		if args.lead_min_n > 1:
			lead_min_n = args.lead_min_n
		else:
			lead_min_n = 1
		if args.tail_min_n > 1:
			tail_min_n = args.tail_min_n
		else:
			tail_min_n = 1

	if all([lead_min_n, tail_min_n]) == 0:
		no_min_n = True
	else:
		no_min_n = False

	# Process fasta
	with open(args.fasta, "r") as f:
		with open(args.output_fasta, "w") as o:
			# the following text for parsing fasta based on Brent Pedersen's
			# response here:  https://www.biostars.org/p/710/
			faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
			for header in faiter:
				# drop the ">"
				header_id = header.next()[1:].strip()
				# join all sequence lines to one.
				seq = "".join(s.strip() for s in faiter.next())
				# trim sequence
				new_seq = trim_sequence(
					seq, lead_hard_buffer, tail_hard_buffer,lead_soft_buffer,
					tail_soft_buffer, lead_min_n, tail_min_n)

				# Write fasta record
				o.write(">{}\n".format(header_id))
				length = len(new_seq)
				if wrap is not None:
					for i in range(0, length, wrap):
						o.write(new_seq[i: i + wrap] + "\n")
				else:
					o.write(new_seq + "\n")


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser(
		description="This program trims Ns from the beginning and/or ends "
		"of fasta records.  It also allows for an arbitrarily long hard or soft "
		"buffer to handle 'orphan' bases before Ns at the beginning of a "
		"sequence or after Ns at the end of a sequence.")

	parser.add_argument(
		"--fasta", required=True,
		help="Fasta reference file to be processed")

	parser.add_argument(
		"--output_fasta", required=True,
		help="Output file to write new reference")

	parser.add_argument(
		"--wrap_length", default=50,
		help="Provide either intger length of each line to use "
		"when writing sequence data or None to prevent wrapping and write "
		"sequence on a single line. Default = 50.")

	parser.add_argument(
		"--hard_buffer", default=0, type=int,
		help="If set, length of sequence to remove from the beginning and end "
		"of a fasta sequence before looking for Ns. This has the effect "
		"of hard clipping a sequence *whether or not there are Ns*.")

	parser.add_argument(
		"--soft_buffer", default=0, type=int,
		help="If set, script will remove *up to/including this length* "
		"from the beginning and end of a fasta sequence before looking for Ns. "
		"This is the dynamic version of --hard_buffer, and will not clip "
		"if Ns aren't discovered.")

	parser.add_argument(
		"--lead_hard_buffer", default=0, type=int,
		help="If set, length of sequence to remove from the BEGINNING "
		"of a fasta sequence before looking for Ns.  This has the effect "
		"of hard clipping a sequence *whether or not there are Ns*.")

	parser.add_argument(
		"--lead_soft_buffer", default=0, type=int,
		help="If set, script will remove *up to/including this length* "
		"from the BEGINNING of a fasta sequence before looking for Ns. "
		"This is the dynamic version of --hard_buffer, and will not clip "
		"if Ns aren't discovered.")

	parser.add_argument(
		"--tail_hard_buffer", default=0, type=int,
		help="If set, length of sequence to remove from the END "
		"of a fasta sequence before looking for Ns.  This has the effect "
		"of hard clipping a sequence *whether or not there are Ns*.")

	parser.add_argument(
		"--tail_soft_buffer", default=0, type=int,
		help="If set, script will remove *up to/including this length* "
		"from the END of a fasta sequence before looking for Ns. "
		"This is the dynamic version of --hard_buffer, and will not clip "
		"if Ns aren't discovered.")

	parser.add_argument(
		"--min_n", default=0, type=int,
		help="The minimum number of consecutive Ns for a sequence to be trimmed. "
		"This number is applied separately at the beginning and end of a sequence. "
		"This count starts AFTER buffers are applied.  Trimming Ns starts from "
		"the beginning of a sequence first and then moves to the end. Values of "
		"0 and 1 are treated equivalently (trim even if only a single N present)")

	parser.add_argument(
		"--lead_min_n", default=0, type=int,
		help="The minimum number of consecutive Ns for a sequence to be trimmed. "
		"This number is applied to the BEGINNING of a sequence. "
		"This count starts AFTER buffers are applied.  Trimming Ns starts from "
		"the beginning of a sequence first and then moves to the end. Values of "
		"0 and 1 are treated equivalently (trim even if only a single N present)")

	parser.add_argument(
		"--tail_min_n", default=0, type=int,
		help="The minimum number of consecutive Ns for a sequence to be trimmed. "
		"This number is applied to the END of a sequence. "
		"This count starts AFTER buffers are applied.  Trimming Ns starts from "
		"the beginning of a sequence first and then moves to the end. Values of "
		"0 and 1 are treated equivalently (trim even if only a single N present)")

	args = parser.parse_args()

	# Validate command line flags
	if any(
		[args.hard_buffer, args.soft_buffer, args.lead_hard_buffer,
			args.lead_soft_buffer, args.tail_hard_buffer, args.tail_soft_buffer,
			args.min_n, args.lead_min_n, args.tail_min_n]) < 0:
		sys.exit(
			"Command line parameters cannot be negative for any buffer or "
			"min_n flags.")

	if args.hard_buffer > 0:
		if any(
			[args.soft_buffer, args.lead_hard_buffer, args.lead_soft_buffer,
				args.tail_hard_buffer, args.tail_soft_buffer]) > 0:
			sys.exit(
				"--hard_buffer set, so no other buffer flags can be set.")

	if args.soft_buffer > 0:
		if any(
			[args.hard_buffer, args.lead_hard_buffer, args.lead_soft_buffer,
				args.tail_hard_buffer, args.tail_soft_buffer]) > 0:
			sys.exit(
				"--soft_buffer set, so no other buffer flags can be set.")

	if args.lead_hard_buffer > 0:
		if any(
			[args.hard_buffer, args.soft_buffer, args.lead_soft_buffer]) > 0:
			sys.exit(
				"--lead_hard_buffer set, so --soft_buffer and --lead_soft_buffer "
				"cannot be set.")

	if args.lead_soft_buffer > 0:
		if any(
			[args.hard_buffer, args.soft_buffer, args.lead_hard_buffer]) > 0:
			sys.exit(
				"--lead_soft_buffer set, so --soft_buffer and --lead_hard_buffer "
				"cannot be set.")

	if args.tail_hard_buffer > 0:
		if any(
			[args.hard_buffer, args.soft_buffer, args.tail_soft_buffer]) > 0:
			sys.exit(
				"--tail_hard_buffer set, so --soft_buffer and --tail_soft_buffer "
				"cannot be set.")

	if args.tail_soft_buffer > 0:
		if any(
			[args.hard_buffer, args.soft_buffer, args.tail_hard_buffer]) > 0:
			sys.exit(
				"--tail_hard_buffer set, so --soft_buffer and --tail_hard_buffer "
				"cannot be set.")

	if args.min_n > 0:
		if any(
			[args.lead_min_n, args.tail_min_n]) > 0:
			sys.exit(
				"--min_n set, so --lead_min_n and tail_min_n cannot be set.")

	if args.lead_min_n > 0:
		if any(
			[args.min_n]) > 0:
			sys.exit(
				"--lead_min_n set, so --min_n cannot be set.")

	if args.tail_min_n > 0:
		if any(
			[args.tail_min_n]) > 0:
			sys.exit(
				"--tail_min_n set, so --min_n cannot be set.")

	if args.wrap_length != "None":
		if int(args.wrap_length) <= 0:
			sys.exit("--wrap_length must either be 'None' or a positive integer")

	return args


def trim_sequence(
	input_sequence, l_hard_buffer, t_hard_buffer,
	l_soft_buffer, t_soft_buffer, l_min_n, t_min_n):
	"""
	Trim Ns from a sequence.

	Inputs
	------
	input_sequence : the sequence to be trimmed
	l_hard_buffer : Number of bases to be hard clipped (no matter what) from
		beginning of sequence before looking for Ns. Note that even if Ns are
		not present, these bases will be clipped.
	t_hard_buffer : Number of bases to be hard clipped (no matter what) from
		end of sequence before looking for Ns. Note that even if Ns are
		not present, these bases will be clipped.
	l_soft_buffer : Number of bases to be dynamically clipped from
		beginning of sequence before looking for Ns. Function will move base
		by base through sequence up to this value. If a N is hit, all preceeding
		bases are removed and then trimming of Ns commences. Otherwise, if an N
		is not found, then the sequence is left as is.
	t_soft_buffer : Number of bases to be dynamically clipped from
		end of sequence before looking for Ns. Function will move base
		by base backwards through sequence up to this value. If a N is hit, all preceeding
		bases are removed and then trimming of Ns commences. Otherwise, if an N
		is not found, then the sequence is left as is.
	l_min_n : The minimum number of consecutive Ns at the beginning of a
		sequence before trimming occurs.  Count happens after the application
		of any buffers.
	t_min_n : The minimum number of consecutive Ns at the end of a
		sequence before trimming occurs.  Count happens after the application
		of any buffers.

	Output
	------
	trimmed_sequence : The processed sequence output as a string

	Notes
	-----
	buffer_set and min_n_set, despite being a bit cumbersome to enter, exist
	simply to speed up computation.
	"""
	if l_hard_buffer > 0:
		print("leading hard buffer = {}".format(l_hard_buffer))
		# Handle beginning of sequence
		trimmed_sequence = input_sequence[l_hard_buffer:]
		n_count = 0
		for base in trimmed_sequence:
			if base == "N" or base == "n":
				n_count += 1
			else:
				break
		if n_count >= l_min_n:
			trimmed_sequence = trimmed_sequence[n_count:]
		# Handle end of sequence
		if t_hard_buffer > 0:
			t_hard_buffer = t_hard_buffer * -1
			trimmed_sequence = trimmed_sequence[:t_hard_buffer]
			n_count = 0
			for base in trimmed_sequence[::-1]:
				if base == "N" or base == "n":
					n_count += 1
				else:
					break
			if n_count >= t_min_n:
				n_count = n_count * -1
				return trimmed_sequence[:n_count]
			else:
				return trimmed_sequence
		elif t_soft_buffer > 0:
			b_seen = 0
			while b_seen < t_soft_buffer:
				b_index = -1 * (b_seen + 1)
				base = trimmed_sequence[b_index]
				if base == "N" or base == "n":
					break
				else:
					b_seen += 1
			if b_seen == t_soft_buffer:
				return trimmed_sequence
			else:
				if b_seen > 0:
					trimmed_sequence = trimmed_sequence[:b_index]
				n_count = 0
				for base in trimmed_sequence[::-1]:
					if base == "N" or base == "n":
						n_count += 1
					else:
						break
				if n_count >= t_min_n:
					n_count = n_count * -1
					return trimmed_sequence[:n_count]
				else:
					return trimmed_sequence
		else:
			return trimmed_sequence
	elif l_soft_buffer > 0:
		print("leading soft buffer = {}".format(l_soft_buffer))
		# Handle beginning of a sequence
		b_seen = 0
		while b_seen < l_soft_buffer:
			base = input_sequence[b_seen]
			if base == "N" or base == "n":
				break
			else:
				b_seen += 1
		if b_seen == l_soft_buffer:
			trimmed_sequence = input_sequence
		else:
			if b_seen > 0:
				trimmed_sequence = input_sequence[b_seen:]
			else:
				trimmed_sequence = input_sequence
			n_count = 0
			for base in trimmed_sequence:
				if base == "N" or base == "n":
					n_count += 1
				else:
					break
			if n_count >= l_min_n:
				trimmed_sequence = trimmed_sequence[n_count:]
		# Handle end of sequence
		if t_hard_buffer > 0:
			t_hard_buffer = t_hard_buffer * -1
			trimmed_sequence = trimmed_sequence[:t_hard_buffer]
			n_count = 0
			for base in trimmed_sequence[::-1]:
				if base == "N" or base == "n":
					n_count += 1
				else:
					break
			if n_count >= t_min_n:
				n_count = n_count * -1
				return trimmed_sequence[:n_count]
			else:
				return trimmed_sequence
		elif t_soft_buffer > 0:
			b_seen = 0
			while b_seen < t_soft_buffer:
				b_index = -1 * (b_seen + 1)
				base = trimmed_sequence[b_index]
				if base == "N" or base == "n":
					break
				else:
					b_seen += 1
			if b_seen == t_soft_buffer:
				return trimmed_sequence
			else:
				if b_seen > 0:
					trimmed_sequence = trimmed_sequence[:b_index]
				n_count = 0
				for base in trimmed_sequence[::-1]:
					if base == "N" or base == "n":
						n_count += 1
					else:
						break
				if n_count >= t_min_n:
					n_count = n_count * -1
					return trimmed_sequence[:n_count]
				else:
					return trimmed_sequence
		else:
			return trimmed_sequence

	else:
		print("no leading buffer")
		# No leading buffer, head directly to trailing buffers
		trimmed_sequence = input_sequence
		if t_hard_buffer > 0:
			t_hard_buffer = t_hard_buffer * -1
			trimmed_sequence = trimmed_sequence[:t_hard_buffer]
			n_count = 0
			for base in trimmed_sequence[::-1]:
				if base == "N" or base == "n":
					n_count += 1
				else:
					break
			if n_count >= t_min_n:
				n_count = n_count * -1
				return trimmed_sequence[:n_count]
			else:
				return trimmed_sequence
		elif t_soft_buffer > 0:
			b_seen = 0
			while b_seen < t_soft_buffer:
				b_index = -1 * (b_seen + 1)
				base = trimmed_sequence[b_index]
				if base == "N" or base == "n":
					break
				else:
					b_seen += 1
			if b_seen == t_soft_buffer:
				return trimmed_sequence
			else:
				if b_seen > 0:
					trimmed_sequence = trimmed_sequence[:b_index]
					n_count = 0
					for base in trimmed_sequence[::-1]:
						if base == "N" or base == "n":
							n_count += 1
						else:
							break
					if n_count >= t_min_n:
						n_count = n_count * -1
						return trimmed_sequence[:n_count]
					else:
						return trimmed_sequence
		else:
			# No leading buffer AND no trailing buffer
			# Handle beginning of sequence
			n_count = 0
			for base in input_sequence:
				if base == "N" or base == "n":
					n_count += 1
				else:
					break
			if n_count >= l_min_n:
				trimmed_sequence = input_sequence[n_count:]
			#
			n_count = 0
			for base in trimmed_sequence[::-1]:
				if base == "N" or base == "n":
					n_count += 1
				else:
					break
			if n_count >= t_min_n:
				n_count = n_count * -1
				return trimmed_sequence[:n_count]
			else:
				return trimmed_sequence

if __name__ == '__main__':
	main()
