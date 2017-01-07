# NextGenUtilities

A collection of scripts to perform various operations on high throughput sequencing data.  Currently organized by data type, there's otherwise no real unifying theme to these utilities (other than I wasn't able to easily find or implement this functionality using existing programs).  I'll create them as I need them, and add descriptions to the README as the scripts are complete.  Please feel free to contact me with any questions, either by filing an Issue or by sending an email.

## Fasta

### Reduce_contigs.py
This script takes a fasta file and compresses all records (generally contigs/scaffolds) smaller than a certain length in to a supercontig OR compresses all records not belonging to a list of contigs.  It outputs a final fasta containing all records larger than said length OR belonging to contigs on the include list intact, and a single final record containing all of the concatenated smaller records.  It also outputs a bed file with the coordinates of the original contigs in the new file. I created this script because I was working with a reference genome that contained about 800,000 contigs/scaffolds and couldn't be easily used by a number of programs (e.g., GATK).  I later expanded it to handle IDs, when working with a different reference with a large number of unplaced scaffolds in addition to well-defined chromosomes.

Usage:

To filter by length:
```
python Reduce_contigs.py --mode LENGTH --fasta <input_fasta> --output_fasta <output_fasta_name> --output_bed <output_bed_name> --size <threshold_size> --supercontig_name <name> --wrap_length <letters_per_output_line>
```

To filter by contig ID:
```
python Reduce_contigs.py --mode ID --fasta <input_fasta> --output_fasta <output_fasta_name> --output_bed <output_bed_name> --ids <text_file_of_ids_to_keep> --supercontig_name <name> --wrap_length <letters_per_output_line>
```

For a full list of available options, use the command
```
python Reduce_contigs.py -h
```

Requirements:
```
bioawk

python 2.7
```
Notes/Recommendations:
Right now Python 2.7 is required because of the changed .next() syntax in Python 3.  You need to manually change ```.next()``` in two lines of the "for header in faiter" loop to ```__next__()``` to make the script Python 3 compatible.

You might have some blank lines left scattered throughout (Depending on how bioawk handles your records).  You can easily remove them with the following ```sed``` command that will edit your file in place: ```sed -i '/^$/d' <output_fasta_file>```.

Also, note that this script doesn't clean up after itself, and will leave 5 temporary files in your directory.  All begin with "tmp" and will take up approximiately 2x your reference genome fasta file size in space.  These are left in case you run into any issues, as they'll allow you to check each step of the process.
