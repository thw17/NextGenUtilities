#NextGenUtilities

A collection of scripts to perform various operations on high throughput sequencing data.  Currently organized by data type, there's otherwise no real unifying theme to these utilities (other than I wasn't able to easily find or implement this functionality using existing programs).  I'll create them as I need them, and add descriptions to the README as the scripts are complete.  Please feel free to contact me with any questions, either by filing an Issue or by sending an email.

## Fasta

### Reduce_contigs.py
This script takes a fasta file and compresses all records (generally contigs/scaffolds) smaller than a certain length in to a supercontig.  It outputs a final fasta containing all records larger than said length intact with a single final record containing all of the concatenated smaller records.  It also outputs a bed file with the coordinatesI created this script because I was working with a reference genome that contained about 800,000 contigs/scaffolds and couldn't be easily used by a number of programs (e.g., GATK).

Usage:
```
python Reduce_contigs.py --fasta <input_fasta> --output_fasta <output_fasta_name> --output_bed <output_bed_name> --size <threshold_size> --supercontig_name <name> --wrap_length <letters_per_output_line>
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
Because it leaves the records with long enough sequences intact, if you want wrap length to be uniform throughout the final fasta, use the line wrap length in the original file for ```--wrap_length```
