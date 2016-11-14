# NextGenUtilities

A collection of scripts to perform various operations on high throughput sequencing data.  Currently organized by data type, there's otherwise no real unifying theme to these utilities (other than I wasn't able to easily find or implement this functionality using existing programs).  I'll create them as I need them, and add descriptions to the README as the scripts are complete.  Please feel free to contact me with any questions, either by filing an Issue or by sending an email.

## Fasta

### Reduce_contigs.py
This script takes a fasta file and compresses all records (generally contigs/scaffolds) smaller than a certain length in to a supercontig.  It outputs a final fasta containing all records larger than said length intact with a single final record containing all of the concatenated smaller records.  It also outputs a bed file with the coordinates. I created this script because I was working with a reference genome that contained about 800,000 contigs/scaffolds and couldn't be easily used by a number of programs (e.g., GATK).

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

python 2.7 (python 3 should also work, but I haven't tested this yet)
```
Notes/Recommendations:
You might have some blank lines left scattered throughout (Depending on how bioawk handles your records).  You can easily remove them with the following ```sed``` command that will edit your file in place: ```sed -i '/^$/d' <output_fasta_file>```.

Also, note that this script doesn't clean up after itself, and will leave 5 temporary files in your directory.  All begin with "tmp" and will take up approximiately 2x your reference genome fasta file size in space.  These are left in case you run into any issues, as they'll allow you to check each step of the process.
