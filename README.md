# TraDIS ONT scripts

[![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Transposon-directed insertion-site sequencing (TraDIS) is a method for identifying genes involved in bacterial survival under specific conditions. This repository provides command-line scripts for processing Oxford Nanopore Technologies (ONT) amplicon reads generated from TraDIS experiments.

These scripts extract reads containing transposon-genome junctions, trim them and convert alignment data into genome plot files which can be visualised using [Diana](https://diana.wytamma.com), a web-based tool for TraDIS data exploration.

The method is described in detail in the following publication:<br>
[**Krause AL, Wirth W, Turner AM, Judd L, Li L, van Schaik W, Howden BP, Carter GP, Seemann T, Wick R, Stinear TP, Buultjens AH, Monk IR. Transposon-directed insertion-site sequencing (TraDIS) analysis of Enterococcus faecium using nanopore sequencing and a WebAssembly analysis platform. Microbiology Spectrum. 2025. doi:10.1128/spectrum.00628-25.**](https://doi.org/10.1128/spectrum.00628-25)

To run these scripts via a web portal (instead of the command line), visit the [Centre for Pathogen Genomics Bioinformatics Analysis Portal](https://portal.cpg.unimelb.edu.au).




## Installation

This repo contains stand-alone Python scripts, so no installation is required. You can run the scripts directly from this repo, or copy them to somewhere in your `PATH` (e.g. `~/.local/bin`) for easier access. Python 3.7 or later is required to run the scripts.

There is one external requirement: [minimap2](https://github.com/lh3/minimap2). If you can run `minimap2 -h` from the command line, you should be good to go.




## Quick usage

```bash
extract_tradis_ont_reads.py --reads ont.fastq.gz --start start.fasta --end end.fasta > trimmed.fastq
minimap2 -c -t 16 reference.fasta trimmed.fastq > alignments.paf
create_plot_files.py --ref reference.fasta --alignments alignments.paf --out_dir plot_files
```




## `extract_tradis_ont_reads.py`

This script finds ONT reads which match expectations for TraDIS sequencing, and it outputs those reads trimmed down to their genomic sequence.

Assuming this genomic configuration:
```
<-------------|---------------------------|------------->
     genome            transposon             genome
```

ONT reads will be made by amplifying DNA from the end of the transposon (start seq) to a ligated
adapter (end seq):
```
                                      >>>>>-------------<<<<<
                                      start              end
```

This tool will output reads which match this expectation, with the start/end seqs trimmed off and
(if necessary) flipped to the positive strand:
```
                                           -------------
                                            output read
```

In order for a read to be included in the output, the following criteria must be met:
* The read must have exactly one start seq alignment and one end seq alignment. These alignments
  need to meet the thresholds set by the `--min_id` and `--max_gap` settings.
* The start seq alignment and end seq alignment must be on the same strand.
* The start seq alignment must come before the end seq alignment on the read.
* If the `--neg` option is used, the read must not align to any sequence in that file. These
  alignments need to meet the thresholds set by the `--neg_id` and `--neg_gap` settings.


### Full usage
```
usage: extract_tradis_ont_reads.py --reads READS --start START [--end END] [--neg NEG]
                                   [--min_id MIN_ID] [--max_gap MAX_GAP] [--neg_id NEG_ID]
                                   [--neg_gap NEG_GAP] [--trim TRIM] [--threads THREADS] [-h]

Extract TraDIS ONT reads

Required inputs:
  --reads READS      Input ONT reads
  --start START      FASTA file containing the expected read-start sequence

Optional inputs:
  --end END          FASTA file containing the expected read-end sequence
  --neg NEG          FASTA file containing negative sequences (matching reads will be discarded)

Settings:
  --min_id MIN_ID    Minimum alignment percent identity for start/end seqs (default: 95)
  --max_gap MAX_GAP  Maximum allowed gap in bp for start/end seqs (default: 5)
  --neg_id NEG_ID    Minimum alignment percent identity for negative seqs (default: 95)
  --neg_gap NEG_GAP  Maximum allowed gap in bp for negative seqs (default: 5)
  --trim TRIM        Trim output reads so they do not exceed this length (default: do not trim
                     to a target length)
  --threads THREADS  Threads to use for alignment (default: 8)

Help:
  -h, --help         Show this help message and exit
```




## `create_plot_files.py`

This script create plot files compatible with [Artemis](https://www.sanger.ac.uk/tool/artemis) and [Diana](https://diana.wytamma.com). One plot file is made for each sequence in the reference genome, each line in the file corresponds to a position in that sequence, and there are two numbers on each line: the insertion count for the forward strand and the insertion count for the reverse strand.

The following alignments will be ignored:
* Secondary alignments (those with a `tp:A:S` tag).
* Low-identity alignments (as determined by the `--min_id` setting).
* Alignments too far from the start of the read (as determined by the `--max_gap` setting). It is therefore crucial that the reads have been properly trimmed (by `extract_tradis_ont_reads.py`) so they align from the start of their sequence.

When run, this script will also display some stats to stderr, such as the number of `TA` sites in the reference genome and the number of insertion sites found.


### Full usage
```
usage: create_plot_files.py --alignments ALIGNMENTS --ref REF --out_dir OUT_DIR
                            [--min_id MIN_ID] [--max_gap MAX_GAP] [--exclude_non_ta]
                            [--exclude_sites_below EXCLUDE_SITES_BELOW] [-h]

Generate plot files from alignments

Required inputs:
  --alignments ALIGNMENTS
                        Read alignments in PAF format
  --ref REF             FASTA file containing the reference sequence
  --out_dir OUT_DIR     Output directory (will be created if necessary)

Settings:
  --min_id MIN_ID       Minimum alignment identity for start/end seqs (default: 95)
  --max_gap MAX_GAP     Maximum allowed unaligned bases at the start of a read (default: 5)
  --exclude_non_ta      Exclude all insertions at non-TA sites
  --exclude_sites_below EXCLUDE_SITES_BELOW
                        Sites with fewer than this many insertions will be rounded down to 0
                        (default: 2, i.e. exclude sites with only 1 insertion)

Help:
  -h, --help            Show this help message and exit
```




## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
