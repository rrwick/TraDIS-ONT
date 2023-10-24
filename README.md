# Extract TraDIS ONT reads

[![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This tool finds ONT reads which match expectations for TraDIS sequencing, and it outputs those reads trimmed down to their genomic sequence.

It can be run like this:
```bash
extract_tradis_ont_reads.py --reads ont.fastq.gz --start start.fasta --end end.fasta > out.fastq
```

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




## Installation

This repo contains a stand-alone Python script, so no installation is required. You can run it directly from this repo, or copy it to somewhere in your `PATH` (e.g. `~/.local/bin`) for easier access.

There is one external requirement: [minimap2](https://github.com/lh3/minimap2). If you can run `minimap2 -h` from the command line, you should be good to go.




## Usage

```
usage: extract_tradis_ont_reads.py --reads READS --start START [--end END] [--neg NEG]
                                   [--min_id MIN_ID] [--max_gap MAX_GAP] [--neg_id NEG_ID]
                                   [--neg_gap NEG_GAP] [--threads THREADS] [-h]

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
  --threads THREADS  Threads to use for alignment (default: 8)

Help:
  -h, --help         Show this help message and exit
```



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
