#!/usr/bin/env python3
"""
Copyright 2023 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/TraDIS-ONT

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import argparse
import collections
import gzip
import pathlib
import subprocess
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Extract TraDIS ONT reads', add_help=False)

    input_args = parser.add_argument_group('Required inputs')
    input_args.add_argument('--reads', type=pathlib.Path, required=True,
                            help='Input ONT reads')
    input_args.add_argument('--start', type=pathlib.Path, required=True,
                            help='FASTA file containing the expected read-start sequence')
    
    optional_args = parser.add_argument_group('Optional inputs')
    optional_args.add_argument('--end', type=pathlib.Path, required=False,
                               help='FASTA file containing the expected read-end sequence')
    optional_args.add_argument('--neg', type=pathlib.Path, required=False,
                               help='FASTA file containing negative sequences (matching reads '
                                    'will be discarded)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--min_id', type=float, default=95.0,
                              help='Minimum alignment percent identity for start/end seqs '
                                   '(default: 95)')
    setting_args.add_argument('--max_gap', type=int, default=5,
                              help='Maximum allowed gap in bp for start/end seqs (default: 5)')
    setting_args.add_argument('--neg_id', type=float, default=95.0,
                              help='Minimum alignment percent identity for negative seqs '
                                   '(default: 95)')
    setting_args.add_argument('--neg_gap', type=int, default=5,
                              help='Maximum allowed gap in bp for negative seqs (default: 5)')
    setting_args.add_argument('--threads', type=int, default=8,
                              help='Threads to use for alignment (default: 8)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help',
                           help='Show this help message and exit')
    
    return parser.parse_args()


def main():
    args = get_arguments()
    extract_reads(args.reads, args.start, args.end, args.neg, args.threads,
                  args.min_id, args.max_gap, args.neg_id, args.neg_gap)


def extract_reads(reads, start, end, neg, threads, min_id, max_gap, neg_id, neg_gap):
    """
    This function contains the tool's main functionality. It's factored out of the main function
    for easier testing.
    """
    start_alignments = align_reads(reads, start, threads, min_id, max_gap)
    end_alignments = {} if end is None else align_reads(reads, end, threads, min_id, max_gap)
    neg_alignments = {} if neg is None else align_reads(reads, neg, threads, neg_id, neg_gap, True)

    total, negative, missing_start, missing_end, multiple_start, multiple_end, \
        bad_strand, bad_pos, zero_length, passed = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    print(f'\nLoading {reads}:', file=sys.stderr)
    for name, header, seq, qual in iterate_fastq(reads):
        total += 1

        if len(start_alignments[name]) == 0:
            missing_start += 1
            continue
        if len(start_alignments[name]) > 1:
            multiple_start += 1
            continue
        if end is not None and len(end_alignments[name]) == 0:
            missing_end += 1
            continue
        if end is not None and len(end_alignments[name]) > 1:
            multiple_end += 1
            continue
        start_a = start_alignments[name][0]
        end_a = None if end is None else end_alignments[name][0]

        if end is not None and start_a.strand != end_a.strand:
            bad_strand += 1
            continue

        start_trim = start_a.query_end if start_a.strand else start_a.query_start
        end_trim = None if end is None else end_a.query_start if start_a.strand else end_a.query_end

        if not start_a.strand:  # - strand
            seq, qual = reverse_complement(seq), qual[::-1]
            start_trim = len(seq) - start_trim
            end_trim = None if end is None else len(seq) - end_trim

        if end is not None and end_trim <= start_trim:
            bad_pos += 1
            continue

        if name in neg_alignments:
            negative += 1
            continue

        if end is None:
            trimmed_seq = seq[start_trim:]
            trimmed_qual = qual[start_trim:]
        else:
            trimmed_seq = seq[start_trim:end_trim]
            trimmed_qual = qual[start_trim:end_trim]
        if len(trimmed_seq) == 0:
            zero_length += 1
            continue

        print(f'{header}\n{trimmed_seq}\n+\n{trimmed_qual}')
        passed += 1

    plural = '' if total == 1 else 's'
    print(f'  {total} read{plural}', file=sys.stderr)
    print(f'\nExcluding reads:', file=sys.stderr)
    print(f'  {missing_start} for lacking a start alignment', file=sys.stderr)
    print(f'  {multiple_start} for having multiple start alignments', file=sys.stderr)
    if end is not None:
        print(f'  {missing_end} for lacking an end alignment', file=sys.stderr)
        print(f'  {multiple_end} for having multiple end alignments', file=sys.stderr)
        print(f'  {bad_strand} for start/end alignments on inconsistent strands', file=sys.stderr)
        print(f'  {bad_pos} for having the end alignment at or before the start alignment',
              file=sys.stderr)
    if neg is not None:
        print(f'  {negative} for aligning to {neg}', file=sys.stderr)
    print(f'  {zero_length} for having a trimmed length of zero', file=sys.stderr)
    plural = '' if passed == 1 else 's'
    print(f'\n{passed} read{plural} outputted\n', file=sys.stderr)


def align_reads(reads, target, threads, min_id, max_gap, negative=False):
    print(f'\nAligning {reads} to {target}:', file=sys.stderr)

    minimap2_command = ['minimap2', '-c', '-x', 'map-ont', '-t', str(threads),
                        str(target), str(reads)]
    alignments = []
    with subprocess.Popen(minimap2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            alignments.append(Alignment(line.strip()))
        error_message = p.stderr.read()
    if p.returncode != 0:
        sys.exit(f'\nError: minimap2 failed with the following error:\n{error_message}')

    plural = '' if len(alignments) == 1 else 's'
    print(f'  {len(alignments)} alignment{plural}', file=sys.stderr)

    alignments = [a for a in alignments if a.start_gap <= max_gap and a.end_gap <= max_gap]
    arg_name = 'neg_gap' if negative else 'max_gap'
    print(f'  {len(alignments)} pass --{arg_name} {max_gap}', file=sys.stderr)

    alignments = [a for a in alignments if a.percent_identity >= min_id]
    arg_name = 'neg_id' if negative else 'min_id'
    print(f'  {len(alignments)} pass --{arg_name} {min_id}', file=sys.stderr)

    alignments_by_read_name = collections.defaultdict(list)
    for a in alignments:
        alignments_by_read_name[a.query_name].append(a)
    return alignments_by_read_name


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        assert len(line_parts) >= 11

        self.query_name = line_parts[0]
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4] == '+'

        self.start_gap = int(line_parts[7])
        self.end_gap = int(line_parts[6]) - int(line_parts[8])

        self.percent_identity = 100.0 * int(line_parts[9]) / int(line_parts[10])


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def iterate_fastq(filename):
    with get_open_func(filename)(filename, 'rt') as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith('@'):
                continue
            name = line[1:].split()[0]
            header = line
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            assert len(sequence) == len(qualities)
            yield name, header, sequence, qualities


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


if __name__ == '__main__':
	main()


def get_test_files():
    test_dir = pathlib.Path(__file__).parent.resolve() / 'tests'
    start = test_dir / 'start.fasta'
    end = test_dir / 'end.fasta'
    neg = test_dir / 'negative.fasta'
    return test_dir, start, end, neg


def test_reverse_complement():
    assert reverse_complement('ACGACTACG') == 'CGTAGTCGT'
    assert reverse_complement('AXX???XXG') == 'CNN???NNT'


def test_extract_reads_01a(capfd):
    """
    This read has neither the start seq nor the end seq, so it shouldn't be in the output.
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_01.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_01b(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_01.fastq'
    extract_reads(reads, start, None, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_02a(capfd):
    """
    This read has:
    * The entire start seq (`ACTTG...TGTTA`)
    * A middle bit (`GAACA...TTCAT`)
    * The entire end seq (`AGATC...AGAAA`)
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_02.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GAACA')
    assert out_lines[1].endswith('TTCAT')
    assert len(out_lines[1]) == len(out_lines[3]) == 1000
    assert '1 read outputted' in err


def test_extract_reads_02b(capfd):
    """
    This read has:
    * The entire start seq (`ACTTG...TGTTA`)
    * A middle bit (`GAACA...TTCAT`)
    * The entire end seq (`AGATC...AGAAA`)
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_02.fastq'
    extract_reads(reads, start, None, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GAACA')
    assert out_lines[1].endswith('AGAAA')
    assert len(out_lines[1]) == len(out_lines[3]) == 1060
    assert '1 read outputted' in err


def test_extract_reads_03a(capfd):
    """
    This read is identical to test_02.fastq but on the other strand. It should be in the output and
    have it strand flipped to match test_02.fastq.
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_03.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GAACA')
    assert out_lines[1].endswith('TTCAT')
    assert len(out_lines[1]) == len(out_lines[3]) == 1000
    assert '1 read outputted' in err


def test_extract_reads_03b(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_03.fastq'
    extract_reads(reads, start, None, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GAACA')
    assert out_lines[1].endswith('AGAAA')
    assert len(out_lines[1]) == len(out_lines[3]) == 1060
    assert '1 read outputted' in err


def test_extract_reads_04a(capfd):
    """
    This read has:
    * The entire start seq (`ACTTG...TGTTA`)
    * Additional seq in the negative file (`GCATGCAAGGTTATGCTGCT`)
    * A middle bit (`CCTAT...GCTAA`)
    * The entire end seq (`AGATC...AGAAA`)

    It should not be in the output if the negative file is provided. If no negative file is
    provided, it should be in the output.
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_04.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_04b(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_04.fastq'
    extract_reads(reads, start, end, None, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GCATG')
    assert out_lines[1].endswith('GCTAA')
    assert len(out_lines[1]) == len(out_lines[3]) == 1020
    assert '1 read outputted' in err


def test_extract_reads_04c(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_04.fastq'
    extract_reads(reads, start, None, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_04d(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_04.fastq'
    extract_reads(reads, start, None, None, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GCATG')
    assert out_lines[1].endswith('AGAAA')
    assert len(out_lines[1]) == len(out_lines[3]) == 1080
    assert '1 read outputted' in err


def test_extract_reads_05a(capfd):
    """
    This read is identical to test_02.fastq but is missing 4 bases at the start. It should be in
    the output unless a `--max_gap` of less than 4 is provided.
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_05.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 4, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GAACA')
    assert out_lines[1].endswith('TTCAT')
    assert len(out_lines[1]) == len(out_lines[3]) == 1000
    assert '1 read outputted' in err


def test_extract_reads_05b(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_05.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 3, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_06a(capfd):
    """
    This read is identical to test_02.fastq but has 4 mismatches in the start sequence. It should
    be in the output unless a `--min_id` of 99 or greater is provided.
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_06.fastq'
    extract_reads(reads, start, end, neg, 8, 98, 4, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GAACA')
    assert out_lines[1].endswith('TTCAT')
    assert len(out_lines[1]) == len(out_lines[3]) == 1000
    assert '1 read outputted' in err


def test_extract_reads_06b(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_06.fastq'
    extract_reads(reads, start, end, neg, 8, 99, 3, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_07a(capfd):
    """
    This read has the start seq and end seq on opposite strands, so it shouldn't be in the output
    if an end seq is used.
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_07.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_07b(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_07.fastq'
    extract_reads(reads, start, None, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('GCCGG')
    assert out_lines[1].endswith('GATCT')
    assert len(out_lines[1]) == len(out_lines[3]) == 1060
    assert '1 read outputted' in err


def test_extract_reads_08a(capfd):
    """
    This read has the end seq before the start seq, so it shouldn't be in the output if an end seq
    is used.
    """
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_08.fastq'
    extract_reads(reads, start, end, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    assert out == ''
    assert '0 reads outputted' in err


def test_extract_reads_08b(capfd):
    test_dir, start, end, neg = get_test_files()
    reads = test_dir / 'test_08.fastq'
    extract_reads(reads, start, None, neg, 8, 95, 5, 95, 5)
    out, err = capfd.readouterr()
    out_lines = out.splitlines()
    assert out_lines[1].startswith('TTGTT')
    assert out_lines[1].endswith('CGTAT')
    assert len(out_lines[1]) == len(out_lines[3]) == 296
    assert '1 read outputted' in err
