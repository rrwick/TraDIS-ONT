#!/usr/bin/env python3
"""
This script takes in a PAF file (created by minimap2) and produces an insertion site plot file
for each of the reference sequences. 

Copyright 2024 Ryan Wick (rrwick@gmail.com)
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
import gzip
import pathlib
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Generate plot files from alignments', add_help=False)

    input_args = parser.add_argument_group('Required inputs')
    input_args.add_argument('--alignments', type=pathlib.Path, required=True,
                            help='Read alignments in PAF format')
    input_args.add_argument('--ref', type=pathlib.Path, required=True,
                            help='FASTA file containing the reference sequence')
    input_args.add_argument('--out_dir', type=pathlib.Path, required=True,
                            help='Output directory for plot files (will be created if necessary)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--min_id', type=float, default=95.0,
                              help='Minimum alignment identity for start/end seqs (default: 95)')
    setting_args.add_argument('--max_start_gap', type=int, default=1,
                              help='Maximum allowed unaligned bases at the start of a read '
                                   '(default: 1)')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help',
                           help='Show this help message and exit')
    
    return parser.parse_args()


def main():
    args = get_arguments()
    ta_site_counts, ref_lengths = count_ta_sites(args.ref)
    forward_counts, reverse_counts = get_per_site_counts(ref_lengths, ta_site_counts,
                                                         args.alignments, args.min_id,
                                                         args.max_start_gap)
    create_plot_files(ref_lengths, forward_counts, reverse_counts)


def count_ta_sites(filename):
    print(f'\nCounting TA sites in {filename}:', file=sys.stderr)
    ta_counts = {}
    ref_lengths = {}
    for name, seq in iterate_fasta(filename):
        ta_count = seq.count('TA')
        print(f'  {name}: {ta_count}', file=sys.stderr)
        ta_counts[name] = ta_count
        ref_lengths[name] = len(seq)
    ta_total = sum(ta_counts.values())
    print(f'  total TA sites: {ta_total}', file=sys.stderr)
    return ta_counts, ref_lengths


def get_per_site_counts(ref_lengths, ta_site_counts, alignments_filename, min_id, max_start_gap):
    print(f'\nTallying unique insertion sites:', file=sys.stderr)
    forward_counts, reverse_counts = {}, {}
    for name, length in ref_lengths.items():
        forward_counts[name] = [0] * length
        reverse_counts[name] = [0] * length

    with open(alignments_filename, 'rt') as paf_file:
        for line in paf_file:
            a = Alignment(line)
            if a.secondary or a.percent_identity < min_id or a.query_start > max_start_gap:
                continue
            try:
                if a.strand == '+':
                    forward_counts[a.ref_name][a.ref_start] += 1
                else:  # - strand
                    reverse_counts[a.ref_name][a.ref_end-1] += 1
            except KeyError:
                sys.exit(f'\nError: {a.ref_name} not in reference genome')
            except IndexError:
                sys.exit(f'\nError: incorrect sequence length for {a.ref_name}')

    unique_site_counts = {}
    for name, length in ref_lengths.items():
        unique_sites = sum(1 if forward_counts[name][i] > 0 or reverse_counts[name][i] > 0 else 0
                           for i in range(length))
        unique_site_counts[name] = unique_sites
        percent = 100.0 * unique_sites / ta_site_counts[name]
        print(f'  {name}: {unique_sites} ({percent:.2f}%)', file=sys.stderr)
    total_sites = sum(unique_site_counts.values())
    percent = 100.0 * total_sites / sum(ta_site_counts.values())
    print(f'  total unique insertion sites: {total_sites} ({percent:.2f}%)', file=sys.stderr)

    return forward_counts, reverse_counts


def create_plot_files(ref_lengths, forward_counts, reverse_counts):
    print(f'\nCreating plot files:', file=sys.stderr)
    for ref in ref_lengths.keys():
        filename = f'{ref}.plot'
        print(f'  {filename}', file=sys.stderr)
        with open(filename, 'wt') as f:
            for forward_count, reverse_count in zip(forward_counts[ref], reverse_counts[ref]):
                f.write(f'{forward_count} {reverse_count}\n')
    print(file=sys.stderr)


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in PAF format')

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.query_cov = 100.0 * (self.query_end - self.query_start) / self.query_length

        self.secondary = 'tp:A:S' in paf_line

    def __repr__(self):
        return self.query_name + ':' + str(self.query_start) + '-' + str(self.query_end) + \
               '(' + self.strand + '), ' + \
               self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
               ' (' + ('%.3f' % self.percent_identity) + '%)'


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


def iterate_fasta(filename):
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    yield name.split()[0], ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            yield name.split()[0], ''.join(sequence)


if __name__ == '__main__':
    main()
