#!/usr/bin/env python3

'''Script to compare fasta file to reference to get overall sequence characteristics'''

import argparse
import io
import itertools
import tempfile
import os
import pandas as pd
import logging
import sys

from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from glob import glob
from pathlib import Path
from typing import Tuple

# Set the loglevel with the env variable LOGLEVEL=debug, LOGLEVEL=info, etc.
LOGLEVEL = os.environ.get('LOGLEVEL', 'info').upper()
logging.basicConfig(
    level=LOGLEVEL,
    stream=sys.stdout,
    format='%(asctime)s %(funcName)4s %(levelname)4s: %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S"
)

# 0-indexed gene position
#  Using the locations of the WGS minus termini references
GENES_DICT = {
    'N': range(0,1578),
    'N450': range(1127,1578),
    'P/V/C': range(1699,3233),
    'M': range(3330,4338),
    'MF-NCR': range(4338,5350),
    'F': range(5350,7003),
    'H': range(7163,9017),
    'L': range(9126,15678),
}

# For tracking sites with Ns and how many are there for the overall dataset
N_TRACK_DICT = {i:0 for i in range(1,15679)}

# Lowercase IUPAC bases and what they mean
IUPAC_BASES_DICT = {
    'r': ['a', 'g'],
    'y': ['c', 't'],
    's': ['g', 'c'],
    'w': ['a', 't'],
    'k': ['g', 't'],
    'm': ['a', 'c'],
    'b': ['c', 'g', 't'],
    'd': ['a', 'g', 't'],
    'h': ['a', 'c', 't'],
    'v': ['a', 'c', 'g'],
}


def find_N(sequence: str, start=True) -> int:
    '''
    PURPOSE:
        Find the position of N's on either the 5' and 3' end of the sequence
            based on if start is set to true or not
    INPUTS:
        sequence: String of the lowercase genomic nucleotide sequence
        start: Boolean switch to determine if were counting 5' or 3' N's

    RETURNS:
        Integer index of first non N from either the start (start = True) or end (start = False) of the sequence
    '''

    if start:
        # Nucleotide location of first non-N is equal to number of N's (or -) at the start
        for position, char in enumerate(sequence):
            if (char != 'n') and (char != '-'):
                return position
    else:
        # Iteration from the end of the sequence, position of first N is equal to the length
        # of the sequence minus the first location that is not an N when iterating from the end
        for position, char in enumerate(sequence[::-1]):
            if (char != 'n') and (char != '-'):
                return position
    return position


# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
def intervals_extract(iterable):
    """Create list of intervals with sequential numbers"""
    iterable = sorted(set(iterable))
    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]


def count_iupac(sequence: str) -> int:
    '''
    PURPOSE:
        Counts and returns integer number of IUPAC bases that are not N's
            Used in the process of setting up the score class
    INPUTS:
        sequence: String of the lowercase genomic nucleotide sequence

    RETURNS:
        Integer of the count of IUPAC bases
    '''
    count = 0
    for iupac_base in IUPAC_BASES_DICT.keys():
        count = count + sequence.count(iupac_base)
    return count


def _get_first_last_base_index(sequence: str) -> Tuple[int, int]:
    '''
    PURPOSE:
        Return the first and last nucleotide in the sequence
    INPUTS:
        sequence: String of the lowercase genomic nucleotide sequence

    RETURNS:
        Tuple[start, end] containing the starting covered base and the ending covered base 
    '''
    for position, char in enumerate(sequence):
            if (char != 'n') and (char != '-'):
                start = position
                break
    for position, char in enumerate(sequence[::-1]):
        if (char != 'n') and (char != '-'):
            end = position
            break
    return start, end


class sample_info:
    '''
    Class used to keep track of the sequence and alignment metrics we are looking at reporting out

    INPUTS:
        name: (Str) Any string name
    '''
    def __init__(self, name):
        # Initialize Class
        self.name = name
        self.sequence = ''

        # Basic sequence statistics
        self.sequence_length = 0
        self.iupac_bases = 0
        self.total_n = 0

        # Placeholders calculated using alignment
        #  As the reference is going to be shorter than most sequences input
        #  Based on how its been trimmed
        self.start_n = 0
        self.end_n = 0
        self.internal_n = 0
        self.genome_completeness = 0
        self.insertions = 0
        self.deletions = 0
        self.n_locs = []

    # Functions
    def set_seq_from_aln(self, alignment) -> None:
        '''
        Finding if the given sample is longer or not than the reference sequence
        '''
        start_offset, end_offset = _get_first_last_base_index(alignment[0].seq)
        end = alignment.get_alignment_length() - end_offset

        # Insertions as they are based on the ref
        self.insertions = alignment[0].seq[start_offset:end].count('-')

        # Set proper sequence to the class
        self.sequence = alignment[1].seq[start_offset:end]

    def calculate_metrics(self) -> None:
        '''
        Calculate metrics based on the seq in self.sequence
        '''
        self.sequence_length = len(self.sequence)

        # Find and count the Ns at the start or end based on the alignment
        self.start_n = find_N(self.sequence)
        self.end_n = find_N(self.sequence, start=False)
        from_covered_seq = self.sequence[self.start_n:self.sequence_length - self.end_n]

        self.internal_n = from_covered_seq.count('n')

        # Count total Ns
        self.total_n = self.internal_n + self.start_n + self.end_n

        # N locations
        for i, base in enumerate(self.sequence, start=1):
            if base == 'n':
                self.n_locs.append(i)
                N_TRACK_DICT[i] += 1

        n_intervals = list(intervals_extract(self.n_locs))
        flat_n_interval = [f'{x[0]}-{x[1]}' for x in n_intervals]
        self.n_locs = flat_n_interval

        # IUPACs
        self.iupac_bases = count_iupac(from_covered_seq)

        # Dels
        self.deletions = from_covered_seq.count('-')

        # Completeness
        self.genome_completeness = '{:.2f}'.format(((self.sequence_length - self.total_n) / self.sequence_length) * 100)


def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', required=True)
    parser.add_argument('--reference',  required=True)
    parser.add_argument('--out', '-o', required=False, default='fasta_seq_stats.tsv')
    parser.add_argument('--output_alignment', default=False, action='store_true', required=False)

    return parser


def concat_and_align_sequences(sequence_list: list):
    '''
    PURPOSE:
        Create AlignIO sequence alignment object between both input fasta sequences by generation of a
        temporary concatenated file and subsequent alignment with MAFFT. Concatenated file is subsequently used
        for Pangolin before removal

    INPUTS:
        sequence_list: List containing both sequences as SeqRecords structured as:
            [reference, sample]

    RETURNS:
        AlignIO sequence alignment object
    '''

    # Need a temp file to concatenate sequences fasta files and feed them into MAFFT (needs file-type obj)
    # Could just write a real output file but this works well as long as propanol will allow the temp file
    fd, concat_temp = tempfile.mkstemp()

    try:
        # Create multi-sequence fasta file for alignment
        with open(concat_temp, 'w') as handle:
            SeqIO.write(sequence_list, handle, 'fasta')

        # MAFFT options and inputs (can be adjusted if needed)
        # Found at https://biopython.org/docs/dev/api/Bio.Align.Applications.html#Bio.Align.Applications.MafftCommandline)
        mafft_cline=MafftCommandline(input=concat_temp)
        mafft_cline.adjustdirection = True

        # MAFFT output goes to stdout 
        stdout, stderr = mafft_cline()

    # Remove temp file as MAFFT output is in stdout and we no longer need the file
    finally:
        os.close(fd)
        os.remove(concat_temp)

    return AlignIO.read(io.StringIO(stdout), 'fasta')


def main():
    """Main entry point"""
    logging.debug("Begin")
    parser = init_parser()
    args = parser.parse_args()

    logging.debug("Parsing reference")
    # Generate SeqRecords from input files
    # Lowercase sequences are used as the output of mafft is lowercase so want to keep that consistent
    reference = SeqIO.read(args.reference, format='fasta')
    reference.seq = reference.seq.lower()

    final_dict_list = []
    for f in glob(f'{args.directory}/**/*.fa*', recursive=True):
        # Get the files we expect and the data we need
        logging.debug(f"Parsing {f}")
        local_path = Path(f)
        sample = SeqIO.read(local_path, format='fasta')
        sample.seq = sample.seq.lower()

        # Generate class of the input sequence to keep track of wanted metrics
        name = sample.name
        logging.info(f"Initializing initial sample info - {name}")
        sample_data = sample_info(name)

        # Create alignment
        logging.debug(f"Creating alignment - {name}")
        alignment = concat_and_align_sequences([reference, sample])

        # Calculate final values based on alignment analysis
        logging.debug(f"Calculating final values - {name}")
        sample_data.set_seq_from_aln(alignment)
        sample_data.calculate_metrics()

        if args.output_alignment:
            fname = 'alignment.fasta'
            if name:
                fname = f'{name}.alignment.fasta'
            logging.debug(f"Writing alignment - {fname}")
            with open(fname, 'w') as handle:
                SeqIO.write(alignment, handle, 'fasta')

        # Final parsable output, removing the sequence to not make the final tsv massive
        d = sample_data.__dict__
        del(d["sequence"])
        final_dict_list.append(d)
        logging.debug("Done")

    # Output
    df = pd.DataFrame(final_dict_list)
    df.to_csv(args.out, sep='\t', index=False)
    df = pd.DataFrame.from_dict(N_TRACK_DICT, orient="index")
    df.to_csv('consensus_n_summary.tsv', sep='\t')

if __name__ == "__main__":
    main()
