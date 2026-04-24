#!/usr/bin/env python3

'''Script to compare 2 sequences to each other based on using a reference to guide the alignment and starting points'''

import argparse
import io
import tempfile
import os
import json
import logging
import sys

from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Set the loglevel with the env variable LOGLEVEL=debug, LOGLEVEL=info, etc.
LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
logging.basicConfig(
    level=LOGLEVEL,
    stream=sys.stdout,
    format='%(asctime)s %(funcName)20s %(levelname)8s: %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S"
)

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


def count_iupac(sequence: str, iupac_bases: list) -> int:
    '''
    PURPOSE:
        Counts and returns integer number of IUPAC bases that are not N's
            Used in the process of setting up the score class
    INPUTS:
        sequence: String of the lowercase genomic nucleotide sequence
        iupac_bases: List of lowercase IUPAC bases

    RETURNS:
        Integer of the count of IUPAC bases
    '''
    count = 0
    for iupac_base in iupac_bases:
        count = count + sequence.count(iupac_base)
    return count


class score:
    '''
    Class used to keep track of the sequence and alignment metrics we are looking at reporting out

    INPUTS:
        name: (Str) Any string name
        sequence: (Str) Nucleotide sequence as a string
    '''
    iupac_codes=['r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v']

    # Genomic locations (1-index) with one added to the end to capture end location
    orf_range_list = [ 
                        range(108, 1685),
                        range(1807, 2705),
                        range(1807, 3330),
                        range(1829, 2389),
                        range(3438, 4445),
                        range(5458, 7110),
                        range(7271, 9124),
                        range(9234, 15785)
                    ]

    def __init__(self, name, sequence, solution_length):
       # Input data
       self.name = name
       self.sequence = sequence
       self.solution_length = solution_length

       # Basic sequence statistics
       self.sequence_length = len(self.sequence)
       self.start_n = 0
       self.end_n_position = 0
       self.end_n = 0
       self.internal_n = 0
       self.total_n = self.sequence.count('n')
       self.iupac_bases = count_iupac(self.sequence, self.iupac_codes)

       # Alignment statistics
       self.alignment_length = 0
       self.matches = 0
       ##### Ambiguous matches occur where an iupac base is in the input sequence and the expansion of that base is found
       # in the reference. For example a K expands to G and T so if a G is found in the reference sequence, this
       # would qualify as an ambiguous match and be awared 0.5 as a score (due to 1/2 = 0.5)
       self.ambiguous_matches = 0
       self.ambiguous_matches_score = 0
       #####
       self.mismatches = 0
       self.insertions = 0
       self.deletions = 0
       # If there are any unexpected indels they will be added during the analyze_alignment function with a subtype
       #    ex. {5000: {'type': 'insertion', 'length': 1, 'subtype': 'errant'}}
       self.indel_position_dict = None

       # Calculated statistics
       self.pairwise_agreement = 0
       self.genome_completeness = 0
       self.frameshift = False

       # Other
       self.reference_iupac = 0
       self.ncbi_ns = 0
       self.mismatches_list = []
       self.ambiguous_position_list = []
       self.reference_start_n = 0
       self.reference_end_n = 0


    # Functions
    def set_n_locations(self, aligned_seq):
        '''
        Finds the starting called nucleotide base locations from the aligned sequence on both the 5' and 3' ends of the genome
            Assumes that `-` at the start and end are unknown data just like an N would be
        '''
        # Find and count the Ns at the start or end based on the alignment
        self.start_n = find_N(aligned_seq)
        self.end_n = find_N(aligned_seq, start=False)
        end_n_alignment = self.alignment_length - self.end_n
        self.internal_n = str(aligned_seq)[self.start_n:end_n_alignment].count('n')

        # End N position is based on the input genome not the alignment so have to use the unaligned seq
        self.end_n_position = self.sequence_length - find_N(self.sequence, start=False)

        # For when there is no Ns at the start and the end, need to account for that missing data
        self.total_n = self.internal_n + self.start_n + self.end_n
    
    def set_reference_n_locations(self, aligned_ref_seq):
        '''
        Set where the reference starts and ends due to how the wgs minus termini
            Makes some of the input sequences longer than the ref
        '''
        self.reference_start_n = find_N(aligned_ref_seq)
        self.reference_end_n = find_N(aligned_ref_seq, start=False)

    def get_pairwise_agreement(self):
        pid = (self.matches + self.ambiguous_matches_score - self.deletions) / (self.matches + self.ambiguous_matches + self.mismatches + self.insertions) * 100
        self.pairwise_agreement = '{:.4f}'.format(pid)

    def get_genome_completeness(self):
        gc = (1 - ((self.total_n + self.deletions - self.reference_start_n - self.reference_end_n) / (self.solution_length - self.insertions))) * 100
        self.genome_completeness = '{:.4f}'.format(gc)

    def check_frameshift(self):
        for genomic_loc_start in self.indel_position_dict:
            # Check if errant, missing shouldn't cause issues with frameshifts as they'd match the reference which doesn't have any
            if self.indel_position_dict[genomic_loc_start]['subtype'] == 'errant':
                genomic_loc_end = genomic_loc_start + self.indel_position_dict[genomic_loc_start]['length'] -1 # Minus one due to the first position part of length
                # Check if start or end is in an ORF
                orf_contained_range = [x for x in self.orf_range_list if (genomic_loc_start in x) or (genomic_loc_end in x)]
                if orf_contained_range != []:
                    overlap = range(max(orf_contained_range[0][0], genomic_loc_start), min(orf_contained_range[0][-1], genomic_loc_end)+1)
                    if len(overlap) %3 != 0:
                        self.frameshift = True

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--external', required=True)
    parser.add_argument('--reference',  required=True)
    parser.add_argument('--name',  required=False)
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


def analyze_alignment(alignment, iupac_bases_dictionary: dict, score_obtained: score) -> None:
    '''
    PURPOSE:
        Compare the pairwise alignment position by position to the Wuhan-1 reference and the expected sequence for the sample to determine how well it matches

    INPUTS:
        alignment: AlignIO object of sequence alignment
        iupac_bases_dictionary: Dictionary of IUPAC codes (other than N)
        score_obtained: Class that stores the found values of the samples

    RETURNS:
        None, adds to the score_obtained class
    '''
    # Find starting spot
    #  Alignment[0] is the reference
    #  Alignment[1] is what we expect to be correct
    #  alignment[2] is the sample input
    start_pos = score_obtained.start_n
    solution_n_loc = find_N(alignment[1])
    if solution_n_loc > start_pos:
        start_pos = solution_n_loc

    reference_start_n_loc = find_N(alignment[0])

    # Keep track of reference position and sample position separately
    #   Done as insertions don't increment the alignment position and deletions don't increment the sample position
    #   The s_to_r_array contains the sample position genomic index (+1) while the index of the array relates to the wuhan-1 reference genomes index
    ref_position = start_pos

    # Populate the sample array with the starting Ns as that is where we start looping and we need all positions mapped
    s_to_r_array = [i for i in range(0, start_pos+1)] # Add one to match genomic indexs
    sample_position = start_pos

    # Indel dicts for tracking
    errant_indel_loc = {}

    # Just set the current_del_ref_position variable
    current_del_ref_position = -1

    # Go through alignment
    for i in range(start_pos, alignment.get_alignment_length() - score_obtained.end_n):
        ### Sample to Reference Array creation ###
        #   Note that both the ref_position and the sample position are genomic (index + 1)
        # Insertion
        if alignment[0][i] == '-':
            logging.debug(f"Insertion: global={i+1}, ref={ref_position}")
            sample_position += 1
        # Deletion - only care for submitted sample
        elif alignment[2][i] == '-':
            logging.debug(f"Deletion: global={i+1}, sample={sample_position}")
            ref_position +=1
            s_to_r_array.append(sample_position)
        else:
            ref_position += 1
            sample_position += 1
            s_to_r_array.append(sample_position)

        ### Alignment checks ###
        # Checking for gaps
        # Errant Insertion
        if (alignment[0][i] == '-' and alignment[1][i] == '-'):
            if (alignment[0][i-1] == '-' and alignment[1][i-1] == '-'):
                logging.debug(f"Errant insertion extends: global={i+1}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[ref_position - reference_start_n_loc]['length'] += 1
            else:
                logging.debug(f"Errant insertion begins: i={i}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[ref_position - reference_start_n_loc] = {'type': 'insertion', 'length': 1, 'subtype': 'errant'}
            score_obtained.insertions += 1
        # Missing Insertion
        elif (alignment[0][i] == '-' and alignment[2][i] == '-'):
            if (alignment[0][i-1] == '-' and alignment[2][i-1] == '-'):
                logging.debug(f"Missing insertion extends: global={i+1}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[ref_position - reference_start_n_loc]['length'] += 1
            else:
                logging.debug(f"Missing insertion begins: global={i+1}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[ref_position - reference_start_n_loc] = {'type': 'insertion', 'length': 1, 'subtype': 'missing'}
            
            score_obtained.deletions += 1 # Relative to our expected sequence
        # Errant Deletion
        elif (alignment[1][i] != '-' and alignment[2][i] == '-'):
            if (alignment[1][i-1] != '-' and alignment[2][i-1] == '-'):
                logging.debug(f"Errant deletion extends: global={i+1}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[current_del_ref_position]['length'] += 1
            else:
                logging.debug(f"Errant deletion begins: global={i+1}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[ref_position - reference_start_n_loc] = {'type': 'deletion', 'length': 1, 'subtype': 'errant'}
                current_del_ref_position = ref_position - reference_start_n_loc

            score_obtained.deletions += 1
        # Missing Deletion
        elif (alignment[1][i] == '-' and alignment[2][i] != '-'):
            if (alignment[1][i-1] == '-' and alignment[2][i-1] != '-'):
                logging.debug(f"Missing deletion extends: global={i+1}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[current_del_ref_position]['length'] += 1
            else:
                logging.debug(f"Missing deletion begins: global={i+1}, ref={ref_position - reference_start_n_loc}")
                errant_indel_loc[ref_position - reference_start_n_loc] = {'type': 'deletion', 'length': 1, 'subtype': 'missing'}
                current_del_ref_position = ref_position - reference_start_n_loc

            score_obtained.insertions += 1 # Relative to our expected sequence

        # Expected Insertion - pass as we expect it and have covered other insertions
        elif alignment[0][i] == '-':
            pass
        # Expected Deletion - pass as we expect it
        elif (alignment[1][i] == '-' and alignment[2][i] == '-'):
            pass

        # Solution should not have internal N's or iupac bases so we are passing on these locations
        elif alignment[1][i] in iupac_bases_dictionary.keys():
            if alignment[1][i] == alignment[2][i]:
                score_obtained.ambiguous_matches += 1
                score_obtained.ambiguous_matches_score += 1
            elif alignment[2][i] in iupac_bases_dictionary[alignment[1][i]]:
                score_obtained.ambiguous_matches += 1 # Count of this occurance
                score_obtained.ambiguous_matches_score += 1/len(iupac_bases_dictionary[alignment[1][i]])

            score_obtained.reference_iupac += 1
            score_obtained.ambiguous_position_list.append(f'NCBI:{alignment[1][i]}{ref_position - reference_start_n_loc}')

        # N in External NCBI data
        elif alignment[1][i] == 'n':
            score_obtained.ncbi_ns += 1

        # Internal Ns already accounted for in genome completeness calculation and can be skipped as
        # we cannot be sure if they are ambiguity or missing data
        elif alignment[2][i] == 'n':
            pass

        # Matches
        elif alignment[1][i] == alignment[2][i]:
            score_obtained.matches += 1

        # Mismatches
        else:
            # Scoring changes based on the expansion of the iupac base found. Ex. K = ['G', 'T'] which scores 0.5
            if alignment[2][i] in iupac_bases_dictionary.keys() and alignment[1][i] in iupac_bases_dictionary[alignment[2][i]]:
                score_obtained.ambiguous_matches += 1 # Count of this occurance
                score_obtained.ambiguous_matches_score += 1/len(iupac_bases_dictionary[alignment[2][i]]) # decimal score
                score_obtained.ambiguous_position_list.append(f'MeaSeq:{alignment[2][i]}{ref_position - reference_start_n_loc}')
            else:
                score_obtained.mismatches += 1
                score_obtained.mismatches_list.append(f'MeaSeq:{alignment[2][i]} {ref_position - reference_start_n_loc} NCBI:{alignment[1][i]}')

    # Check for frameshifts
    if errant_indel_loc != {}:
        score_obtained.indel_position_dict = errant_indel_loc
        score_obtained.check_frameshift()


def main():

    parser = init_parser()
    args = parser.parse_args()
    logging.debug("Begin")

    iupac_bases_dictionary = {  'r': ['a', 'g'],
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

    logging.debug("Parsing sequences.")
    # Generate SeqRecords from input files
    #  Lowercase sequences are used as the output of mafft is lowercase so want to keep that consistent
    sequences = [
        SeqIO.read(args.reference, format='fasta'),
        SeqIO.read(args.external, format='fasta'),
        SeqIO.read(args.sample, format='fasta')
    ]
    sequences[0].seq = sequences[0].seq.lower()
    sequences[1].seq = sequences[1].seq.lower()
    sequences[2].seq = sequences[2].seq.lower()

    # Just going to ignore really low coverage samples
    #  Measles genome is ~15,000 bases so at 33% completeness we can just axe the sample
    if sequences[2].seq.count('n') > 10000:
        print(' Too many Ns')
        exit(0)

    # Generate score of the input sequence to keep track of wanted metrics
    logging.debug("Running score.")
    name = 'samplex'
    if args.name:
        name = args.name
    score_obtained = score(name, str(sequences[2].seq), len(sequences[0].seq))

    # Create alignment
    logging.debug("Creating alignment.")
    alignment = concat_and_align_sequences(sequences)

    # Add simple alignment info to score:
    score_obtained.alignment_length = alignment.get_alignment_length()
    score_obtained.set_n_locations(alignment[2].seq) # The 3rd sequence is the our one we are scoring

    # Analyze alignment
    logging.debug("Analyzing alignment.")
    analyze_alignment(alignment, iupac_bases_dictionary, score_obtained)

    # Calculate final values based on alignment analysis
    logging.debug("Calculating final values.")
    score_obtained.set_reference_n_locations(alignment[0])
    score_obtained.get_genome_completeness()
    score_obtained.get_pairwise_agreement()

    if args.output_alignment:
        fname = 'alignment.fasta'
        if name:
            fname = f'{name}.alignment.fasta'
        logging.debug(f"Writing alignment: {fname}")
        with open(fname, 'w') as handle:
            SeqIO.write(alignment, handle, 'fasta')

    # Final parsable output
    score_obtained.sequence = ""
    score_obtained_json = json.dumps(score_obtained.__dict__)
    logging.debug("Done")
    print(score_obtained_json, end="\n")

if __name__ == "__main__":
    json.dumps(main())
