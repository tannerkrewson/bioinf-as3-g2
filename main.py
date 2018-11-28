from readfasta import readfasta
import glob, os, sys

from alignment import pairwise_alignment, multi_seq_alignment
from upgma import generate_tree

# from output import align

def main():
    print( "*****\nBioinformatics - Assignment 3 - Group 2\n*****\n" )

    sequences = {}

    # scan in all fasta files in the "genes" directory
    os.chdir( os.getcwd() + "/genes/" )
    for file in glob.glob( "*.fasta" ):
        # read all the genes from the fasta file
        print( file )
        filename = file.split(".")[0]

        genes = readfasta( file )
        sequences[filename] = genes

    # grab the selected sequences from the cmd line args
    sequences_to_align = get_sequences_to_align_from_command_line( sequences )

    just_sequences = remove_sequence_names( sequences_to_align )

    if len(just_sequences) == 2:
        result = pairwise_alignment( just_sequences[0], just_sequences[1] )
        aligned_sequences = result[0]
        pairwise_score = result[1]
    else:
        guide_tree = generate_tree( just_sequences )
        aligned_sequences = multi_seq_alignment( just_sequences, guide_tree )

    alignments = replace_sequences_with_alignments( sequences_to_align, aligned_sequences )

    print(alignments)


def get_sequences_to_align_from_command_line( all_sequences ):
    # eg: main.py H1N1 0 1 4 H5N1 0 1
    # would align the first, second, and fifth sequences from the H1N1 file
    # and the first and second sequences from the H5N1 file all in 
    # multisequence alignment

    sequences_to_align = []
    current_file = ""

    # for each cmd line argument (0 is "main.py", so skip that)
    for i in range(1, len(sys.argv)):
        arg = sys.argv[i]
        if arg.startswith("H"):
            current_file = arg
        else:
            the_sequence = all_sequences[current_file][int(arg)]
            sequences_to_align.append(the_sequence)

    return sequences_to_align


def remove_sequence_names( sequences ):
    result = []
    for seq in sequences:
        result.append(seq[1])

    return result

def replace_sequences_with_alignments( sequences, alignments ):
    for i in range(0, len(sequences)):
        sequences[i][1] = alignments[i]

    return sequences

if __name__ == '__main__':
    main()