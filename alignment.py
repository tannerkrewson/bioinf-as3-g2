import numpy

from upgma import generate_tree

class Alignment:
    def __init__( self, seq_list ):
        self.name_list = []
        self.sequence_list = []

        for i in range(len(seq_list)):
            self.name_list[i] = seq_list[i][0]
            self.sequence_list[i] = seq_list[i][1]

        self.aligned_sequences = []

        self.score = 0
        self.percent_identical_sites = 0

    def align( self ):
        if len(self.sequence_list) == 2:
            self.pairwise_alignment()
        else:
            self.multi_alignment()

    def pairwise_alignment( self ):
        result = align_alignments( self.sequence_list[0], self.sequence_list[1] )

        # the first two things in the result will be the alignments
        self.aligned_sequences.append(result[0])
        self.aligned_sequences.append(result[1])

        self.score = result[2]

    def multi_alignment( self ):
        self.aligned_sequences = self.sequence_list.copy()
        
        guide_tree = generate_tree( self.sequence_list )
        self.progressive_alignment( guide_tree )

    def progressive_alignment( self, guide_tree ):
        left_tree = guide_tree[0]
        right_tree = guide_tree[1]

        # recurse if a tree member is a tuple
        # else its a leaf, so don't recurse
        if type( left_tree ) == tuple:
            alignment_pair = self.progressive_alignment( left_tree )
            seq_1 = alignment_pair[0]
        else:
            seq_1 = self.aligned_sequences[left_tree]

        if type( right_tree ) == tuple:
            alignment_pair = self.progressive_alignment( right_tree )
            seq_2 = alignment_pair[1]
        else:
            seq_2 = self.aligned_sequences[right_tree]
        
        return align_alignments( seq_1, seq_2 )

def align_alignments( alignment_1, alignment_2 ):
    #core code written by Kimberlyn, repurposed by Domenic
    gap_penalty = -3.5

    #initialize all of the sequences or alignments with a '-'
    for i in range(len(alignment_1)):
        alignment_1[i][0] = "-" + alignment_1[i][0]
    for i in range(len(alignment_2)):
        alignment_2[i][0] = "-" + alignment_2[i][0]

    scoring_matrix = numpy.zeros((len(alignment_2[0][0]), len(alignment_1[0][0])))

    #initialize first row
    for j in range(0, len(alignment_1[0][0])): 
        scoring_matrix[0, j] = j * (gap_penalty * len(alignment_1[0][0]))

    #initialize first column
    for i in range(0, len(alignment_2[0][0])): 
        scoring_matrix[i, 0] = i * (gap_penalty * len(alignment_2[0][0]))

    #dp matrix calculation
    for i in range(1, len(alignment_2[0][0])):
        if (i % 100) == 0:
            print(i)
        for j in range(1, len(alignment_1[0][0])):
            diagonal_cell = scoring_matrix[i-1, j-1]
            diagonal_cell += calculate_cell_addition( alignment_1, alignment_2, i, j )

            left_cell = scoring_matrix[i, j-1] + (gap_penalty * len(alignment_1)) #value from the left
            above_cell = scoring_matrix[i-1, j] + (gap_penalty * len(alignment_2)) #value from above

            scoring_matrix[i, j] = max(diagonal_cell, left_cell, above_cell)

    #figure out how to trace back
    i = len(alignment_2[0][0]) - 1
    j = len(alignment_1[0][0]) - 1

    alignment_score = scoring_matrix[i, j]
    
    traceBackDirections = ""
    while i > 0 or j > 0:
        diagonal_cell = scoring_matrix[i-1, j-1]
        diagonal_cell += calculate_cell_addition( alignment_1, alignment_2, i, j)
        left_cell = scoring_matrix[i, j-1] + (gap_penalty * len(alignment_1)) #value from the left
        above_cell = scoring_matrix[i-1, j] + (gap_penalty * len(alignment_2)) #value from above

        if scoring_matrix[i, j] == diagonal_cell and i >= 0 and j >= 0:
            #trace back to the diagonal
            traceBackDirections = traceBackDirections + "d" 
            i = i - 1
            j = j - 1
        elif scoring_matrix[i, j] == left_cell and i >= 0 and j >= 0:
            #trace back to the left
            traceBackDirections = traceBackDirections + "b"
            j = j - 1
        else: #scoring_matrix[i, j] == above_cell
            #trace back to above
            traceBackDirections = traceBackDirections + "u"
            i = i - 1

    #align with the trace back directions 
    seq1spot = len(alignment_1[0][0]) - 1
    seq2spot = len(alignment_2[0][0]) - 1
    #initialize a new list of alignments
    new_alignment_1 = \
    [["", alignment_1[i][1]] for i in range(len(alignment_1))]
    new_alignment_2 = \
    [["", alignment_2[i][1]] for i in range(len(alignment_2))]
    for i in range(len(traceBackDirections)):
        #choose the direction based on current instruction
        if traceBackDirections[i] == "d":
            for j in range(len(alignment_1)):
                new_alignment_1[j][0] = alignment_1[j][0][seq1spot] + \
                new_alignment_1[j][0]

            for j in range(len(alignment_2)):
                new_alignment_2[j][0] = alignment_2[j][0][seq2spot] + \
                new_alignment_2[j][0]

            seq1spot = seq1spot - 1
            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "u":
            for j in range(len(alignment_1)):
                new_alignment_1[j][0] = '-' + new_alignment_1[j][0]

            for j in range(len(alignment_2)):
                new_alignment_2[j][0] = alignment_2[j][0][seq2spot] + \
                new_alignment_2[j][0]

            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "b":
            for j in range(len(alignment_1)):
                new_alignment_1[j][0] = alignment_1[j][0][seq1spot] + \
                new_alignment_1[j][0]

            for j in range(len(alignment_2)):
                new_alignment_2[j][0] = '-' + new_alignment_2[j][0]

            seq1spot = seq1spot - 1

    #append the lists of alignments to each other
    combined_alignments = new_alignment_1 + new_alignment_2
    new_alignments = [combined_alignments, alignment_score]

    return new_alignments

def calculate_cell_addition( alignment_1, alignment_2, i, j ):
    #calculates the value to add to a cell given the alignments
    match_bonus = 1
    mismatch_penalty = -2
    cell_addition = 0

    for sequence_1 in alignment_1:
        for sequence_2 in alignment_2:
            if (sequence_2[0][i] == sequence_1[0][j]):
                cell_addition += match_bonus #value from diagonal if they match
            else:
                cell_addition += mismatch_penalty #value from diagonal if they don't match

    return cell_addition

def reorder_alignments(aligns):
    #sorts the alignments by their sequence number
    new_align_list = [[] for i in range(len(aligns))]
    
    for align in aligns:
        new_align_list[align[1]] = align

    return new_align_list
