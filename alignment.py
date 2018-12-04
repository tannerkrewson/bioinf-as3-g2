import numpy

def align_sequences( sequence_1, sequence_2 ):
    #core code written by Kimberlyn, repurposed by Domenic
    gap_penalty = -3.5

    #initialize all of the sequences or alignments with a '-'
    sequence_1 = "-" + sequence_1
    sequence_2 = "-" + sequence_2

    scoring_matrix = numpy.zeros((len(sequence_2), len(sequence_1)))

    #initialize first row
    for j in range(0, len(sequence_1)): 
        scoring_matrix[0, j] = j * gap_penalty

    #initialize first column
    for i in range(0, len(sequence_2)): 
        scoring_matrix[i, 0] = i * gap_penalty

    #dp matrix calculation
    for i in range(1, len(sequence_2)):
        if (i % 250) == 0:
            print(int(i/len(sequence_2)*100), "%")

        for j in range(1, len(sequence_1)):
            diagonal_cell = scoring_matrix[i-1, j-1]
            diagonal_cell += calculate_sequence_cell( sequence_1, sequence_2, i, j )

            left_cell = scoring_matrix[i, j-1] + gap_penalty #value from the left
            above_cell = scoring_matrix[i-1, j] + gap_penalty #value from above

            scoring_matrix[i, j] = max(diagonal_cell, left_cell, above_cell)

    #figure out how to trace back
    i = len(sequence_2) - 1
    j = len(sequence_1) - 1

    alignment_score = scoring_matrix[i, j]
    
    traceBackDirections = ""
    while i > 0 or j > 0:
        diagonal_cell = scoring_matrix[i-1, j-1]
        diagonal_cell += calculate_sequence_cell( sequence_1, sequence_2, i, j)
        left_cell = scoring_matrix[i, j-1] + gap_penalty #value from the left
        above_cell = scoring_matrix[i-1, j] + gap_penalty #value from above

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
    seq1spot = len(sequence_1) - 1
    seq2spot = len(sequence_2) - 1
    #initialize a new list of alignments
    new_sequence_1 = ""
    new_sequence_2 = ""
    for i in range(len(traceBackDirections)):
        #choose the direction based on current instruction
        if traceBackDirections[i] == "d":
            new_sequence_1 = sequence_1[seq1spot] + \
            new_sequence_1

            new_sequence_2 = sequence_2[seq2spot] + \
            new_sequence_2

            seq1spot = seq1spot - 1
            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "u":
            new_sequence_1 = '-' + new_sequence_1

            new_sequence_2 = sequence_2[seq2spot] + \
            new_sequence_2

            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "b":
            new_sequence_1 = sequence_1[seq1spot] + \
            new_sequence_1

            new_sequence_2 = '-' + new_sequence_2

            seq1spot = seq1spot - 1

    #append the lists of alignments to each other
    aligned_sequences = [new_sequence_1, new_sequence_2, alignment_score]

    print("done!")
    return aligned_sequences

def align_alignments( alignment_1, alignment_2 ):
    #core code written by Kimberlyn, repurposed by Domenic
    gap_penalty = -3.5

    #initialize all of the sequences or alignments with a '-'
    if type(alignment_1) == list:
        for i in range(len(alignment_1)):
            alignment_1[i] = "-" + alignment_1[i]
    else:
        alignment_1 = "-" + alignment_1
    if type(alignment_2) == list:
        for i in range(len(alignment_2)):
            alignment_2[i] = "-" + alignment_2[i]
    else:
        alignment_2 = "-" + alignment_2

    scoring_matrix = numpy.zeros((len(alignment_2[0]), len(alignment_1[0])))

    #initialize first row
    for j in range(0, len(alignment_1[0])): 
        scoring_matrix[0, j] = j * (gap_penalty * len(alignment_1[0]))

    #initialize first column
    for i in range(0, len(alignment_2[0])): 
        scoring_matrix[i, 0] = i * (gap_penalty * len(alignment_2[0]))

    #dp matrix calculation
    for i in range(1, len(alignment_2[0])):
        if (i % 250) == 0:
            print(int(i/len(alignment_2[0])*100), "%")

        for j in range(1, len(alignment_1[0])):
            diagonal_cell = scoring_matrix[i-1, j-1]
            diagonal_cell += calculate_alignment_cell( alignment_1, alignment_2, i, j )

            left_cell = scoring_matrix[i, j-1] + (gap_penalty * len(alignment_1)) #value from the left
            above_cell = scoring_matrix[i-1, j] + (gap_penalty * len(alignment_2)) #value from above

            scoring_matrix[i, j] = max(diagonal_cell, left_cell, above_cell)

    #figure out how to trace back
    i = len(alignment_2[0]) - 1
    j = len(alignment_1[0]) - 1

    alignment_score = scoring_matrix[i, j]
    
    traceBackDirections = ""
    while i > 0 or j > 0:
        diagonal_cell = scoring_matrix[i-1, j-1]
        diagonal_cell += calculate_alignment_cell( alignment_1, alignment_2, i, j)
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
    seq1spot = len(alignment_1[0]) - 1
    seq2spot = len(alignment_2[0]) - 1
    #initialize a new list of alignments
    new_alignment_1 = \
    ["" for i in range(len(alignment_1))]
    new_alignment_2 = \
    ["" for i in range(len(alignment_2))]
    for i in range(len(traceBackDirections)):
        #choose the direction based on current instruction
        if traceBackDirections[i] == "d":
            for j in range(len(alignment_1)):
                new_alignment_1[j] = alignment_1[j][seq1spot] + \
                new_alignment_1[j]

            for j in range(len(alignment_2)):
                new_alignment_2[j] = alignment_2[j][seq2spot] + \
                new_alignment_2[j]

            seq1spot = seq1spot - 1
            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "u":
            for j in range(len(alignment_1)):
                new_alignment_1[j] = '-' + new_alignment_1[j]

            for j in range(len(alignment_2)):
                new_alignment_2[j] = alignment_2[j][seq2spot] + \
                new_alignment_2[j]

            seq2spot = seq2spot - 1
        if traceBackDirections[i] == "b":
            for j in range(len(alignment_1)):
                new_alignment_1[j] = alignment_1[j][seq1spot] + \
                new_alignment_1[j]

            for j in range(len(alignment_2)):
                new_alignment_2[j] = '-' + new_alignment_2[j]

            seq1spot = seq1spot - 1

    #append the lists of alignments to each other
    new_alignments = new_alignment_1 + new_alignment_2

    print("done!")
    return new_alignments

def calculate_alignment_cell( alignment_1, alignment_2, i, j ):
    #calculates the value to add to a cell given the alignments
    match_bonus = 1
    mismatch_penalty = -2
    cell_addition = 0

    for sequence_1 in alignment_1:
        for sequence_2 in alignment_2:
            if (sequence_2[i] == sequence_1[j]):
                cell_addition += match_bonus #value from diagonal if they match
            else:
                cell_addition += mismatch_penalty #value from diagonal if they don't match

    return cell_addition

def calculate_sequence_cell( sequence_1, sequence_2, i, j ):
    #calculates the value to add to a cell given the alignments
    match_bonus = 1
    mismatch_penalty = -2
    cell_addition = 0

    if (sequence_2[i] == sequence_1[j]):
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
