import math
import numpy
import random
from alignment import align_sequences

def find_distance( gene1, gene2, i, j ):

    if len(gene1) > 2000:
        newgene1 = gene1[:int(len(gene1[1])/2)]
        newgene2 = gene1[int(len(gene1[1])/2):]
        newgene3 = gene1[:int(len(gene2[1])/2)]
        newgene4 = gene1[int(len(gene2[1])/2):]

        distance = find_distance(newgene1, newgene2, i, j)[0] + \
            find_distance(newgene3, newgene4, i, j)[0] / 2
    else:
        print("finding distance between seq", i, "and seq", j)
        aligned_genes = align_sequences(gene1, gene2)
        distance = dK2P( aligned_genes[0], aligned_genes[1] )

    return [ distance, i, j ]

#x and y are the two sequences being compared 
def dK2P( seq1, seq2): 
    count_transition = 0 #the total differ by transition
    count_transversion = 0 #the total that differ by transversion

    transitions = [('G','A'), ('A', 'G'), ('C', 'T'), ('T', 'C')]
    transversions = [('G','T'), ('T', 'G'), ('A', 'T'), ('T', 'A'), \
    ('G', 'C'), ('C', 'G'), ('A', 'C'), ('C', 'A')]

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            site = (seq1[i], seq2[i])
            for j in range(len(transitions)):
                if site == transitions[j]:
                    count_transition = count_transition + 1
            for j in range(len(transversions)):
                if site == transversions[j]:
                    count_transversion = count_transversion + 1


    S = count_transition / len(seq1) 
    #fraction of sites differing by transition
    V = count_transversion / len(seq1) 
    #fraction of sites differing by transversion


    distance = -0.5 * math.log(1 - 2*S - V) - 0.25 * math.log(1 - 2*V) 
    #K2P Formula 

    return distance

'''
def align_gene(sequence1, sequence2):
    gp = -2 #gap penalty
    mb = 1  #match bonus
    mp = -1 #mismatch penalty

    #make sure they aren't the same
    if sequence1 == sequence2:
        return [sequence1, sequence2]

    #take off the parts that are the same at the begining 
    same_begin_index = 0 
    while (sequence1[same_begin_index] == sequence2[same_begin_index]) and \
            (same_begin_index <= min(len(sequence1)-1, len(sequence2)-1)):
        same_begin_index = same_begin_index + 1
    same_seq_begining = sequence1[0:same_begin_index]
    reduced_seq1 = sequence1[same_begin_index:]
    reduced_seq2 = sequence2[same_begin_index:]


    #take off the parts that are the same at the end
    same_end_index1 = len(reduced_seq1)-1
    same_end_index2 = len(reduced_seq2)-1
    while (min(same_end_index1, same_end_index2) > 0) and \
            reduced_seq1[same_end_index1] == reduced_seq2[same_end_index2]:
        same_end_index1 = same_end_index1 - 1
        same_end_index2 = same_end_index2 - 1
    same_seq_ending = reduced_seq1[same_end_index1 + 1:]

    

    #only looking at middle part that is differnt
    seq1_middle = reduced_seq1[:same_end_index1 + 1]
    seq2_middle = reduced_seq2[:same_end_index2+1]
    seq1 = "-" + seq1_middle #seq1 is the sequence on the horizontal
    seq2 = "-" + seq2_middle #seq2 is on the vertical


    A = numpy.zeros((len(seq2), len(seq1)), dtype=int) 
    #dynamic programing matrix wihtout letters


    #initialize first row
    for j in range(0, len(seq1)): 
        A[0, j] = j * gp

    #initialize first column
    for i in range(1, len(seq2)): 
        A[i, 0] = i * gp

    for i in range(1, len(seq2)):
        if (i % 100) == 0:
            print(i)
        for j in range(1, len(seq1)):
            diagonal_option  = 0
            if (seq2[i] == seq1[j]):
                diagonal_option = A[i-1, j-1] + mb 
                #value from diagonal if they match
            else:
                diagonal_option = A[i-1, j-1] + mp 
                #value from diagonal if they don't match
            left_option = A[i, j-1] + gp #value from the left
            upper_option = A[i-1, j] + gp #value from above
            A[i, j] = max(diagonal_option, left_option, upper_option)

    #figure out how to trace back
    i = len(seq2)-1
    j = len(seq1)-1

    trace_back_directions = ""
    while i > 0 or j > 0:
        if (seq2[i] == seq1[j]):
            diagonal_option = A[i-1, j-1] + mb 
            #value from diagonal if they match
        else:
            diagonal_option = A[i-1, j-1] + mp 
            #value from diagonal if they don't match
        left_option = A[i, j-1] + gp #value from the left
        upper_option = A[i-1, j] + gp #value from above

        if A[i, j] == diagonal_option and i >= 0 and j >= 0:
            #trace back to the diagonal
            trace_back_directions = trace_back_directions + "d" 
            i = i - 1
            j = j - 1
        elif A[i, j] == left_option and i >= 0 and j >= 0:
            #trace back to the left
            trace_back_directions = trace_back_directions + "b"
            j = j - 1
        else:
            #trace back to above
            trace_back_directions = trace_back_directions + "u"
            i = i - 1

    #align with the trace back directions 
    seq1_spot = len(seq1_middle) - 1
    seq2_spot = len(seq2_middle) - 1
    new_seq1 = ""
    new_seq2 = ""
    for i in range(len(trace_back_directions)):
        if trace_back_directions[i] == "d":
            new_seq1 = seq1_middle[seq1_spot] + new_seq1
            new_seq2 = seq2_middle[seq2_spot] + new_seq2
            seq1_spot = seq1_spot - 1
            seq2_spot = seq2_spot - 1
        if trace_back_directions[i] == "u":
            new_seq1 = '-' + new_seq1
            new_seq2 = seq2_middle[seq2_spot] + new_seq2
            seq2_spot = seq2_spot - 1
        if trace_back_directions[i] == "b":
            new_seq2 = '-' + new_seq2
            new_seq1 = seq1_middle[seq1_spot] + new_seq1
            seq1_spot = seq1_spot - 1

    return [same_seq_begining + new_seq1 + same_seq_ending, \
    same_seq_begining + new_seq2 + same_seq_ending]
'''

