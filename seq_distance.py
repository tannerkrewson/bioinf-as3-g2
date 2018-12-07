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

