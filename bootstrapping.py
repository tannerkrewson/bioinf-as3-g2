import random
import numpy
import multiprocessing

from seq_distance import dK2P
from upgma import calculate_upgma

def generate_bootstrap_genes( sequences ):
    new_sequences = ["" for i in sequences]

    #randomly pick columns from the progressively aligned sequences
    #and append to the new list of sequences
    for i in range(len(sequences[0])):
        random_index = random.randint(0, len(sequences[0]) - 1)
        for j in range(len(sequences)):
            new_sequences[j] += sequences[j][random_index]

    return new_sequences

def find_boots_distance( gene1, gene2, i, j):
    
    distance = dK2P( gene1, gene2 )

    return [ distance, i, j]

def generate_boots_tree( genes ):
    # find the distance between each gene
    distance_matrix = numpy.zeros((len(genes), len(genes)), dtype=float)

    def store_distance(result):
        distance = result[0]
        i = result[1]
        j = result[2]

        distance_matrix[i, j] = distance

    pool = multiprocessing.Pool()

    for i in range( 0, len( genes ) ):
        for j in range( i+1, len( genes ) ):
            pool.apply_async(find_boots_distance, args = \
                (genes[i], genes[j], i, j), callback = store_distance)
    
    pool.close()
    pool.join()

    # use upgma to generate a tree from the distance matrix
    return calculate_upgma( distance_matrix )