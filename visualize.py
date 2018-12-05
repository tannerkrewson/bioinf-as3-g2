def alignment_output( alignment ):
  if( len(alignment.aligned_sequences) == 2 ):
    seq1 = alignment.aligned_sequences[0]
    seq2 = alignment.aligned_sequences[1]

    seq_1_name = alignment.name_list[0]
    seq_2_name = alignment.name_list[1]

    #shorten name if too long
    seq_1_name = seq_1_name[:19]
    seq_2_name = seq_2_name[:19]
    
    max_per_line = 50

    dash_count_1 = 0
    dash_count_2 = 0
    
    start_place = 0
    end_place = max_per_line
    place_mark = 0
    while (len(seq1) // end_place) >= 1:
        print(f'{seq_1_name:20} {(start_place + 1 - dash_count_1):5} {seq1[start_place : end_place]}')
        seq_2_print_index = start_place + 1 - dash_count_2

        #create the middle part
        alignment = ""
        while place_mark < (end_place):
            #print(place_mark)
            if (seq1[place_mark] == '-'):
                alignment += ' '
                dash_count_1 += 1
            elif (seq2[place_mark] == '-'):
                alignment += ' '
                dash_count_2 += 1
            elif seq1[place_mark] == seq2[place_mark]:
                alignment += '|'
            else:
                alignment += '.'
            place_mark += 1
        #print(alignment)
        print((alignment).rjust(27 + max_per_line))

    
    #print( seq_1_name + " " + seq_1[0:50] )
    #print( seq_2_name + " " + seq_2[0:50] )

        print(f'{seq_2_name:20} {(seq_2_print_index):5} {seq2[start_place : end_place]}')
        print()
        
        start_place += max_per_line
        end_place += max_per_line
        alignment = ""

    #print the leftovers if needed
    if (len(seq1) % max_per_line) != 0:  
        print(f'{seq_1_name:20} {(start_place + 1 - dash_count_1):5} {seq1[start_place:]}')

        alignment = ""
        for place_mark in range (start_place, len(seq1)):
            if (seq1[place_mark] == '-'):
                alignment += ' '
                dash_count_1 += 1
            elif (seq2[place_mark] == '-'):
                alignment += ' '
                dash_count_2 += 1
            elif seq1[place_mark] == seq2[place_mark]:
                alignment += '|'
            else:
                alignment += '.'
        print((alignment).rjust(27 + len(seq1) - start_place))
                
        print(f'{seq_2_name:20} {(start_place + 1 - dash_count_2):5} {seq2[start_place:]}')

  else:
    sequences = alignment.aligned_sequences
    sequences_name = alignment.name_list

    for i in range( len(sequences_name) ):
        sequences_name[i] = sequences_name[i][:19]

    max_line = 50
    start_printed = 0
    start_index = 0
    end = max_line
    #The index of the longest string in the list if needed
    LONG_INDEX = 0 
    #List of the amount of dashes in each sequenes initialize to 0
    dash_count = [0] * len(sequences) 

    print("Multisequence Alignment")

    #find index for longest string in list
    #(only matters if alignments are different size)
    for i in range ( len(sequences) ): 
        if ( max( len(sequences[LONG_INDEX]), len(sequences[i]) ) == len(sequences[i]) ):
          LONG_INDEX = i

    while ( len(sequences[LONG_INDEX]) != 0 ): #print the sequences
        print() #for spacing
        for i in range ( len(sequences) ):
          print(f'{sequences_name[i]:24} {(start_printed + 1 - dash_count[i]):4} {sequences[i][start_index : end]}')
          for j in range ( len(sequences[i][start_index:end]) ):
            if( sequences[i][j] == '-' ): #check to see it there's a gap
                dash_count[i] += 1
          #take off string that has been writen
          sequences[i] = sequences[i][max_line:] 
        start_printed += max_line


def identical_site_percent( sequence1, sequence2 ):
    total_sites = max(len(sequence1), len(sequence2))

    identical_sites = 0
    for i in range(total_sites):
        if sequence1[i] == sequence2[i]:
            identical_sites += 1

    return identical_sites / total_sites
