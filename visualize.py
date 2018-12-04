def alignment_output( alignment ):
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

        print(f'{seq_2_name:20} {(start_place + 1 - dash_count_2):5} {seq2[start_place : end_place]}')
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
