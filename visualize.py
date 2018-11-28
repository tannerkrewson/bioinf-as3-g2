def alignment_output( sequence1, sequence2 ):
    IDENTIFIER_LENGTH = 16
    LINE_LENGTH = 50
    SEQ_POSITION_LENGTH = 4

    longer_sequence = max(len(sequence1[1]), len(sequence2[1]))
    sequence1_pos = 0
    sequence2_pos = 0
    pad_string = " " * (IDENTIFIER_LENGTH + SEQ_POSITION_LENGTH)

    while longer_sequence > 0:
        print(f"{sequence1[0]:IDENTIFIER_LENGTH}{sequence1_pos:SEQ_POSITION_LENGTH}{sequence1[1]:LINE_LENGTH} ")
        print(pad_string, middle_match_line(sequence1[:LINE_LENGTH], sequence2[:LINE_LENGTH], sequence1_pos, sequence2_pos))
        print(f"{sequence2[0]:IDENTIFIER_LENGTH}{sequence2_pos:SEQ_POSITION_LENGTH}{sequence2[1]:LINE_LENGTH} ")
        sequence1 = sequence1[len(sequence1) - LINE_LENGTH:]
        sequence2 = sequence2[len(sequence1) - LINE_LENGTH:]

def middle_match_line( sequence1, sequence2, sequence1_pos, sequence2_pos ):
    middle_line = ""

    for i in range(len):
        if sequence1[i] == '-' or sequence2[i] == '-':
            middle_line += " "
            if sequence1[i] == '-':
                sequence1_pos += 1
            if sequence2[i] == '-':
                sequence2_pos += 1
        elif sequence1[i] == sequence2[i]:
            middle_line += "|"
        else:
            middle_line += "."

    return middle_line

def identical_site_percent( sequence1, sequence2 ):
    total_sites = max(len(sequence1), len(sequence2))

    identical_sites = 0
    for i in range(total_sites):
        if sequence1[i] == sequence2[i]:
            identical_sites += 1

    return identical_sites / total_sites