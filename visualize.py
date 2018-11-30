def alignment_output( alignment ):
    seq_1 = alignment.alignment_list[0]
    seq_2 = alignment.alignment_list[1]

    seq_1_name = alignment.name_list[0]
    seq_2_name = alignment.name_list[1]

    IDENTIFIER_LENGTH = 16
    LINE_LENGTH = 50
    SEQ_POSITION_LENGTH = 4

    longer_sequence = max(len(seq_1), len(seq_2))
    seq_1_pos = 0
    seq_2_pos = 0
    pad_string = " " * (IDENTIFIER_LENGTH + SEQ_POSITION_LENGTH)

    while longer_sequence > 0:
        print(f"{seq_1_name:IDENTIFIER_LENGTH}{seq_1_pos:SEQ_POSITION_LENGTH}{seq_1:LINE_LENGTH} ")
        print(pad_string, middle_match_line(seq_1[:LINE_LENGTH], seq_2[:LINE_LENGTH], seq_1_pos, seq_2_pos))
        print(f"{seq_2_name:IDENTIFIER_LENGTH}{seq_2_pos:SEQ_POSITION_LENGTH}{seq_2:LINE_LENGTH} ")
        seq_1 = seq_1[len(seq_1) - LINE_LENGTH:]
        seq_2 = seq_2[len(seq_1) - LINE_LENGTH:]

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