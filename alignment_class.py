from alignment import align_alignments, align_sequences
from upgma import generate_tree

class Alignment:
    def __init__( self, seq_list ):
        self.name_list = []
        self.sequence_list = []

        for i in range(len(seq_list)):
            self.name_list.append(seq_list[i][0])
            self.sequence_list.append(seq_list[i][1])

        self.aligned_sequences = []
        self.phylo_tree = ()

        self.score = 0
        self.percent_identical_sites = 0

        self.align()

    def align( self ):
        if len(self.sequence_list) == 2:
            self.pairwise_alignment()
        else:
            self.multi_alignment()

    def pairwise_alignment( self ):
        result = align_sequences( self.sequence_list[0], self.sequence_list[1] )

        # the first two things in the result will be the alignments
        self.aligned_sequences.append(result[0])
        self.aligned_sequences.append(result[1])

        self.score = result[2]

    def multi_alignment( self ):
        self.phylo_tree = generate_tree( self.sequence_list )
        self.aligned_sequences = self.progressive_alignment( self.phylo_tree )

    def progressive_alignment( self, guide_tree ):
        left_tree = guide_tree[0]
        right_tree = guide_tree[1]

        # recurse if a tree member is a tuple
        # else its a leaf, so don't recurse
        if type( left_tree ) == tuple:
            alignments = self.progressive_alignment( left_tree )
            seq_1 = alignments
        else:
            seq_1 = [self.sequence_list[left_tree]]

        if type( right_tree ) == tuple:
            alignments = self.progressive_alignment( right_tree )
            seq_2 = alignments
        else:
            seq_2 = [self.sequence_list[right_tree]]
        
        return align_alignments( seq_1, seq_2 )