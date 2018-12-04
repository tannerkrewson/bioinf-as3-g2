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
        print("running pairwise alignment")
        result = align_sequences( self.sequence_list[0], self.sequence_list[1] )

        # the first two things in the result will be the alignments
        self.aligned_sequences.append(result[0])
        self.aligned_sequences.append(result[1])

        self.score = result[2]

    def multi_alignment( self ):
        print("generating tree")
        self.phylo_tree = generate_tree( self.sequence_list )

        print("running multi sequence alignment")
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
    
    # the percent of sites identical, the overall score, 
    # the parameters used, etc
    def print_summary( self ):

        for seq in self.aligned_sequences:
            print(seq[:120])

        print("\n*** SUMMARY ***")
        print("Phylogenetic tree:", self.phylo_tree if len(self.aligned_sequences) != 2 else "n/a for pairwise alignment")
        print("Percent of identical sites:", self.get_percent_sites_identical())
        print("Overall score:", self.score if len(self.aligned_sequences) == 2 else "n/a for multi seq alignment")
        print("Parameters used:", "what does this mean")
        print("*** END SUMMARY ***\n")
        
    def get_percent_sites_identical( self ):
        equal_sites = 0
        total_sites = len(self.aligned_sequences[0])

        for i in range(total_sites):
            base_of_first_seq = self.aligned_sequences[0][i]

            if base_of_first_seq == '-':
                continue

            all_equal = True
            for seq in self.aligned_sequences:
                if base_of_first_seq != seq[i]:
                    all_equal = False
                    break

            if all_equal:
                equal_sites += 1
    
        return equal_sites / total_sites

            
