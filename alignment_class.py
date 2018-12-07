from alignment import align_alignments, align_sequences
from upgma import generate_tree
from bootstrapping import generate_bootstrap_genes, find_boots_distance, generate_boots_tree
from tree_analysis import build_clade_count_dict, clade_search, calculate_confidences

class Alignment:
    def __init__( self, seq_list ):
        self.name_list = []
        self.sequence_list = []
        self.aligned_sequences = []
        self.phylo_tree = ()
        self.clade_confidences = {}

        for i in range(len(seq_list)):
            self.name_list.append(seq_list[i][0])
            self.sequence_list.append(seq_list[i][1])

        self.score = 0
        self.percent_identical_sites = 0

        self.align()
        self.bootstrap()

    def align( self ):
        if len(self.sequence_list) == 2:
            self.pairwise_alignment()
        else:
            self.multi_alignment()

    def bootstrap( self ):
        BOOTSTRAP_TIMES = 20
        clade_count_dict = {}
        build_clade_count_dict( self.phylo_tree, clade_count_dict )

        for i in range ( 0, BOOTSTRAP_TIMES) :
            boots_genes = generate_bootstrap_genes( self.aligned_sequences )

            this_tree = generate_boots_tree( boots_genes )

            print("Bootstrap Tree ", i)
            print(this_tree)

            clade_search( this_tree, clade_count_dict )

        # return a dict containing the clades as keys mapped to their confidence
        self.clade_confidences = calculate_confidences( clade_count_dict, BOOTSTRAP_TIMES )
        print(self.clade_confidences)

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

        seq_1_single = False
        seq_2_single = False

        # recurse if a tree member is a tuple
        # else its a leaf, so don't recurse
        if type( left_tree ) == tuple:
            alignments = self.progressive_alignment( left_tree )
            seq_1 = alignments
        else:
            seq_1 = [self.sequence_list[left_tree]]
            seq_1_single = True

        if type( right_tree ) == tuple:
            alignments = self.progressive_alignment( right_tree )
            seq_2 = alignments
        else:
            seq_2 = [self.sequence_list[right_tree]]
            seq_2_single = True

        if seq_1_single and seq_2_single:
            return align_sequences( seq_1[0], seq_2[0] )[:2]
        else:
            return align_alignments( seq_1, seq_2 )
    
    # the percent of sites identical, the overall score, 
    # the parameters used, etc
    def print_summary( self ):
        print("\n*** SUMMARY ***")
        print("Percent of identical sites:", self.get_percent_sites_identical())
        print("Overall score:", self.score if len(self.aligned_sequences) == 2 else "n/a for multi seq alignment")
        print("Parameters used:", "mm -3; mb +2, gap -5, gap2 -3, gap3 -1")
        print("")

        if len(self.aligned_sequences) != 2:
            print("Phylogenetic tree:", self.phylo_tree)
            for i in range(len(self.name_list)):
                print(i, "=", self.name_list[i])

        print("\n*** END SUMMARY ***\n")
        
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

    def get_name( self, index ):
        return self.name_list[index]