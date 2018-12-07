def build_clade_count_dict( original_tree, clade_count_dict ):
    if type( original_tree[0] ) == tuple:
        build_clade_count_dict( original_tree[0], clade_count_dict )
    if type( original_tree[1] ) == tuple:
        build_clade_count_dict( original_tree[1], clade_count_dict )

    clade_count_dict[ original_tree ] = 0

def clade_search( boots_tree, clade_count_dict ):
    if type( boots_tree[0] ) == tuple:
        clade_search( boots_tree[0], clade_count_dict )
    if type( boots_tree[1] ) == tuple:
        clade_search( boots_tree[1], clade_count_dict )

    if boots_tree in clade_count_dict:
        clade_count_dict[ boots_tree ] += 1

def calculate_confidences( clade_count_dict, bootstrap_num ):
    #calculates the confidence of each clade based on counts and number of
    #bootstraps run
    confidences = {}

    for i, j in clade_count_dict.items():
        confidences[ i ] = j / bootstrap_num

    return confidences