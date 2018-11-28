def alignment_output():
    return False

def identical_site_percent( sequence1, sequence2 ):
    total_sites = max(len(sequence1), len(sequence2))

    identical_sites = 0
    for i in range(total_sites):
        if sequence1[i] == sequence2[i]:
            identical_sites += 1

    return identical_sites / total_sites