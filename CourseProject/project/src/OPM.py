def compute_ranks(seq):
    sorted_seq = sorted(set(seq))
    ranks = {}
    i = 1
    for item in sorted_seq:
        ranks[item] = i
        i += 1
    return [ranks[num] for num in seq]

def order_preserving_match(text, L,k):
    n = len(text)
    # Dictionary to store motifs:
    # key = rank tuple (order-preserving representation)
    # value = tuple(count, list of start positions, list of windows)
    motifs = {}
    for i in range(n - L + 1):
        window = text[i:i+L]
        # Compute the order-preserving rank representation of the window
        prefix_ranks = tuple(compute_ranks(window))
        # If this motif has already been seen, update count, positions, and windows
        if prefix_ranks in motifs:
            count,positions,windows = motifs[prefix_ranks]
            motifs[prefix_ranks] = (count + 1, positions + [i], windows + [window])
        else:
            motifs[prefix_ranks] = (1,[i],[window])
    # Filter motifs that appear at least k times and return them as a list of tuples
    recurring_motifs = [(prefix_ranks, (count, positions, windows))
                        for prefix_ranks, (count, positions, windows) in motifs.items()
                        if count >= k]
    # Return the list of recurring motifs with their counts, positions, and windows
    return recurring_motifs

