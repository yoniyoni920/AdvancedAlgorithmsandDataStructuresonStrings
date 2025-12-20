from collections import defaultdict


def compute_prev(sliding_window, constants):
    last_seen = {}
    prev = []
    for i, char in enumerate(sliding_window):
        if char in constants:
            prev.append(char)
        else:
            if char in last_seen:
                prev.append(i - last_seen[char])
            else:
                prev.append(0)
            last_seen[char] = i
    return tuple(prev)  # will allow to use this as index that has alot of
    # positions in it python can use this to look at it as a single index,
    # so it can be used to say if two windows has the exact same 'indexes' its probably the same structure
    #  in conclution it allows us to instead of compareing
    #  the letters we can now compare signatures in o(1) - if in the same location


def get_mapping(example_str, target_str, constants):
    mapping = {}
    for char_ex, char_target in zip(example_str, target_str):
        if char_ex not in constants:
            mapping[char_ex] = char_target
    return mapping


def discover_parameterized_results(sequences, length, constants):
    """ defaultdict(list) : The default factory is called without arguments to produce a new
        value when a key is not present, in __getitem__ only.
         A defaultdict compares equal to a dict with the same items.
          All remaining arguments are treated the same as if they were passed to the dict constructor,
          including keyword arguments"""
    current_context = f" {sequences}"
    signatures_idx = defaultdict(list)
    motifs = {}
    for s_idx in range(len(sequences)):
        current_seq = sequences[s_idx]
        for i in range(len(current_seq) - length + 1):
            window = current_seq[i: i + length]
            sig = compute_prev(window, constants)  # if two different windows has the same sig they may be the same
            signatures_idx[sig].append((s_idx, i))  # sig is the signature look at the prev function for
            # further explanation,
            # idx is where we found the sequence ,
            # i is the locaion inside the sequence where the substring corresponds with sig.
            if sig not in motifs:
                motifs[sig] = window
    best_sig = None
    max_count = -1

    for sig in signatures_idx:
        current_count = len(signatures_idx[sig])
        if current_count > max_count:
            max_count = current_count
            best_sig = sig
    if best_sig is None:
        return None
    motifs_pattern = motifs[best_sig]
    all_matches = []
    positions_list = signatures_idx[best_sig]

    for j in range(len(positions_list)):
        s_idx, pos = positions_list[j]
        found_text = sequences[s_idx][pos: pos + length]
        mapping = get_mapping(motifs_pattern, found_text, constants)

        match_info = {
            "position": pos + 1,
            "sequence_num": s_idx + 1,
            "text": found_text,
            "mapping": mapping
        }
        all_matches.append(match_info)

    result_data = {
        "signature": best_sig,
        "count": max_count,
        "motif": motifs_pattern,
        "results": all_matches
    }


    return {
        "text": current_context,
        "data": result_data
    }