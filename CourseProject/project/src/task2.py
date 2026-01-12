import random
import sys
sys.setrecursionlimit(300000)

from itertools import groupby
def build_lcp_array(text,suffix_array):
    """
        Build the Longest Common Prefix (LCP) array for a given text and its suffix array.

        LCP[i] stores the length of the longest common prefix between the
        suffixes starting at positions suffix_array[i] and suffix_array[i + 1].

        Parameters:
        ----------
        text : str
            The input string (e.g., DNA sequence).
        suffix_array : list of int
            The suffix array of the text.

        Returns:
        -------
        lcp_array : list of int
            The LCP array of length len(text) - 1.

        """
    n = len(suffix_array)
    rank = [0] * n
    lcp_array = [0] * (n - 1)
    for i in range(n):
        rank[suffix_array[i]] = i
    k = 0
    for i in range(n):
        r = rank[i]
        if r == n - 1:
            k = 0
            continue
        j = suffix_array[r + 1]
        #K: longest common prefix between two suffixes
        while i + k < n and j + k < n and text[i + k] == text[j + k]:
            k += 1
        #Longest Common Prefix
        lcp_array[r] = k
        if k > 0:
            k -= 1

    return lcp_array

def to_int_array(text):
    unique = sorted(set(text))
    mapping = {c: i+1 for i, c in enumerate(unique)}
    return [mapping[c] for c in text]


def radix_sort(arr, key_index, max_key=None):
    if max_key is None:
        max_key = max(item[key_index] for item in arr)
    count = [[] for _ in range(max_key + 1)]
    for item in arr:
        count[item[key_index]].append(item)
    result = []
    for bucket in count:
        result.extend(bucket)
    return result


def DC3(text):
    arr12 = []
    if isinstance(text[0], str):  # text is characters
        text = to_int_array(text)
    #create Triplets
    #arr12[i][3] = rank of the triplet at position arr12[i][4]
    #arr12[i][4] = original index in text
    original_len = len(text)
    text = text + [0, 0, 0]
    for i in range(original_len):
        if i % 3 == 1 or i % 3 == 2:
            if i+2 <  len(text) :
                #'i' at the end for  original position
             arr12.append([text[i],text[i+1],text[i+2],0, i])
            elif i + 1 < len(text) :
                arr12.append([text[i], text[i + 1], 0,0,i])
            else :
                arr12.append([text[i], 0, 0,0 , i])

    #radix sort arr12
    arr12 = radix_sort(arr12, 2)
    arr12 = radix_sort(arr12, 1)
    arr12 = radix_sort(arr12, 0)
    #adding rank
    rank = 1
    arr12[0][3] = 1

    for i in range(len(arr12)):
        # Compare current triplet chars (0,1,2) with previous triplet chars
        if i > 0:
            current_triplet = arr12[i][0:3]
            previous_triplet = arr12[i - 1][0:3]
            if current_triplet != previous_triplet:
                rank += 1
            # if the triplets are ident give the same rank
            arr12[i][3] = rank
    #reduced string step, used for the next steps
   # arr12.sort(key=lambda x: x[4])
    arr12.sort(key=lambda x: (x[4] % 3, x[4]))
    #turn arr12 and arr0 to arr from list
    actual_arr12 = []
    for item in arr12:
        # item[3] is the Rank
        actual_arr12.append(item[3])
    #createing recursion inorder to figure out weather the ident triplets order
    sorted_arr12 =[]
    if len(set(actual_arr12)) == len(actual_arr12):
        sa12 = [x[4] for x in sorted(arr12, key=lambda z: z[3])]
    else:

        reduced = actual_arr12
        reduced_SA = DC3(reduced)
        # Use arr12 directly because it matches the 'reduced' string order
        sa12 = [arr12[i][4] for i in reduced_SA]
    # sort b0/sa0
    sa0 = []
    rank_arr = [0] * (len(text) + 3)

    for i, idx in enumerate(sa12):
        rank_arr[idx] = i + 1



    for i in range(0, original_len, 3):
        sa0.append(i)

    sa0.sort(key=lambda idx: (
         text[idx],
         rank_arr[idx + 1],idx

      ))
    #merge sa0 and sa12
    sa = []
    p0 = 0
    p12 = 0

    while p0 < len(sa0) and p12 < len(sa12):
        i = sa0[p0]  # index in text where suffix starts
        j = sa12[p12]  # index in text where suffix starts

        if compare_suffix(i, j, text, rank_arr):
            sa.append(i)
            p0 += 1
        else:
            sa.append(j)
            p12 += 1

    # Append leftovers
    while p0 < len(sa0):
        sa.append(sa0[p0])
        p0 += 1

    while p12 < len(sa12):
        sa.append(sa12[p12])
        p12 += 1
    return sa
import sys

# Critical for large DNA sequences
sys.setrecursionlimit(300000)

def compare_suffix(i, j, text, rank):
    if j % 3 == 1:
        return (text[i], rank[i + 1] if i + 1 < len(rank) else 0) < (text[j], rank[j + 1] if j + 1 < len(rank) else 0)
    else:
        return (text[i], text[i + 1] if i + 1 < len(text) else 0, rank[i + 2] if i + 2 < len(rank) else 0) < \
               (text[j], text[j + 1] if j + 1 < len(text) else 0, rank[j + 2] if j + 2 < len(rank) else 0)


def load_data(filename):
    """
    Loads data as a sequence of tokens.
    For the algorithm, this LIST acts as the 'long string'.
    """
    dataset_sequence = []
    full_string = ""
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Skip comments and empty lines
                if line.startswith('%') or not line.strip():
                    continue

                parts = line.split()
                for p in parts:
                    if p.lower() == 'nan':
                        continue
                    dataset_sequence.append(p)

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []
    full_string = "".join(dataset_sequence)
    return full_string


def find_longest_repeated_substring(data_string):
    if not data_string:
        return 0, [], 0, []

    # 1. Build Suffix Array (data_string is a list of strings)
    suffix_array = DC3(data_string)

    # 2. Build LCP Array
    # Note: We pass the ORIGINAL list of strings to build_lcp_array so equality checks work
    lcp = build_lcp_array(data_string, suffix_array)

    if not lcp:
        return 0, [], 0, []

    # 3. Find max length
    max_len = max(lcp)

    if max_len == 0:
        return 0, [], 0, []

    # 4. Extract the substring
    best_index = lcp.index(max_len)
    start_pos = suffix_array[best_index]

    # This 'substring' is actually a LIST of strings (the sequence of tokens)
    substring = data_string[start_pos: start_pos + max_len]

    # 5. Find occurrences (simple scan of LCP array)
    positions = set()
    # Add the two suffixes that generated this LCP
    positions.add(suffix_array[best_index])
    positions.add(suffix_array[best_index + 1])

    # Scan neighbors in LCP array that have the same length
    # (Simplified scan: looks for adjacent entries in SA with LCP >= max_len)
    # Scan left
    curr = best_index - 1
    while curr >= 0 and lcp[curr] >= max_len:
        positions.add(suffix_array[curr])
        positions.add(suffix_array[curr + 1])
        curr -= 1

    # Scan right
    curr = best_index + 1
    while curr < len(lcp) and lcp[curr] >= max_len:
        positions.add(suffix_array[curr])
        positions.add(suffix_array[curr + 1])
        curr += 1

    return max_len, substring, len(positions), sorted(list(positions))


def find_most_frequent_substring(data_string, suffix_array, lcp, length):
    """
    Finds the most frequent exact motif of a specific length using DC3 & LCP.
    """
    if not data_string or len(data_string) < length:
        return 0, "", 0


    max_freq = 1
    current_freq = 1
    start_index = suffix_array[0]


    # If LCP[i] >= length, it means the suffix at SA[i] and SA[i+1]
    # share a prefix of at least 'length'. This implies an occurrence.

    for i in range(len(lcp)):
        if lcp[i] >= length:
            current_freq += 1
        else:
            # Check if this was a new max.
            if current_freq > max_freq:
                max_freq = current_freq
                # continue from suffix index of the prev item
                start_index = suffix_array[i]
            current_freq = 1

    # check last seq
    if current_freq > max_freq:
        max_freq = current_freq
        start_index = suffix_array[len(lcp)]

    motif = data_string[start_index: start_index + length]
    return length, motif, max_freq



def main():
    data_string = load_data('../dataset/Berkeley Earth global temperature.txt')
    k = 5
    candidates = []

    #print(data_string)
    s_len, sub, count, pos = find_longest_repeated_substring(data_string)
    curr_freq = 0
    best_freq = 0
    best_len = 0
    best_sub = 0
    #Longest Repeated Sequence
    print(f"Longest Repeated Sequence Length: {s_len}")
    print(f"Sequence: {sub}")
    print("-" * 30)

    #ill add task topk repeated into this task
    #The most frequent repeated substring.
    suffix_array = DC3(data_string)
    lcp = build_lcp_array(data_string, suffix_array)
    for length in range(2, len(data_string)):
        curr_len, curr_sub, curr_freq = find_most_frequent_substring(data_string, suffix_array, lcp, length)
        # If we can't find any substring of this length that appears more than once,
        # we  won't find longer ones.
        if curr_freq <= 1:
            break
        if curr_freq > best_freq:
            best_freq = curr_freq
            best_len = curr_len
            best_sub = curr_sub
        elif curr_freq == best_freq and curr_len > best_len:
            best_len = curr_len
            best_sub = curr_sub
        candidates.append({
            'length': curr_len,
            'sequence': curr_sub,
            'frequency': curr_freq
        })

    print(f"Most Frequent subSequence: {best_sub}")
    print(f"Length: {best_len}")
    print(f"Frequency: {best_freq}")
    # we can now continue with the list of canidates to find top k
    #sort

    top_k_candidates = sorted(candidates, key=lambda x: (x['frequency'], x['length']), reverse=True)[:k]
    for i, can in enumerate(top_k_candidates, 1):
        print(f"{i}. Sequence: {can['sequence']} Length: {can['length']} | Frequency: {can['frequency']}")

    print("   " + "-" * 20)
    L = 6
    print(f"maximal repeats of length ≥  {L}")
    #All maximal repeats of length ≥ L (for different values of L).(length >=l && only most Frequency)
    top_candidates = sorted(candidates, key=lambda x: (x['frequency'], x['length']), reverse=True)
    # sort by priority first frequency, then by length , only add the most left one - for each frequency highest length
    maximal_repeats = []

    for freq, group in groupby(top_candidates, key=lambda x: x['frequency']):
        keepers_for_this_freq = []
        # group is a temporary list of all candidates with this specific frequency.
        # Because we sorted by Length Descending previously, the Longest are at the top.
        for cand in group:
            is_redundant = False
        # We ONLY compare against keepers of the SAME frequency.
        # We also clear this list every time the frequency changes, keeping it small.
            for existing in keepers_for_this_freq:
                if cand['sequence'] in existing['sequence']:
                    is_redundant = True
                    break

            if not is_redundant:
                keepers_for_this_freq.append(cand)
                maximal_repeats.append(cand)

    cnt = 1
    for can in maximal_repeats:
        if can['length']>=L:
            print(f"{cnt}. Sequence: {can['sequence']} Length: {can['length']} | Frequency: {can['frequency']}")
            cnt += 1

if __name__ == "__main__":
    main()

