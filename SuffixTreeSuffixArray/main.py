import random
import sys
sys.setrecursionlimit(300000)

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

def all_repeated_sequences_lengthK(text,suffix_array, lcp_array,length):
    """
       Find all unique repeated substrings of a specific length in the text.

       Parameters:
       ----------
       text : str
           The input string.
       suffix_array : list of int
           The suffix array of the text.
       lcp_array : list of int
           The LCP array corresponding to the suffix array.
       length : int
           The length of repeated substrings to find.

       Returns:
       -------
       repeated_sequences : set of str
           A set containing all unique repeated substrings of the given length.
    """
    all_repeated_sequences_lengthK = set()
    for i in range(len(lcp_array)):
        if lcp_array[i] >= length:
            start = suffix_array[i]
            all_repeated_sequences_lengthK.add(text[start: start + length])

    return all_repeated_sequences_lengthK
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


def main():
    """
       Main function to generate a random DNA sequence, compute its suffix and LCP arrays,
       and find all repeated sequences of lengths 20, 30, 40, 50, and 100.

       Steps:
       1. Generate a random DNA sequence of length 1,000,000.
       2. Build the suffix array using DC3.
       3. Build the LCP array using Kasai's algorithm.
       4. For each desired length, find all repeated substrings.
       5. Print the counts and a few sample repeated substrings.
       """
    random.seed(42)
    dna_bases = ['A','C','G','T']
    n = 1_000_000
    text = ''.join([random.choice(dna_bases) for i in range(n)])
    text = text + "$"
    lengths = [20,30,40,50,100]
    all_repeated_sequences = {}

    suffix_array = DC3(text)
    lcp_array = build_lcp_array(text, suffix_array)
    for length in lengths:
     all_repeated_sequences[length] = all_repeated_sequences_lengthK(text,suffix_array,lcp_array,length)
     print(f"Found {len(all_repeated_sequences[length])} repeated sequences of length {length}")
     print(all_repeated_sequences[length])

    # Print some of them for inspection
    for k in lengths:
     print(f"\nSample repeated sequences of length {k}:")
     sample = list(all_repeated_sequences[k])[:5]
     for seq in sample:
         print(seq)


if __name__ == "__main__":
    main()

