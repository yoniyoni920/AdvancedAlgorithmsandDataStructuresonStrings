from collections import deque, defaultdict

from stack import Stack



def parent_distance(size,s):
    pd = [0] * size
    st = Stack()
    for i  in range (size):
        while not st.is_empty():
            value, index = st.top()
            if value <= s[i]:
                break
            st.pop()
        if st.is_empty():
            pd[i] = 0
        else:
            #if the stack is initially empty, the while loop never runs, and the variable index is never created.
            # parent_index makes sure the program wont crush in that case
            parent_value, parent_index = st.top()
            pd[i] = i - parent_index
        st.push((s[i],i))
    return pd
def compute_failure_function(size,pattern):
    m = len(pattern)
    pd = parent_distance(m,pattern)
    length = 0
    #ff = FAILURE fucntionm
    ff = [0]*m
    for i in range(1, m):
        while length > 0:
            # PD(P[i-len..i])[len+1]
            if pd[i] <= length:
                val_suffix = pd[i]
            else:
                val_suffix = 0
            #Pseudo: PD(P[1..len+1])[len+1]
            val_prefix = pd[length]
            if val_suffix == val_prefix:
                break
            else:
                length = ff[length- 1]
        if pd[i] <= length:
            val_suffix = pd[i]
        else:
            val_suffix = 0
        val_prefix = pd[length]

        if val_suffix == val_prefix:
            length = length + 1
        else:
            length = 0
        ff[i] = length
    return ff

def cartesian_tree_match(text, pattern):
    n = len(text)
    m = len(pattern)
    pd = parent_distance(m, pattern)
    pi = compute_failure_function(m, pattern)
    length = 0
    dq = deque()
    for i in range(n):
        while dq and dq[-1][0] > text[i]:
            dq.pop()
        while length != 0:
            current_window_start = i - length
            if dq:
                parent_index = dq[-1][1]
                if parent_index < current_window_start:
                    val_suffix = 0
                else:
                    val_suffix = i - parent_index
            else:
                val_suffix = 0
            val_prefix = pd[length]
            if val_suffix == val_prefix:
                break
            else:
                length = pi[length - 1]
        while dq and dq[0][1] < (i - length):
            dq.popleft()
        length += 1
        dq.append((text[i], i))
        if length == m:
            print(f"Match at position {i - m + 1}")
            length = pi[length - 1]
            while dq and dq[0][1] <= (i - length):
                dq.popleft()


def discover_ctm_motifs(text, length, min_count=2):
    motifs_counts = defaultdict(int)
    pd_example = {}
    n = len(text)
    if n < length:
        return {}, {}
    for i in range(n - length + 1):
        window = text[i: i + length]
        pd_vector = tuple(parent_distance(length, window))
        motifs_counts[pd_vector] += 1
        if pd_vector not in pd_example:
            pd_example[pd_vector] = window
    filtered_motifs = {k: v for k, v in motifs_counts.items() if v >= min_count}

    return filtered_motifs, pd_example

def read_text_file(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.read().strip()
            if not content: return []
            import re
            numbers = [int(x) for x in re.split(r'[,\s]+', content) if x]
            return numbers
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []

def read_patterns_file(filename):
    patterns = []
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    pat = [int(x) for x in line.split()]
                    patterns.append(pat)
        return patterns
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []
def print_group(title, signatures, all_originals):
    print(f"\n{title} (Total: {len(signatures)})")
    if not signatures:
        print("  (None)")
        return

    sorted_sigs = sorted(list(signatures), key=lambda x: (len(x), x))

    for sig in sorted_sigs:
        motif_values = all_originals.get(sig, 'No Original Found')
        print(f"  Len={len(sig):<2} | Motif: {motif_values}")
def main():
    text = read_text_file('text.txt')
    patterns = read_patterns_file('patterns.txt')
    all_origin = {}
    if not text or not patterns:
        print("Error: Missing text or patterns data.")
        return

    LENGTHS = [2 ,3, 5, 7, 10, 20]

    task5_signatures = set()
    for pat in patterns:
        pd_sig = tuple(parent_distance(len(pat), pat))
        task5_signatures.add(pd_sig)
        if pd_sig not in all_origin:
            all_origin[pd_sig] = pat

    task6_signatures = set()


    for L in LENGTHS:
        motifs, origin = discover_ctm_motifs(text, L, min_count=2)

        for sig in motifs.keys():
            task6_signatures.add(sig)
            all_origin[sig] = origin[sig]

    common_motifs = task5_signatures.intersection(task6_signatures)
    unique_to_task6 = task6_signatures - task5_signatures


    print_group("1. Motifs from Task 5", task5_signatures, all_origin)
    print_group("2. Motifs Found in Task 6", task6_signatures, all_origin)
    print_group("3. Common Motifs", common_motifs, all_origin)
    print_group("4. Unique to Task 6", unique_to_task6, all_origin)


if __name__ == '__main__':
    main()


