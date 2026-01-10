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

    motifs = {}
    n = len(text)
    if n < length:
        return []

    for i in range(n - length + 1):
        window = text[i: i + length]
        pd_vector = tuple(parent_distance(length, window))

        if pd_vector in motifs:
            count, positions, windows = motifs[pd_vector]
            motifs[pd_vector] = (count + 1, positions + [i], windows + [window])
        else:
            # Initialize with count 1, list of positions, list of windows
            motifs[pd_vector] = (1, [i], [window])

    # Filter by min_count and return list
    recurring = [
        (pd, (count, pos, win))
        for pd, (count, pos, win) in motifs.items()
        if count >= min_count
    ]
    return recurring

def read_text_file(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.read().strip()
            if not content: return []
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



    """
        Find recurring motifs in a sequence using Cartesian Tree Matching (CTM).

        Args:
            text (list[int]): The input sequence of integers.
            L (int): Window length to extract motifs.
            k (int): Minimum number of occurrences for a motif to be considered recurring.

        Returns:
            list[tuple]: List of recurring motifs, each as (PD_tuple, (count, positions, windows))
                         - PD_tuple: tuple representing the parent-distance array of the motif
                         - count: number of occurrences
                         - positions: starting indices of occurrences
                         - windows: actual sublists corresponding to the motif
        """
