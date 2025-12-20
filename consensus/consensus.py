import string
alphabet = list(string.digits + string.ascii_letters)


def GetConsensus(motifs):
    if not motifs:
        return ""
    l = len(motifs[0])
    consensus_str = ""
    # Analyze column by column
    for i in range(l):
        # Count frequencies of characters at position i
        counts = {}
        for m in motifs:
            char = m[i]
            counts[char] = counts.get(char, 0) + 1
        # Find the character with the maximum frequency (the locally optimal choice)
        best_char = max(counts, key=counts.get)
        consensus_str += best_char
    return consensus_str

def HammingDistance(s1, s2):
    # Counts the number of positions where characters differ
    distance = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            distance += 1
    return distance
def Score(motifs):
    # Lower score is better(sum of Hamming distances to consensus)
    consensus = GetConsensus(motifs)
    total_distance = 0
    for M in motifs:
        total_distance += HammingDistance(M, consensus)
    return total_distance

def consensus(s, l):
    # S = {s₁, s₂, ..., sₜ} strings of length n
    # l = motif length
    t = len(s)
    n = len(s[0])
    best_motifs = [s1[0:l] for s1 in s]
    best_score = Score(best_motifs)
    #bestScore = ∞
    # Try each l - mer from first sequence as starting point
    for i in range (n - l + 1):
        motifs = [s[0][i: i + l]]# Start with one l-mer
    # Greedily add best l - mer from each remaining sequence
        for j in range (1,t):
            profile = BuildProfile(motifs , l)
            best_lmer = FindBestLmer(s[j], l, profile)
            motifs.append(best_lmer)
        # Evaluate this collection
        score = Score(motifs)
        if score < best_score:
            best_score = score
            best_motifs = list(motifs)
    return best_motifs


def BuildProfile(motifs,l):
    # Create | Σ | × l count matrix profile = matrix of size | Σ | × l initialized to 0
    profile = {char: [0] * l for char in alphabet}
    for m in motifs:
     for i in range (l):
         if m[i] in profile:
             profile[m[i]][i] += 1
    return profile


def FindBestLmer(s , l , profile):
    # bestScore = -∞

    best_score = -1
    best_lmer = s[0:l]
    for i in range(len(s) - l + 1):
        score = 0
        lmer = s[i: i + l]
        for j in range (l):
            score += profile[lmer[j]][j]
            if score > best_score:
                best_score = score
                best_lmer = lmer
    return best_lmer

