import string

alphabet = list(string.digits + string.ascii_letters)


def count_frequencies(text):
    letters_bucket = {}
    for c in alphabet:  # create a bucket for all possible letters in alphabet
        letters_bucket[c] = letters_bucket.get(c, 0)
    for c in text:
        # Get current count or 0, then add 1
        letters_bucket[c] = letters_bucket.get(c, 0) + 1
    return letters_bucket


def alebian(t, p):

    current_context = f" {t}"
    matches = []
    n = len(t)
    m = len(p)
    freq_P = count_frequencies(p)
    freq_W = count_frequencies(t[0:m])

    # Count initial mismatches
    diff = 0
    for c in alphabet:
        if freq_P[c] != freq_W[c]:
            diff += 1
    if diff == 0:  # report position i
            matches.append(0)
    #Sliding - O(n)
    for i in range(1, n - m + 1):
        # Check if current window matches

        if i > n - m:  # don't get out of bounds
            break
        # Update window: remove T[i], add T[i + m]
        left_char = t[i - 1]
        right_char = t[i + m - 1]

        # Handle character leaving window
        if freq_W[left_char] == freq_P[left_char]:
            diff += 1  # was correct, becomes incorrect
        freq_W[left_char] -= 1
        if freq_W[left_char] == freq_P[left_char]:
            diff -= 1  # becomes correct

        # Handle character entering window
        if freq_W[right_char] == freq_P[right_char]:
            diff += 1  # was correct, becomes incorrect
        freq_W[right_char] += 1
        if freq_W[right_char] == freq_P[right_char]:
            diff -= 1  # becomes correct
        if diff == 0:  # report position i
            matches.append(i)
    result_data = {
        "matches": matches,
        "count": len(matches),
        "pattern_used": p
    }
    return {
        "text": current_context,
        "data": result_data
    }