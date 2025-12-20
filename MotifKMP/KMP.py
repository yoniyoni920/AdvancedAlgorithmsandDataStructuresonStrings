import sys
import time
import matplotlib.pyplot as plt
# 1. New function to read from text file
def read_file(filename):
    try:
        with open(filename, 'r') as f:
            lines = [line.strip().upper() for line in f.readlines() if line.strip()]

        # Parse based on your format:
        # Line 0: Num sequences, Line 1: Motif Length, Lines 2+: Sequences
        return int(lines[0]), int(lines[1]), lines[2:]
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)


def main():
    #use cd to go to relevent foldir
    #write: python KMP.py test2.txt
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file.txt>")
        return
    filename = sys.argv[1]  # Get the filename from arguments
    #filename = "test2.txt"# use this instead of cmd if you feel like it
    num_lines, motif_length, string_arr = read_file(filename)
    start_time = time.time()

    pattern_arr = []
    seen_motifs = set()  # set to track processed motifs
    k = 0
    j = 0




    for i in range(len(string_arr[0]) - motif_length + 1):
        current_motif= string_arr[0][i:i + motif_length]
        # This is a 3d arr, the motif is first d,
        # The string is the second dimensions,
        # The location is the third
        # so its location of the pattern in each string in each motif

        if current_motif not in seen_motifs:
            pattern_arr.append(find_all_pattern_in_arr(num_lines,string_arr,current_motif))#I give you the motif you give me the answer
            seen_motifs.add(current_motif)
        else:
            pattern_arr.append([])
    end_time = time.time()
    runtime = end_time - start_time

    letters = Alphabet(string_arr)
    print_string(pattern_arr,motif_length,letters,string_arr[0],num_lines,runtime)

    #print(pattern_arr)
def print_string(pattern_arr,motif_length,letters,string_arr,num_lines,runtime):
    k = 1
    motif_count = 0
    print_motif = ""
    motif_exists_flag= False
    print("Exact Motif Discovery Results")
    print("Number of sequences: " + str(num_lines))
    print("Motif length: " +  str(motif_length) )
    print("Alphabet:"+ str(letters))

    # print_motif_and_number_of_patterns_in_each_arr
    for i in range(len(pattern_arr)):
            if (pattern_arr[i]):# motif candidate exists
                # check if this motif appears in ALL sequences
                appears_in_all = all(len(pattern_arr[i][j]) > 0 for j in range(num_lines))
                if appears_in_all:
                    print_motif += "Motif " + str(k) + " " + string_arr[i:i + motif_length] + "\n"
                    k = k + 1
                    for j in range(num_lines):
                        print_motif += "s" +str(j+1)+" positions: "+ str(pattern_arr[i][j]) + "\n"
                    motif_count +=1
                    motif_exists_flag = True
    print(f"Exact motifs found: {motif_count}")
    if motif_exists_flag:
        print(print_motif)
    if not motif_exists_flag:
        print("Motif not found")
    print(f"Runtime: {runtime:.4f} seconds")







def count_motifs_found(pattern_arr):
    count = 0
    for pos in pattern_arr:
        if pos:
            count += 1
    return count


def Alphabet(arr):
    letters = set()  #hash table will help find all letters in string
    for letter in arr:
        for ch in letter:
            if ch.isalpha():      # only letters
                letters.add(ch)   # add to the set (duplicates ignored)

    return sorted(letters)
#do KMP check pattern on all strings in stringArr
def find_all_pattern_in_arr(numberOfLines,stringArr,pattern):
    pattern_arr = []
    for i in range(numberOfLines):
        temp_arr = kmp(stringArr[i],pattern)#I give String you give location on motif in it
        pattern_arr.append(temp_arr)

    return pattern_arr



def kmp(text,pattern):
    n= len(text)
    m= len(pattern)
    k=0
    pattern_index_arr = []
    #get LPS array
    lps=get_lps(pattern)
    i = 0 # index for text
    j = 0 # index for pattern
    #step 2: search
    while i < n:
        if text[i] == pattern[j]:
            i = i+1
            j = j+1
        if j == m:
            #pattern found at index(i - j)
           # print("pattern found at",(i-j))
            pattern_index_arr.append(i-j)
            k = k+1
            j = lps[j-1]   #contunue searching for more occurences
        elif (i < n) and (text[i] != pattern[j]):
            if j != 0:
                j = lps[j-1]#use lps to return to closest relevent character
            else:
                i = i + 1 #move to next character in text
    return pattern_index_arr




def get_lps(pattern):
    m = len(pattern)
    lps = [0] * len(pattern)#First character has no proper prefix
    length = 0 # length of previous longest prefix-suffix,
    # also used as index to know what location we need to check next if prefix == suffix prefix
    i = 1# start from second character
    while i < m:
        if pattern[i] == pattern[length]:
            length = length + 1
            lps[i] = length
            i = i + 1
        else:
            if length != 0:
                #mismatch try shorter prefix
                length = lps[length -1]
            else:
                #no prefix possible
                lps[i] = 0
                i = i + 1
    return lps








if __name__ == "__main__":
    main()