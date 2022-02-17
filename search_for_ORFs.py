#!/usr/bin/python3
'''
Will find ORFs in a given sequence by first checking all frames forward and backward.
'''
# --------------------------------------------------------------------
# Get the sequence
f = open('AF516335.txt')
comment = f.readline()
print(comment)
seq = f.read()
print(seq)
n = len(seq)
print('sequence length is:', n)
print('----------------------------------------')
# ---------------------------------------------------
# a function to count matches when aligning s1 and s2 which are the same length
# ---------------------------------------------------
def matches(s1, s2, threshold=0) :
    '''
    Assumes that s1 and s2 are equal length sequences.
    Returns True if the number of mismatches when aligning them directly
    is less than or equal to threshold, False otherwise.
    '''
    n1 = len(s1)
    n2 = len(s2)
    if n1 != n2 :
        print( 'Error: the input sequences are not the same length:', s1, s2)
        exit(0)
    count = 0
    for i in range(n1) :
        if s1[i] != s2[i] :
            count += 1
    return count <= threshold
# ---------------------------------------------------
# a function to find a pattern in a given string
# Find the first instance of the pattern in the range start to end
# Threshold is the number of letters than can be wrong
# Threshold 0 means no mismatches are allowed
# ---------------------------------------------------
def findpattern(pattern, seq, start, end, increment, threshold=0) :
    '''
    Assumes that pattern is a string and seq is a string,
    start, end, increment, and threshold are integers.
    Finds the first match of the string pattern in seq
    when looking at positions in range(start, end, increment),
    where a match means mismatches <= threshold.
    Returns the index of that first match, or -1 if no match.
    '''
    n = len(pattern)
    for i in range(start, end, increment) :
        # try to match the pattern to seq[i:i + npattern]
        if matches(pattern, seq[i:i+n], threshold) :
            return i
    return -1
# ---------------------------------------------------
def findstop(seq, start, end, increment) :
    '''
    Assumes seq is a string, start, end, increment are integers.
    Finds the first exact match to a stop codon in seq
    when looking at positions in range(start, end, increment).
    Returns the index of that first match, or -1 if no match.
    '''
    stops = ('TAA', 'TAG','TGA')
    for i in range(start, end, increment) :
        triplet = seq[i:i+3]
        if triplet in stops :
            return i
    return -1
# ---------------------------------------------------
# Look for all of the ORFs in a single frame
# ---------------------------------------------------
def findorfs(seq, frame) :
  '''
  Assumes that seq is an uppercase DNA sequence
  and that frame is either 0, 1 or 2, indicating which frame to look in.
  Prints ORFs defined by starting with ATG and ending with STOP codon 
  in the specified frame.
  '''
  length = len(seq)
  i = frame
  countORF = 0
  while i < length :
    atg = findpattern('ATG', seq, i, length - 3, 3)
    if atg == -1 :
        break # since there are no more start codons in this frame
    else :
      stop = findstop(seq, atg + 3, length - 3, 3)
      if stop == -1 :
        break # since there are no more stop codons in this frame
      else :  # found an ORF
        orf = seq[atg : stop + 3]
        if len(orf) >= 100 :
            countORF += 1
            print( 'ORF' + str(countORF), 'at position', atg, 'in frame', frame, 'is:')
            print( orf, 'length', len(orf), 'bp')
            print('-------------------------------')
        i = stop + 3 # to look for another ORF in this frame
# -----------------------------------------------------------------
findorfs(seq, 0)
