#Finding the approximate origin of the bacterial genome.
def skewArray(genome):
  skew = [0]
  score = 0
  for i in genome:
    if i == 'C':
      score -= 1
    if i == 'G':
      score +=1
    skew.append(score)
  return skew

#The minimum skew (least amount of C and most amount of G is the origin
def minimumSkew(genome):
  arr = skewArray(genome)
  return(arr.index(min(arr)))

#Returns a complementary to the reverse of our OG pattern
def reverseComplement(pattern):
  bases = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
  }
  newStrand = ""
  for base in pattern:
    newStrand += bases[base]
  return newStrand[::-1]

#Returns the # of mismatches between two patterns
def hammingDistance(window, pattern):
  x = 0
  for i in range(len(window)):
    if window[i] != pattern[i]:
      x += 1
  return x

#Returns the most frequent pattern in the Ori region
def frequentPatterns(window, maxMismatch):
  freq = []
  for i in window:
    pattern = window[i:i+9]
    for j in window:
      if window[j:j+9] == pattern:
        freq[pattern] += 1
      elif hammingDistance(window[j:j+9], pattern) <= maxMismatch:
        freq[pattern] += 1
      elif window[j:j+9] == reverseComplement(pattern):
        freq[pattern] += 1
      elif hammingDistance(window[j:j+9], reverseComplement(pattern)) <= maxMismatch:
        freq[pattern] += 1
  return freq

#Determines where else in the genome the selected pattern appears
def patternInRestOfGenome(oriStart, oriEnd, genome, pattern):
  x = []
  window = genome[0:oriStart] + genome[oriEnd:]
  for i in window:
    if (window[i:i + 9] == pattern):
      x.append(i)
  return x

#Executes the rest of the code above:)
def executeProgram(genome):

  oriStart = minimumSkew(genome)-250
  oriEnd = minimumSkew(genome)+249

  window = genome[oriStart:oriEnd]
  freq = frequentPatterns(window, 1)
  max = freq.index(max(freq))
  n = patternInRestOfGenome(oriStart, oriEnd, genome, max)
  if len(n) > 0:
    print("The DNAa box seems to be ", max, "however, this gene also appears in other parts of the genome.")
    print("This gene is present at the following position(s):")
    for i in n:
      print(i)
  else:
    print("The DNAa box is ", max)

  print("The following is the graph for the C:G ratio of your genome. The minimum value is the apporximate origin.")
  import matplotlib.pyplot as plt
  skew1 = skewArray(genome)
  plt.plot(skew1)
  plt.xlabel('genome position')
  plt.ylabel('skew')
  plt.show()

with open('e_coli.txt') as file:
  genome = file.read()

executeProgram(genome)