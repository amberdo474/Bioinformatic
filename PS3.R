## BCH 339N Problem Set
## Sequence Alignment including BLAST. 

## Q1 - 12 pt
## In class, we have discussed dynamic programming for global alignment of two sequences. 
## Please implement this algorithm. 
## Specifically, the function should take two strings in additon to scores for match, mismatch, gap penalty
## Note that all mismatches will have the same score for this question. 
## We will use a linear gap penalty (no separate penalty for gap open). 
## The output should be the optimal alignment and the associated score

## Expected Output Example: 
## ATT-AGC
## ATTCAGG
## Score: 18 
alignment_function("ATTAGC", "ATTCAGG", 6, -4, -8)

# if match_score = 6, mismatch = -4, gap_penalty = -8
alignment_function = function (str1, str2, match_score, mismatch_score, gap_penalty) { 
  
  score = matrix(nrow = nchar(str1) + 1, ncol = nchar(str2) + 1) # numeric
  path = matrix(nrow = nchar(str1) + 1, ncol = nchar(str2) + 1) # string ("left","up","diagonal")
  score[1,1] = 0
  path[1,1] = 0
  for (i in 2:ncol(score)) {
    score[1,i] = score[1,i-1] + gap_penalty
    path[1,i] = 1
  }
  for (i in 2:nrow(score)) {
    score[i,1] = score[i-1,1] + gap_penalty
    path[i,1] = 2
  }
  str1 = unlist(strsplit(paste("x",str1,sep = ""), ""))
  str2 = unlist(strsplit(paste("y", str2, sep = ""), ""))
 # str1 = unlist(strsplit("xATTAGC", ""))
  #str2 = unlist(strsplit("yATTCAGG", ""))
  # dynamic programming, doing this backward
  dynamic_function = function(x, y) {

    if ( is.na(score[x,y]) ) {
     point = mismatch_score
      if (str1[x] == str2[y]) {
        point = match_score
      } 
      direction = ""
      score[x,y] <<- (max( dynamic_function(x, y-1) + gap_penalty, 
                           dynamic_function(x-1,y) + gap_penalty, 
                           dynamic_function(x-1, y-1) + point) )
   
      direction = c(score[x,y-1] + gap_penalty, 
                    score[x-1,y] + gap_penalty, 
                    score[x-1, y-1] + point)
      direction = which.max(direction)
      path[x,y] <<- direction
    }
    return ( score[x,y] )
  }

  dynamic_function(length(str1), length(str2))
  
  # we don't know actual length until iterate through path, hit 0, 
  # starting point can be change 
  x = length(str1)
  y = length(str2)
  dna1 = ""  # x, left
  dna2 = "" # y, top  
  while (x != 1 && y != 1) {
    direction = path[x,y]
    if (direction == 1) {  #left is -
      dna1 = paste("-", dna1, sep = "")
      dna2 = paste( str2[y], dna2, sep = "")
      y = y - 1
    } else if (direction == 2) {  # top is -
      dna2 = paste("-", dna2, sep = "")
      dna1 = paste(str1[x], dna1, sep = "")
      x = x - 1
    } else {
      dna1 = paste(str1[x], dna1, sep = "")
      dna2 = paste(str2[y], dna2, sep = "")
      x = x-1
      y = y-1
    }
  }
  return (paste( dna1, dna2, score[length(str1),length(str2)], sep = ", ") )
}


    
alignment_function("ATTAGC", "ATTCAGG", 6, -4, -8)
alignment_function("ATTTT", "ATCGCG", 6, -4, -8)




## Q2 - 3pt

## The above algorithm implements a global alignment. Which line in your code would you need to modify if you wanted to implement local alignment?



## Q2 -6pt 
## Enumerate k-mers
## BLAST uses heuristics to speed up the task of searching a query sequence against a large database. 
## For example, given a word size (default 3), BLAST will create a table of all possible short words (k-mers) contained in the query. 
## Write a function that will create this table of all possible words of given size. 
## For example, given  a a word size 3 and query sequence (LRITSLRI): 
## we will have LRI, RIT, ITS, TSL, SLR, LRI => Hence 5 distinct words are possible. 

test_sequence = "LRITSLRI"
test_sequence2 = "LRITSLRIK"

enumerate_words = function(query_seq, word_size = 3) { 
  
}

enumerate_words (test_sequence) 
## Expected Output
## [1] 5
enumerate_words (test_sequence2) 
## Expected Output
## [1] 6

