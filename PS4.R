#PS4
## BCH 339N Problem Set
## Sequence Alignment including BLAST. 


## Q1 - 6pt
## Sequences for the 2019-nCoV are available through NCBI: https://www.ncbi.nlm.nih.gov/genbank/2019-ncov-seqs/ Links to an external site.
## Wuhan_nCov19.fasta (see Wuhan_nCov19-1.fasta attached Download attached) contains the complete sequence of one isolate of the 2019 Coronavirus
## An open reading frame (ORF) is a stretch of sequence that starts with a start codon (ATG) 
## and ends with a stop codon (TAG, TAA, TGA) without any other stop codons in between. 
## By enumerating ORFs, we can create a list of candidate proteins encoded by this coronavirus. 
## Write a function that will find all ORFs with at least 59 aa. 
## Your function should output start and stop position of all such ORFs. 

ncov = readLines('~/Downloads/Wuhan_nCov19.fasta')
## We can collapse the multiline FASTA as follows: 
sequence = paste(ncov[2:length(ncov)], collapse= "" )

ncov_orf_finder = function (query_sequence, orf_length_threshold = 60) { 
  # vector to stored codon
  codon = c()
  reuse = c() ## if the stop codon index already use, 
  #don't add that sequence into codon(), store stop index that been used. (int)
  
  # ok , so only put the sequence into ... once you find the start codon. 
  for (i in 1: (nchar(query_sequence) - 2)) {
    codon[length(codon) + 1] = substr(query_sequence, i, i + 2)
  }
  
  orf = c()
  for (i in 1: length(codon)) {
    
    # if we find start
    if (codon[i] == "ATG") {
      
      # if we find stop , subsection it, and then find stop codon, 
      # make sure to check if it's a valid length (n >= 60)
      n = 1
      j = i + 3. # let's say j is 4 + 3 = 7, 10 

      found = FALSE
      while (found == FALSE && j <= length(codon)) {
        check = codon[j]
        if ((check == "TAG") || (check == "TAA") || (check == "TGA") ) {
          found = TRUE
          use = FALSE
          if (n >= orf_length_threshold) {
            # add to found orf
            k = 1
            while (use == FALSE && k <= length(reuse)) {
              if ( reuse[k] == (j + 2) ) {
                use = TRUE
              }
              k  = k + 1
            
            }
            ## only added the sequence if the stop codon index is not used.
            if (use == FALSE) {
              orf[length(orf) + 1] = paste(i, "-", (j + 2))
              reuse[length(reuse) + 1] = j + 2
            }
          }
        }
        j = j + 3 # lap for every 3 char is 1 codon
        n = n + 1. # length of at least 60 codons
      }
    }

  }
  return (orf) 
  #return (cat(orf, sep = "\n"))
}

ncov_orf_finder(sequence)

## Example Output: 

## Note that the order doesn't need to be identical to the example output.
# [1] "28734-28955"
# [1] "28284-28577"
# [1] "27894-28259"
# [1] "26523-27191"
# [1] "21936-22199"
# [1] "10215-10400"
# [1] "6156-6350"
# [1] "2958-3206"
# [1] "28274-29533"
# [1] "21536-25384"
# [1] "15461-15667"
# [1] "266-13483"
# [1] "27394-27759"
# [1] "27202-27387"
# [1] "26245-26472"
# [1] "25393-26220"
# [1] "13768-21555"

#81-65 + 1 = 17


## Q2 - 6pt
## A potential ORF we identified in Q1 spans the nucleotides "26523-27191"
substr(sequence, 26523, 27191)
## This is one of the proteins encoded by the new 2019 coronavirus genome. 
## Use blast to search for similar sequences.

## For database, select: Reference RNA sequences (refseq_rna). 
## You will notice that the first few hits match the Severe acute respiratory syndrome coronavirus 2.
## What other organism(s) did you find? 

## Hint: Use the "exclude" option. 
## Click on at least one such alignment and explain in your own words the output. Include a description of 
## Score, Expect and Identities and include a screenshot. 


## Q3 - 3pt
## We will use the same sequence for this question
substr(sequence, 26523, 27191)
## Use blastx to search the database. 
## Describe how blastx differs from blastn.
## Given the results page, what do you think is the most likely function of this protein? 



## Q4 - 6pt

## We talked about the intuition behind assessing statistical significance in our discussion of BLAST. 

## In this problem, we will see another example of how to think about statistical significance. 

## One of the potential ORF we identified in Q1 spans the nucleotides "266-13483"
candidate_orf = substr(sequence, 266, 13483)

## Install and load the Biostrings package. 
install.packages("Biostring")
library(BiocManager)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library (Biostrings)



## We can calculate the nucleotide frequencies with the following code
nuc_freq = letterFrequency( DNAString(candidate_orf), letters = c("A", "C", "G", "T"), as.prob = T )

nuc_freq


## Similarly, we can count trinucleotide frequencies 
tri_freq_orf = trinucleotideFrequency(DNAString(candidate_orf), as.prob =  T )
tri_freq_orf



## For a random sequence,
## we would expect its trinucleotide frequencies to be a simple multiplication of the individual base frequencies. 
## For example, Freq (ATG) = Freq(A) * Freq(T) * Freq(G) 
freq_ATG_expected = nuc_freq[1] * nuc_freq[4] * nuc_freq[3]
## One way to quantify the difference between expected and observed trinucleotide frequencies is to use the sum of squared of the differences. 
( tri_freq_orf['ATG'] - freq_ATG_expected ) ^2

## The following function quantifies this across all codons. 

diff_trinuc_freq = function (nuc_freq, tri_freq) { 
  expected_freqs = numeric(0)
    for (first in 1:4) { 
      for (second in 1:4) { 
        for (third in 1:4) { 
          expected_freqs = c(expected_freqs, nuc_freq[first] * nuc_freq[second] * nuc_freq[third] ) 
        }
      }  
    }
  return ( sum ( (tri_freq - expected_freqs)^2 ) ) 
}

candidate_orf_dif  = diff_trinuc_freq (nuc_freq,  tri_freq = tri_freq_orf)



## Please generate 1000 random sequences of length 13218  using nuc_freq.

## In other words, each sequence will be a string of length 13218 such that the probability that each nucleotide is (A,C,G,T) is equal to nuc_freq.  

#random_sequences = YOUR_CODE
random_sequences_fun = function() {
  nucleotides         <- c("A", "C", "G", "T") 
  one_thousand_sequences = c()
  for (n in 1: 1000) {
    temp13218 = ""
    for (length in 1: 13218) {
      temp13218 = paste(temp13218, sample(nucleotides, 1, prob = nuc_freq), sep = "")
    }
    one_thousand_sequences[n] = temp13218
  }
  return (one_thousand_sequences)
}

random_sequences = random_sequences_fun()

str(random_sequences)  chr [1:1000] "GCGTCCTTTCCAGGTGAGGCATTATTCCTAATATATTGGTTATCCCTTTATCACAATAGTGTCTCATGCCAATTTAAGCTTTCGTGACTCTTAGAAGTATTGACTTATCGG"| __truncated__ ...

## Note that the sequences need not have the exact single nucleotide frequencies as they are randomly generated. 

## For example

letterFrequency( DNAString(random_sequences[2]), letters = c("A", "C", "G", "T"), as.prob = T )        

A         C         G         T 

0.2972462 0. 1751400 0.2040399 0.3235739 

nuc_freq        

A         C         G         T 

0.2988349 0.1763504 0.2009381 0.3238765 


## Plot the diff_trinuc_freq of these 1000 sequences and compare it to candidate_orf_dif.
## Write a paragraph to interpret the results and how this exercise relates to statistical significant calculation in BLAST. 

random_sequences_trifreq = sapply (random_sequences, function(x){ trinucleotideFrequency (DNAString( x ), as.prob = T ) } )

random_diffs = apply (random_sequences_trifreq, 2, function(x) {diff_trinuc_freq(nuc_freq, x) } )

hist( random_diffs, xlim = c(0, 0.002))
abline ( v = candidate_orf_dif, col = "red")

## so the red vertical line (candidateorf_dif) is the difference between the probability of
# tri_freq_orf (observed in the real ORF of virus) and the multiplication of probability (expected using prob of individual base from orf) 
## and all the grays data (random_diffs) is the difference between the probability of 
# random_sequences sample of length 12318 (1000 samples based on nuc_freq) and the multiplication of probability (expected using prob. from orf)
### conclude that: for large database, if results between a real DNA, probability of (specific amino acids) is very unique, 
# vs. that of the expected probability of codons (specific amino acids). 
# The random of generating artificial sequences (even based on pro. of real DNA), came out to be bias, since 
# results are very similar to the expected probability. 
# Therefore, if BLAST search using a large query (more codons), it would increases chance that it's not random, 
# b/c of uniqueness of real DNA of large target bases. 

