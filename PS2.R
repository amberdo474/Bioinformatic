## BCH 339N Problem Set 2

## Please add your code in between the comments
## For this assignment please only use base R and do not rely on any packages. 


## Q1 - 6pt
## GC-content of a piece of DNA is defined as the percentage of Gs + Cs. 
## For example "ACCTGCA" has a 57.1% GC-content
## Write a function that will take a text file in FASTA format as input. 
## FASTA format: A sequence record in a FASTA format consists of a single-line description (sequence name), followed by line(s) of sequence data. 
## The first character of the description line is a greater-than (">") symbol.

example_dna

# Example input/output:
# The file "example_dna.fa Download example_dna.fa" contains:
# >T7_promoter_1
# TAATACGACTCACTATAGGG
# >T7_promoter_2
# TAATACGACTCACTATAGGGG


gc_calculator("example_dna.fa")
# [1] T7_promoter_1
# [1] 0.4

gc_calculator('~/Downloads/example_dna.fa')

## The output should be the sequence id with the lowest GC content followed by the calculated value 
gc_calculator = function ( fasta_file) { 
  fasta_file = readLines(fasta_file)
  # If code deals with multi-line fasta, that is great! However, it is not required.
  seq_id = c() # String of name, vector
  gc_content = c() # String of DNA  for now, vector
  
  for (i in 1:length(fasta_file)) {
    line = fasta_file[i]
    if ( substring(line, 1, 1) == ">"  ) { 
      seq_id[length(seq_id) + 1] = line 
    } else {
      gc_content[length(gc_content) + 1] = line
    }
  }
  count = sapply( gc_content, count_gc )
  
  print ( paste ( seq_id[which.min(count)], min(count), sep = ", ") )
}

count_gc = function(promoter1) {
  promoter1 = strsplit(promoter1, "")[[1]]
  CG_promoter1 = ( sum(promoter1 == "C") + sum(promoter1 == "G")) / length(promoter1)
  return ( CG_promoter1 )
}

## Q2 - 6pt
## Every amino acid in a protein is encoded by three nucleotides. 
## Execute the following two lines to get a list of all codons and corresponding amino acids
codons = c('UUU','UUC','UUA','UUG','UCU','UCC','UCA','UCG','UAU','UAC','UAA','UAG','UGU','UGC','UGA','UGG','CUU','CUC','CUA','CUG','CCU','CCC','CCA','CCG','CAU','CAC','CAA','CAG','CGU','CGC','CGA','CGG','AUU','AUC','AUA','AUG','ACU','ACC','ACA','ACG','AAU','AAC','AAA','AAG','AGU','AGC','AGA','AGG','GUU','GUC','GUA','GUG','GCU','GCC','GCA','GCG','GAU','GAC','GAA','GAG','GGU','GGC','GGA','GGG')
amino_acids = c('F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G' )

## Write a function that will take a coding region sequence as input. You can assume the sequence is starting with AUG and is a multiple of three nucleotides. ## mean whatt?
## The output should be the corresponding protein sequence. Report up to the first stop codon. (stop codon are: UAA, UAG, UGA)

# Example input/output:
translate_rna_to_protein("AUGCUGGUGUAGUCGUGGUUAUUCUUU")
# [1] "MLV"

rna_seq = "AUGCUGGUGUAGUCGUGGUUAUUCUUU"

# "MLV*SWLFF" will be acceptable but
#  the protein sequence should ideally end at stop codons.

'# version 1'
translate_rna_to_protein = function ( rna_seq) { 
  protein_seq = ""
  position = 1
   while (position <= nchar(rna_seq)) {
     codon = substring(rna_seq, position, position + 2)

    if (codon == "UAA" || codon == "UAG" || codon == "UGA")  
      break 
     
    for (i in 1:length(codons)) {
      if (codon == codons[i]) {
        protein_seq = paste(protein_seq, amino_acids[i], sep = "")
        break
      }
    }
     
     position = position + 3
   }
  return (protein_seq)
}

# version 2, assume we don't know code for stop codon, except "*"
translate_rna_to_protein = function ( rna_seq) { 
  protein_seq = ""
  position = 1
  while (position <= nchar(rna_seq)) {
    codon = substring(rna_seq, position, position + 2)
    
    for (i in 1:length(codons)) {
      if (codon == codons[i] ) {
        
        if (amino_acids[i] == "*") {
          codon = "STOP"
          break
        }
        
        protein_seq = paste(protein_seq, amino_acids[i], sep = "")
        break
      }
    }
    if (codon == "STOP")
      break
    
    position = position + 3
  }
  return (protein_seq)
}
## Q3 - 3pt

## Given a positive integer n,

## write a recursive function that calculates of squares of consecutive integers up to n

## series_sum(4) = 4^2 + 3^2 + 2^2 + 1^1 = 30

# Example input/output:
 series_sum(20)
# [1] 2870

series_sum <- function(n){
  if (n == 1) {
    return (1)
  } else {
    return ( n^2 + series_sum(n-1) )
  }
}

sum = c()
sum[1] = 1
series_sum = function(n) {
  if ( is.na(sum[n]) ) {
    sum[n] <<- n^2 + series_sum(n-1)
  }
  return (sum[n])
}


## Q4 - 6pt
## Let's define
## F(n) = F(n-1) + 3* F(n-2)
## F(1) = 1
## F(0) = 0
## Write a function that takes any positive integer k as input and 
## prints the value of smallest "n" such that F(n) >= k

# Example input/output:
system.time(smallest_n_finder(30) )
# [1] 6

smallest_n_finder = function ( k) { 
  n = 0
  while ( function4(n) < k ) {
    n = n + 1
  }
  print(n)
}

smallest = c()
smallest[1] = 1
function4 = function(n) {
  if (n == 0) {return (0)}
  if ( is.na(smallest[n]) ) {
    smallest[n] <<- function4(n-1) + 3*function4(n-2)
  }
  return( smallest[n] )
}

function4 = function(n) {
  if (n == 0) {
    return (0)
  } else if (n == 1) {
    return (1)
  } else {
    return ( function4(n-1) + 3*function4(n-2) )
  }
}
