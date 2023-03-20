## BCH 339N
## PS5 - Hidden Markov Models

## We have covered the Viterbi and Forward Algorithms in class. 
## These methods are implemented in the R package - HMM

install.packages("HMM")
library(HMM)
?initHMM
?forward
?viterbi 

## In this homework, we are going to develop an HMM that can predict secondary structure of proteins.
## Kaggle (https://www.kaggle.com) organizes data-science and machine learning challenges for the community
## They have many challenges related to bioinformatics. 
## We will explore the https://www.kaggle.com/alfrandom/protein-secondary-structure/data challenge. 
## Read the associated introduction to familiarize yourself with this dataset. 

## Q1 -6pt
## What is protein secondary structure? 
  #Protein secondary structure is the spatial confornmation of the polypeptide backbone, is bonded through hydrogen bond between polypeptide. 
  #The 2 most common forms are alpha helices and beta sheets

## We are going to work with the three state secondary structure encoding (sst3). 
    
## Explain the biochemical features of the three  that we will model (C, E, and H ).
    #Q3 - sst3: the three-state (Q3) secondary structure
    #biochemical features of the three states : (three (Q3) by merging (E, B) into E, (H, G, I) into E, and (C, S, T) into C.)
    # C - loops, E - beta-strand, H- alpha-helix
#C- C-loops are short. They are typically found on the surface of proteins and are important for protein-protein interactions and binding to ligands.
#E- Beta sheets are secondary structure elements that are formed by hydrogen bonding between adjacent strands of polypeptide chains. E-beta sheets are a type of beta sheet that are parallel, meaning that the strands run in the same direction. They are often found in the cores of proteins and are important for stabilizing the protein structure.
#H-Alpha helices are another type of secondary structure element that are formed by hydrogen bonding between adjacent amino acid residues within a polypeptide chain. H-alpha helices are a type of alpha helix that have a right-handed coil and are often found on the surface of proteins. They are important for protein-protein interactions and for forming channels through which small molecules can pass.




  
## Q2 -6pt

## We will use the data with the strict quality controls for this homework. Here Download Here is the file for convenience. 

## You will have to change the next line of code to specify the folder where the file is. 
## First, we will remove all examples with non-standard amino acids 
#path = file.path("~", "datasets", "states.csv")
# "~/datasets/states.csv"
prot_sec_data = read.csv("~/Downloads/2018-06-06-pdb-intersect-pisces.csv", stringsAsFactors = F)

str(prot_sec_data)
## only get the standard_ amino acids
prot_sec_data_std = prot_sec_data[as.logical(prot_sec_data$has_nonstd_aa) == F, ] 



## Next, we are going to split our dataset into two parts. 
## The first part should contain 10% of the sequences and will be used for training our HMM
## The remaining 90% will be used as a test set to assess how well our HMM does in predicting secondary structure. 
## To ensure reproducibility of your code, we are going to use the set.seed function. 
?set.seed
set.seed(3)

train  = sample(nrow(prot_sec_data_std), size = floor(nrow(prot_sec_data_std)*.1 ) )
test =  setdiff(1:nrow(prot_sec_data_std), train  )

train_data = prot_sec_data_std[train, ]
test_data = prot_sec_data_std[test, ]

##  We will use the training data to infer the parameters of the HMM model. 
## Specifically, use sst3 and seq columns of the data and determine the 
## transition, emission and initial probabilities.
## Write a function that will take the training data as input. 
## The output should be a list with names
## c(“initial_probs”, “emission_probs”, “transition_probs”)
## Note that emission_probs and transition_probs should be matrices of the form defined in initHMM. 
symbols = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") 
states =  c("C", "E", "H")



get_hmm_probs<- function(train_data, symbols, states){
  # the initial state , initialProb
  countC = 0
  countE = 0
  countH = 0
  
  
  for (i in 1: length(train_data$sst3)) {
    Q3 = substr(train_data$sst3[i], 1, 1)
    if (Q3 == "C") { countC = countC + 1 }
    else if (Q3 == "E") { countE = countE + 1}
    else { countH = countH + 1 }
  }
  initial_probs = c(countC, countE, countH)
  initial_probs = initial_probs / length(train_data$sst3)
  
  #transition matrix
  C = c(0,0,0)# where previous state was C 
  E = c(0,0,0) # where previous state was E
  H = c(0,0,0)  # where previous state was H

  for (i in 1: length(train_data$sst3)) {
    # archive each line
    line = train_data$sst3[i]

    for (j in 2: nchar(line) ) {
      
      current = substr(line, j, j)
      previous = substr(line, j -1 , j -1 )
      
      
      if (previous == "C" ) {
        if (current == "C") { C[1] = C[1] + 1 }
        else if (current == "E") { C[2] = C[2] + 1 }
        else { C[3] = C[3] + 1 }
      } else if (previous == "E") {
        if (current == "C") { E[1] = E[1] + 1 }
        else if (current == "E") { E[2] = E[2] + 1}
        else { E[3] = E[3] + 1 } 
      } else {
        if (current == "C") { H[1] = H[1] + 1 }
        else if (current == "E") {  H[2] = H[2] + 1 }
        else {  H[3] = H[3] + 1 }
      }
    }
  }
  
  transition_probs = matrix( c( C / sum(C) ,E / sum(E) , H / sum(H) ), nrow = 3, ncol = 3, byrow = TRUE )
  
  
  
  #emission matrix
  # 60 track , 3 for each is 20 observation
  Cs = rep(0, 20)
  Es = rep(0, 20)
  Hs = rep(0, 20)

  for (i in 1: length(train_data$sst3)) {
    # grab both line seq and sst3 of the same sample i 
    line_sst3 = train_data$sst3[i]
    line_seq = train_data$seq[i]
    
    for (j in 1: nchar(line_seq) ) {
      # grab each char of these 2 sequences to compare
      sst3_current = substr(line_sst3, j, j)
      seq_current = substr(line_seq, j, j)
      
      # check for states first, then check for observation 
      if (sst3_current == "C")  {
        
        # match with the symbol in vector symbols . 
        # which(my_vector == match_element) , give me the index to change my element in rows (C)
        index = which(symbols == seq_current) 
        Cs[index] = Cs[index] + 1

        
      } else if (sst3_current == "E") {
        index = which(symbols == seq_current) 
        Es[index] = Es[index] + 1
        
      } else {
        index = which(symbols == seq_current) 
        Hs[index] = Hs[index] + 1
      }
      
    }
  }
  
  emission_probs = matrix( c( Cs / sum(Cs) , Es / sum(Es) , Hs / sum(Hs) ), nrow = 3, byrow = TRUE )
  
  return(list(initial_probs = initial_probs, "emission_probs" = emission_probs,  "transition_probs" = transition_probs)) 
}


hmm_values =  get_hmm_probs(train_data, symbols, states)
hmm_values

## The expected output should be as follows. 
## Do not worry if yours is slightly different. 
# $initial_probs
# [1] 1 0 0
# 
# $emission_probs
# [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]
# [1,] 0.06585424 0.01165331 0.07881398 0.05734681 0.03039266 0.11675254 0.03626634
# [2,] 0.06147861 0.01684346 0.03523090 0.04828457 0.05855106 0.04896633 0.02452327
# [3,] 0.11237240 0.01019533 0.05272974 0.09006702 0.04173868 0.03504955 0.02049012
# [,8]       [,9]      [,10]      [,11]      [,12]      [,13]      [,14]
# [1,] 0.03188721 0.05509976 0.06162143 0.02146717 0.05971928 0.07833321 0.03431193
# [2,] 0.09332077 0.04435443 0.10398829 0.02221732 0.02648834 0.02263841 0.02845341
# [3,] 0.05793930 0.06470303 0.11973293 0.02671922 0.03325915 0.02469259 0.04735854
# [,15]      [,16]      [,17]      [,18]      [,19]      [,20]
# [1,] 0.04505597 0.07720446 0.05752448 0.04383315 0.01026327 0.02659880
# [2,] 0.04784344 0.05181368 0.06474705 0.12821078 0.01862806 0.05341782
# [3,] 0.05923237 0.04983277 0.04145271 0.05986647 0.01582762 0.03674048
# 
# $transition_probs
# [,1]        [,2]       [,3]
# [1,] 0.80925703 0.110084193 0.08065877
# [2,] 0.20296365 0.783140502 0.01389585
# [3,] 0.09893198 0.004737097 0.89633093


## Q3 - 6pt
## Look at how hmms are defined in the HMM package
## We will Use the inferred parameters from Q3 to define an HMM. 
?initHMM

sec_struct_hmm = initHMM(States = states, Symbols = symbols,
                         startProbs = hmm_values$initial_probs, 
                         transProbs = hmm_values$transition_probs, 
                         emissionProbs = hmm_values$emission_probs)


## Next, We are going to assess the performance on the test data
## For each example in the test data, use the given sequence (seq column) to predict the most likely path of hidden states. 
## You can use the appropriate function in the hmm package for this step. 
## For each example, compare the predicted most likely hidden state path with  the experimentally identified values (sst3 column). 
## Output a vector named percent_correct containing the percentage of amino acids whose secondary structure was correctly predicted for each example in test data.

percent_correct = c()
lessThan1 = c()
predicted_States = c()

## YOUR CODE HERE
 # seq is the observation, and sst3 is the hidden states
 # predict sst3 using the observation, and compare that to real sst3 outcome
for (i in 1: length(test_data$seq)) {
  
  line_seq = unlist(strsplit(test_data$seq[i], ""))
  actual_sst3 = unlist(strsplit(test_data$sst3[i], ""))
  
  predicted_sst3 = viterbi(sec_struct_hmm, line_seq)
  print(predicted_sst3)
  # check each letter for correct number of predicted sst3 vs real outcome
  correct = 0
  total = length(actual_sst3)
  
  for (index in 1: total ) {
    if ( actual_sst3[index] == predicted_sst3[index] ) { correct = correct + 1 }
  }
  
  ratio = correct/total * 100
  if ( ratio < 1 ) { 
    lessThan1 = c(lessThan1, test_data$pdb_id[i]) 
    predicted_States = c(predicted_States, paste(predicted_sst3, collapse = ""))
    
  }
  percent_correct = c(percent_correct, ratio)
  
}

## Plot the distribution of percent_correct 
hist(percent_correct, 40)

## Determine the mean, median and 90th percentile of this distribution 
mean(percent_correct)
median(percent_correct)
quantile(percent_correct, .90)
#[1] 42.40327
#[1] 41.62304
#90% 
#53.94821 

###Write a sentence about your interpretation of these results.
# the predicted state of secondary structure (sst3), on average, roughly only 40% correct 
# out of the actual outcome of the train_data sst3, given the probabilty was drawn from test_data 
# the viterbi algorithm to predict secondary structure is not very accurate. 


## Q4 - 3pt 
## Identify the training examples where we succeed less than 1%. Output the pdb_ids of these examples
## Examine the predicted hidden state for these and the actual secondary structure. 
## You will notice that these test examples have a high percentage of "E" in their secondary structure. 


## Explain why you think our HMM doesn't work for this certain class of proteins?
## Hint: Think about the biology of the "E" state and the definition of a Markov chain

for (i in 1: length(lessThan1)) {
  # print out actual secondary structure 
  print(test_data$sst3[test_data$pdb_id == lessThan1[i]])
  # the predicted structure
  print(predicted_States[i])
}

## Explain why you think our HMM doesn't work for this certain class of proteins
## Hint: Think about the biology of the "E" state and the definition of a Markov chain

#Complex sequence patterns: HMMs are particularly good at recognizing simple sequence patterns, 
#such as repeated C (loops). However, if a particular class of proteins contains complex sequence patterns 
#that are difficult to model, an HMM may not be effective, and E is (beta-sheet) so it could be more complex and unique, 
# Markov chain is bias with various probability, so it will be difficult for it to stay in E state, or switch to E state

Using the professor data for initial, emission, and transition matrix. Identified, 3 of the samples with less than 1% succes. 


initial = c(1, 0, 0)
emission = matrix( c(0.06585424, 0.01165331, 0.07881398, 0.05734681, 0.03039266, 0.11675254, 0.03626634,
                     0.03188721, 0.05509976, 0.06162143, 0.02146717, 0.05971928, 0.07833321, 0.03431193,
                     0.04505597, 0.07720446, 0.05752448, 0.04383315, 0.01026327, 0.02659880,
                     
                     0.06147861, 0.01684346, 0.03523090, 0.04828457, 0.05855106, 0.04896633, 0.02452327,
                     0.09332077, 0.04435443, 0.10398829, 0.02221732, 0.02648834, 0.02263841, 0.02845341,
                     0.04784344, 0.05181368, 0.06474705, 0.12821078, 0.01862806, 0.05341782,
                     
                     
                     0.11237240, 0.01019533, 0.05272974, 0.09006702, 0.04173868, 0.03504955, 0.02049012,
                     0.05793930, 0.06470303, 0.11973293, 0.02671922, 0.03325915, 0.02469259, 0.04735854,
                     0.05923237, 0.04983277, 0.04145271, 0.05986647, 0.01582762, 0.03674048) , nrow = 3, byrow = TRUE)

transition = matrix( c(0.80925703, 0.110084193, 0.08065877,
                       0.20296365, 0.783140502, 0.01389585,
                       0.09893198, 0.004737097, 0.89633093), nrow = 3, byrow = TRUE)


struct_hmm = initHMM(States = states, Symbols = symbols,
                     startProbs = initial, 
                     transProbs = transition, 
                     emissionProbs = emission)

percent_correct_2 = c()
Es = 0
lessThan1 = c()
predicted_States = c()

for (i in 1: length(test_data$seq)) {
  
  line_seq = unlist(strsplit(test_data$seq[i], ""))
  actual_sst3 = unlist(strsplit(test_data$sst3[i], ""))
  
  predicted_sst3 = viterbi(struct_hmm, line_seq)
  # check each letter for correct number of predicted sst3 vs real outcome
  correct = 0
  total = length(actual_sst3)
  
  for (index in 1: total ) {
    if ( actual_sst3[index] == predicted_sst3[index] ) { correct = correct + 1 }
    if ( predicted_sst3[index] == 'E') { Es = Es + 1 }
  }
  
  ratio = correct/total * 100
  if ( ratio < 1 ) { 
    lessThan1 = c(lessThan1, test_data$pdb_id[i]) 
    predicted_States = c(predicted_States, paste(predicted_sst3, collapse = ""))
    
  }
  percent_correct_2 = c(percent_correct_2, ratio)
  
}

hist(percent_correct_2, 40)

## Determine the mean, median and 90th percentile of this distribution 
mean(percent_correct_2)
median(percent_correct_2)
quantile(percent_correct_2, .90)


## Q4 - 3pt 
## Identify the training examples where we succeed less than 1%. Output the pdb_ids of these examples
## Examine the predicted hidden state for these and the actual secondary structure. 
## You will notice that these test examples have a high percentage of "E" in their secondary structure. 
for (i in 1: length(lessThan1)) {
  # print out actual secondary structure 
  print(test_data$sst3[test_data$pdb_id == lessThan1[i]])
  # the predicted structure
  print(predicted_States[i])
}

## Explain why you think our HMM doesn't work for this certain class of proteins
## Hint: Think about the biology of the "E" state and the definition of a Markov chain

#Complex sequence patterns: HMMs are particularly good at recognizing simple sequence patterns, 
#such as repeated C (loops). However, if a particular class of proteins contains complex sequence patterns 
#that are difficult to model, an HMM may not be effective, and E is (beta-sheet) so it could be more complex and unique, 
# Markov chain is bias with various probability, so it will be difficult for it to stay in E state, or switch to E state



