#PS1 , due 1/18/23. by 12. midnight



## Q1 -3pt
## Write a function that will simulate the number of tails in a series of coin tosses.
## The input to the function is the number of coin tosses
## The output should be a single number equivalent to the number of tails in the series
#input with numerical trial
#output/return of the successful times of tails (numerical, value of 1, or TRUE)

number_of_tails = function(coin_tosses) {
  total_tails <- 0
  for (iter in 1:coin_tosses) {
    random = sample(0:1)
    if (random[1] == 1)
      total_tails = total_tails + 1
  }
  return(total_tails)  
}
#number_of_tails(5)


## Q2 -3pt
## Using the function you wrote in Q2, generate 5000 experiments each with 40 coin tosses
## Plot the distribution of the outputs as a histogram
output_sample_store <- c(1:5000)
for(iter in 1:5000) {
  output_sample_store[iter] = number_of_tails(40)
  print(output_sample_store[iter])
}
hist(output_sample_store, xlab = "number of coins flip", main = "Distribution of number of coins flip across 5000 trials")

#Q3
#dna1 , #dna2, in String
##Values and objects have different data types in R. 
## You are given two DNA sequences of equal length
##Find out how many bases are different between the tw
dna1 = "AGCACAGA"
dna2 = "AGCACCCA"
#function that convert String into DNA

count_different_bases = function(dna1, dna2) {
  ## break string into vector of individual bases 
  split_dna1 <- strsplit(dna1, "")
  vector_dna1 <-unlist(split_dna1)
  
  split_dna2 <- strsplit(dna2, "")
  vector_dna2 <- unlist(split_dna2)
  
  count = 0
  for ( iter in 1:length(vector_dna1) ) {
    if ( vector_dna1[iter] != vector_dna2[iter] ) 
      count = count + 1
  }
  return (count)
}
count_different_bases(dna1, dna2)

## Q4 -3p t## A purine to purine or a pyrimidine to pyrimidine change is called a transition
## Transversion is a change from a pyrimidine to a purine or vice versa.
## As in Q3, you are given two DNA sequences of equal length
## Write a function that will output the ratio of transversions to transitions

dna1 = "AGCACAGA"
dna2 = "AGCACCCC"
transversion_transition_ratio = function (dna1, dna2) {
  split_dna1 <- strsplit(dna1, "")
  vector_dna1 <-unlist(split_dna1)
  
  split_dna2 <- strsplit(dna2, "")
  vector_dna2 <- unlist(split_dna2)
  
  
  vector_P_Y_labeled_dna1 <- c()
  vector_P_Y_labeled_dna2 <- c()
  
  transition = 0
  transversion = 0
  
  ## if A or G, is P (purine)
  for ( iter in 1:nchar(dna1)) {
    if ( vector_dna1[iter] == "A" || vector_dna1[iter] == "G" )
      vector_P_Y_labeled_dna1[iter] <- "P"
    if ( vector_dna2[iter] == "A" || vector_dna2[iter] == "G" )
      vector_P_Y_labeled_dna2[iter] <- "P"

    # if C or T, is Y (pyrimidine)
    if ( vector_dna1[iter] == "C" || vector_dna1[iter] == "T" )
      vector_P_Y_labeled_dna1[iter] <- "Y"
    if ( vector_dna2[iter] == "C" || vector_dna2[iter] == "T" )
      vector_P_Y_labeled_dna2[iter] <- "Y"
    
    # check for transition, P to P, and Y to Y, otherwise, it's a transversion.
   if ( vector_P_Y_labeled_dna1[iter] == vector_P_Y_labeled_dna2[iter] )
     transition = transition + 1
   
    else 
      transversion = transversion + 1
  }
  return (transversion/transition)
}

#transversion_transition_ratio(dna1, dna2)

## Q5 -3pt
## You are given a matrix of unknown dimensions.
## Write a function that will replace all zeros in the matrix with 1s.

## Example solution
#matrixTest <- matrix( c(0,0,0,0,0,1,1,0,3,4) , byrow = TRUE, nrow = 5)
#matrixRep <- matrixTest
replace_zeros = function(matrix) {
  
  #go through by row 1st, then column,
  for ( i in 1: nrow(matrix) ) {
    #in 1st row, with for loops go for x amount of column 
    for ( j in 1: ncol(matrix) ) {
      #change all zeros to 1
      if ( matrix[i, j] == 0 )
        matrix[ i,j ] = 1
    }
  }
  return ( matrix )
}
#replace_zeros(matrixRep)


## Q6 -3pt ## Download the dslabs package 
##(https://cran.r-project.org/web/packages/dslabs/index.html)
## execute "data(brca)" to load the BRCA data into your session
## execute "str(brca)
##Write a few sentence to describe what the output means in your own words as a comment
## Example solution install.packages("dslabs") library(dslabs) data(brca) str (brca)
install.packages("dslabs")
library(dslabs)
data(brca)
str(brca)
## the output of str(brca) provides data frame of list of 2 objects, name x and y. 
##Object x is numeric data type, and y is factor type with order "B", "M"



## Q7 -3pt
## This is an open ended question so many answers are possible.
## Look at the help page related to brca using "?brca"
## Explore the relationship between biopsy features and whether the tissue is malignant or benign
## Your answer will receive full points as long as it has at least one plot.
# Example solution
#?brca
data(brca)
table(brca$y)
table(brca$x)
dim(brca$x)
head(brca$x)
?brca

object_x = brca$x
object_y = brca$y

radius_mean = object_x[ ,1 ]
texture_mean = object_x[ ,2 ]
# just want "B" or benign
benign = sum(object_y == "B")
malignant = sum(object_y == "M")

benign = object_y[1:357]
malignant = object_y[358:569]
radius_mean_benign = radius_mean[1:357]
texture_mean_benign = texture_mean[1:357]

plot(radius_mean_benign, texture_mean_benign, xlab = "radius mean", ylab = "texture mean", main = "Benign breast masses coresponding to radius mean and texture mean")
hist(radius_mean_benign, xlab = "Radius_mean", main = "Distribution of Radius_mean across benign breast masses")

## on average, Distribution of radius_mean of mass tumor is around 12 for benign breast masses.
## on avrage, texture size around 18 and radius around 12 of breast massess are usually benign tumor.


