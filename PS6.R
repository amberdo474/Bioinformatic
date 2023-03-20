## PS6 - Databases



## Q1 - 6pt Complete the assigned course in datacamp about importing data into R from the web.



## Optional - If you would like to learn more, I also recommend the data.table package. Specifically, the function "fread". 

## Q2 - 9pt 
## Identify one database that you find interesting/useful. You can use one of the options that I posted on the course website.
## Describe how the data is organized, how the data can be accessed (web interface, downloads, API?). 
## Write another paragraph about a question that can be answered using the information in the database.
## Provide screenshots as needed to justify your answers. 

The GEO (Gene Expression Omnibus) database if very interesting. This database uses as Gene expression data repository. 
It stores high-throughput gene expression and functional genomics data. 
The data in GEO is organized into series and samples, with each series containing a group of related samples.
A sample is a single experiment or assay, while a series is a set of related samples that share a common biological question or hypothesis.

The GEO database can be accessed through a web interface, which allows users to search, browse, 
and download the data. The search function of the web interface allows users to search for data 
using keywords, such as gene names, diseases, or authors. Users can also browse the data by organism, 
platform, or study type.Ex: doc is included. 

Once a user has identified a dataset of interest, 
they can download the data in various formats, including GEO SOFT format, 
which is a tab-delimited text format that contains both sample and series data. 
The one i download here is txt file. 


The question that can be answered using the information in this database is :
the expression of Axin2 gene through all 6 samples: using the identifier for a 
specific gene identifier, we can see the top 10 gene identifier that is highly express 
compared to the rest in this specific database. So there, we can have a clue that not just Axin2 gene, 
but some other related gene that are helpful toward stem cell development in the brain. (For mouse) 


## Q3 - 6pt 
## Import data from your database of choice into R.  You can use import functions you learned in the datacamp course (Q1) or any other tools. 
## Produce a plot using the data. Write a few sentences about how you would interpret these plots. 
path3 = '~/Downloads/GSE178509_bulk_TPM.txt'
mouse_sequence = read.delim(path3, stringsAsFactors = FALSE)
genes = c()
length_gene = length(mouse_sequence[, 2]) 
for (i in 1: length_gene ) {
  genes = c(genes, mean(c(mouse_sequence$SampleR_1[i], mouse_sequence$SampleR_2[i], mouse_sequence$SampleR_3[i],
                          mouse_sequence$SampleA_4[i], mouse_sequence$SampleA_5[i], mouse_sequence$SampleA_6[i] )) )
}

hist(genes)
hist (log10( genes + 1 ) )  

## this show that a histogram of Skewed data, which mean 
#- few highly expressed genes that are driving the overall expression pattern 
# - can identify patterns and trends in the data and to generate hypotheses about the biological 
# processes or conditions that underlie these expression patterns.


