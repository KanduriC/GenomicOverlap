# GenomicOverlap

##### Author: Chakravarthi Kanduri

### compute_overlap function

This R function computes the overlap of the regions between two files containing genomic intervals. 

### Usage 

`compute_overlap (fileA, fileB)`

### Arguments

- fileA, fileB:        Tab-delimited files with genomic intervals sorted in ascending order.  
- genomeSize:          The size of the reference genome in bases; numeric; defaults to 10000000.  
- startCol, endCol:    The column numbers of start and end coordinates in files; numeric; defaults to 1 and 2.  

### Examples

`source("computeOverlap.R")`  
`compute_overlap("A.txt","B.txt")`  
`compute_overlap(newfile1,newfile2,genomeSize=20000000,startCol=2,endCol=3)`  
