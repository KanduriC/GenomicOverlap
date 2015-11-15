### Author: Chakravarthi Kanduri

### compute_overlap function: To compute the overlap of the regions between two
### files with genomic intervals 
### Usage: compute_overlap(fileA, fileB)

### Arguments: fileA, fileB: tab-delimited files with genomic intervals sorted in
### ascending order. genomeSize: The size of the reference genome in bases; numeric;
### defaults to 10000000. startCol, endCol: The column numbers of start and end
### coordinates in files; numeric; defaults to 1 and 2

compute_overlap <- function(fileA, fileB, genomeSize = 1e+07, startCol = 1, endCol = 2) {
    
    startTime <- Sys.time()
    message("Computing overlap ... ...")
    # read both files containing genomic intervals
    file1 <- read.table(fileA, header = F, sep = "\t")
    file2 <- read.table(fileB, header = F, sep = "\t")
    
    # create a new array (vector) of zeros with size (length) equal to reference
    # genome's size For each interval (x,y) in fileB score_array[x]=1 and to ensure
    # end-exclusiveness, score_array[y-1]= -1
    score_array <- rep(0, genomeSize)
    score_array[file2[, startCol]] = 1
    score_array[file2[, endCol] - 1] = -1
    
    # Compute a cumulative sum vector from score_array (for each i in score_array:
    # cumulative_score[i]=score_array[i]+cumulative_score[i-1])
    cum_score <- vector()
    length(cum_score) = length(score_array)
    cum_score[1] = score_array[1]
    for (i in 2:length(score_array)) {
        cum_score[i] = score_array[i] + cum_score[i - 1]
    }
	# Adjust the overlap size of each interval by 1, to accommodate end-exclusiveness
    for (j in length(cum_score):2) {
        if (cum_score[j - 1] - cum_score[j] > 0) {
            cum_score[j] = 1
        }
    }
    # Given a genomic interval, this function computes its overlap with reference
    # genome
    overlap <- function(interval) {
        ind1 <- interval[1]
        ind2 <- interval[2] - 1
        return(sum(cum_score[ind1:ind2]))
    }
    # Compute the overlap of all the genomic intervals in fileA 
    # Return the cumulative sum of overlapping bases
    cumSum <- apply(file1, 1, overlap)
    endTime <- Sys.time()
    Time.Taken <- endTime-startTime
    cumSum <- cbind(sum(cumSum),Time.Taken)
    colnames(cumSum)=c("Overlap (bases)","Computation time (Secs)")
    return(cumSum)
} 
