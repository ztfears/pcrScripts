#This script takes an multiple sequence alignment (in fasta format) as input and exports a spreadsheet with a ranked list of potential CRISPR-Cas guide RNAs.
library(seqinr)

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}

input_file <- args[1]

z <- args[2] #number of k-mers to find
len <- args[3]


  aln <- read.alignment(input_file, format = "fasta", forceToLower = TRUE)
  l <- as.numeric(getLength(aln$seq[1]))#get length of sequence alignment
  n <- aln$nb #get number of sequences
  con <- consensus(aln, method = "majority")#get consensus sequence
  
  #convert input to a character matrix
  aln1 <- as.matrix(aln)
  #count number of nucleotide characters
  aSum <- colSums(aln1=='a', na.rm=TRUE)
  cSum <- colSums(aln1=='c', na.rm=TRUE)
  gSum <- colSums(aln1=='g', na.rm=TRUE)
  tSum <- colSums(aln1=='t', na.rm=TRUE)
  gapSum <- colSums(aln1=='-', na.rm=TRUE)
  entroP <- vector(mode="numeric", length=l)#initiate vector
  
  #calculate Shannon entropy for each column of the alignment
  for(i in 1:l){
    #first find the frequency of each character
    a <- (unname(aSum[i]))/n
    c <- (unname(cSum[i]))/n
    g <- (unname(gSum[i]))/n
    t <- (unname(tSum[i]))/n
    gap <- (unname(gapSum[i]))/n
    
    #now calculate entropy at for all bases at each position
    al <- a*log2(a)
    cl <- c*log2(c)
    gl <- g*log2(g)
    tl <- t*log2(t)
    gapl <- gap*log2(gap)
    if(is.nan(al)){al<-0}
    if(is.nan(cl)){cl<-0}
    if(is.nan(gl)){gl<-0}
    if(is.nan(tl)){tl<-0}
    if(is.nan(gapl)){gapl<-0}
    
    #calculate the Entropy with either of 2 algorithms
    #entroP[i] <- (log2(5)-(-sum(al,cl,gl,tl,gapl))) #(per Crooks et al, 2004)
    entroP[i] <- (-100*sum(al,cl,gl,tl,gapl)) #(per ViPRbrc.org)
  }
  ##write entropy results for each position to file
  entroP1 <- as.data.frame(entroP,drop = FALSE)
  entropy_out <- paste0(input_file,"-AllEntropyValues.csv")
  write.csv(file=entropy_out,entroP1, row.names=TRUE)
  
  #iterate over results to find 20 adjacent nucleotides with lowest entropy
  m <- l-len
  entroP20 <- vector(mode="numeric", length=m)#initiate vector to hold summed values for each 20-mer
  for(k in 1:m){
    entroP20[k] <- sum(entroP[k],entroP[k+19],entroP[k+18],entroP[k+17],entroP[k+16],entroP[k+15],entroP[k+14],entroP[k+13],entroP[k+12],entroP[k+11],entroP[k+10],entroP[k+9],entroP[k+8],entroP[k+7],entroP[k+6],entroP[k+5],entroP[k+4],entroP[k+3],entroP[k+2],entroP[k+1])
  }
  
  #find the lowest 20 entropy values (i.e. the most conserved 20-mer regions)
  best <- order(entroP20)[1:z]
  
  #initialize data fram
  results.df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Entropy_Score","Rank","Sequence","Aligned_Start_Position","Aligned_End_Position"))
  
  for(q in 1:z){
    #print(q)
    #get entropy value for each of the best 20mers
    aa <- entroP20[best[q]]
    #get rank for each 20mer
    bb <- as.numeric(q)
    #get sequence data for each of the best 20mers
    cc <- paste(con[best[q]],con[best[q]+1],con[best[q]+2],con[best[q]+3],con[best[q]+4],con[best[q]+5],con[best[q]+6],con[best[q]+7],con[best[q]+8],con[best[q]+9],con[best[q]+10],con[best[q]+11],con[best[q]+12],con[best[q]+13],con[best[q]+14],con[best[q]+15],con[best[q]+16],con[best[q]+17],con[best[q]+18],con[best[q]+19], sep="")
    #get aligned start position for each of the best 20mers
    dd <- as.numeric(best[q])
    #get aligned end position for each of the best 20mers
    ee <- (best[q])+19
    if(cc %in% "-") next
    addRow <- data.frame(Entropy_Score=aa,Rank=bb,Sequence=cc,Aligned_Start_Position=dd,Aligned_End_Position=ee)
    results.df <- rbind(results.df,addRow)
  }
  #save results to csv
  tempString <- c(input_file,"_",z,"-best_",len,"mers",".csv")
  outfile <- paste(tempString, collapse="")
  write.csv(file=outfile,results.df, row.names=FALSE)
