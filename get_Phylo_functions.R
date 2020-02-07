library(ape)
require(XLConnect)
############################################################################################################

# A function to compare topology distances of two phylogenetic trees
Tree_TopologyDistance_Compare <- function(t_One, t_Two, method){
  
  #Comparing Topological distance between two phylogenetic trees
  return(dist.topo(t_One, t_Two, method))
}
############################################################################################################

# A function to compare topology of two phylogenetic trees
Tree_Phylogeny_Compare <- function(t_One, t_Two, Use.edge.length=F){
  
  #Do two trees represent the same phylogeny?
  #Phylogeny_match <- all.equal.phylo(t_One, t_Two, use.edge.length)
  return(all.equal.phylo(t_One, t_Two, use.edge.length=F))
}
############################################################################################################

# Defining the ConvertAlignment_Matrix function which reads an alignment of 
# given format(e.g. fasta, stockholm, clustal) and converts the alignment to a matrix
ConvertAlignment_Matrix <- function(path, fileType) {
  
  origMAlign <-readDNAMultipleAlignment(filepath = path, format = fileType)
  return (as.matrix(origMAlign))
  
}
############################################################################################################

# This function replaces base letters in a sequence matrix to 
# numbers (e.g. A=1, T=2, C=3, G=4)
Convert_letters_numbers <- function(seqMatrix){
  
  df <- seqMatrix
  return (matrix(match(df, LETTERS[c(1,3,7,20)]), nrow(df), dimnames=dimnames(df)))
  
}
#############################################################################################################

# Defining "Future_Changes" function to replace bases with following logic;
# If A AND previous nucleotide is G, then 1; if A and  previous is A, T, C  then 2; 
# if not A then 3 bases in first column is labelled as 4.

replace_base <- function(m){
  df <- m
  for (i in 1:nrow(df)){
    #conditional replacement
    if (df[i,2] == "A" && df[i,1] == "G" ){df[i,2] <- 1}
    else if (df[i,2] == "A" && (df[i,1] == "C" || df[i,1] == "A" || df[i,1] == "T")){df[i,2] <- 2}
    else if (df[i,2] != "A") {df[i,2] <- 3}
  }
  out<- df[,2]
  return (out)
}

Future_Changes <- function(seqMatrix){
  df <- seqMatrix
  # pass arguments to replace_base function (arg - two columns at a time from left to right) to replace bases   according to the logic mentioned above
  df_temp <- sapply(ncol(df):2, function(x) {replace_base(df[,(x-1):x])}, simplify = FALSE)
  # Convert df_temp list object to a dataframe
  out_temp <- data.frame(matrix(unlist(df_temp), nrow=length(df_temp[[1]]), byrow =  FALSE),stringsAsFactors=FALSE)
  out_temp <- out_temp[,ncol(out_temp):1]
  
  # Bases in first column is replaced with 4
  first_col <- as.data.frame(df[,1])
  first_col[,1] <- 4
  
  out <- cbind(first_col[,1],out_temp)
  row.names(out) <- row.names(df)
  colnames(out) <- c(1:ncol(df))
  return(out)
}
##############################################################################################################