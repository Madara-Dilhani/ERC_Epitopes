##########################################################################################################
##### 1. Extract pairwise distances values for a given list of sister pairs
library(ape)
library(XLConnect)

out_file = "Output/1/SisPair_BL_90B.xlsx"
wb <- loadWorkbook(out_file, create = TRUE)
createSheet(wb, name = "BL")

#reading pairwise distance values from a tree into a vector and  set relevent taxa names
filenameIN = "Data/1/DENV_NJ-MCL_bootstrap-BLs.nwk"
t_IN <- read.tree(filenameIN)
t_IN$edge.length<-round(t_IN$edge.length,3)
n<-length(t_IN$tip.label)
ee_IN<-setNames(t_IN$edge.length[sapply(1:n,function(x,y) which(y==x),y=t_IN$edge[,2])],t_IN$tip.label)

#read file with sister pair Ids (Pair1, Pair2)
sisp<- read.csv("Data/1/bigtree-0.9-sis.csv", header = FALSE)   

#get branch lengths for pair1 taxa list
j<- sisp[,1]
sp_IN<- ee_IN[c(j)]

# Writing the branch lengths for pair 1 in column 1
writeWorksheet (wb, data= sisp[,1], sheet= "BL", startRow = 2, 
                startCol = 1, header = FALSE, rownames = FALSE)
writeWorksheet (wb, data= sp_IN, sheet= "BL", startRow = 2, 
                startCol = 3, header = FALSE, rownames = FALSE)
saveWorkbook(wb)

#get branch lengths for pair2 taxa list
j<- sisp[,2]
sp_IN<- ee_IN[c(j)]

# Writing the branch lengths for pair 2 in column 2
writeWorksheet (wb, data= sisp[,2], sheet= "BL", startRow = 2, 
                startCol = 2, header = FALSE, rownames = FALSe)
writeWorksheet (wb, data= sp_IN, sheet= "BL", startRow = 2, 
                startCol = 4, header = FALSE, rownames = FALSE)
saveWorkbook(wb)
##############################################################################################################

##### 2. Compare multiple phylogenetic tree datasets to find topological and distance differences
path_out= "Output/2/"
path_in= "Data/2/"
source("get_Phylo_functions.R")

#creating a excel file to write the output
excelname = paste(path_out, "Tree-Topology-Compare", ".xlsx", sep="")
wb <- loadWorkbook(excelname, create = TRUE)
createSheet(wb, name = "Tree-Comapre")

#change length=100 to no of tree comparison you do (*Change this)
compare_phylogeny <- logical(length = 5)

# iterate through each file
for (p in 1:5){
  tryCatch({
    #creating filename for Pol Trees
    filenameOne = paste(path_in, p, "/", "Pol-GTR-out-",p, ".nwk", sep="")
    filenameTwo = paste(path_in,  p, "/", "Pol-out-",p, ".nwk", sep="")
    
    #Reading Pol Trees
    t_One <- read.tree(filenameOne)
    t_Two <- read.tree(filenameTwo)
    
    #compare_distances <- Tree_TopologyDistance_Compare(t_One, t_Two, "PH85")
    compare_phylogeny[p] <- Tree_Phylogeny_Compare(t_One, t_Two)
    
    #Writing the topological distances
    writeWorksheet (wb, data= compare_phylogeny[p], sheet= "Tree-Comapre", startRow = p, 
                    startCol = 1, header = FALSE)
    # Saving the workbook
    saveWorkbook(wb)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")
    
    writeWorksheet (wb, data= "No file", sheet= "Tree-Comapre", startRow = p, startCol = 1,
                    header = FALSE)
    saveWorkbook(wb)
    
  })
}
###############################################################################################################
# 3. Reading an alignment of a given format (e.g. fasta, stockholm, clustal) and 
# converts the alignment into a matrix and then replaces base letters in a sequence matrix to 
# numbers (e.g. A=1, T=2, C=3, G=4)

library("Biostrings")
source("get_Phylo_functions.R")
filepath = "Data/3/Sequence_NA.fas"
Aln_matrix <- ConvertAlignment_Matrix(filepath, "fasta")
Aln_Matrix_Numbers <- Convert_letters_numbers(Aln_matrix)

Aln_changes <- Future_Changes(Aln_matrix)
write.csv(x = Aln_changes, file = "Output/3/DiPeptideChanges.csv")
