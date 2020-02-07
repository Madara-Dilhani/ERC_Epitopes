#######################################################################################################
# Load required packages
library("reshape2", quietly=TRUE)
library("ape", quietly=TRUE)
library("geiger", quietly=TRUE)
library("igraph", quietly=TRUE)
library("plyr", quietly = TRUE)
source("get_Phylo_Indep_Divergent.R")

#######################################################################################################

divergent <-  data.frame()

for(i in 1:10){
  tryCatch({
    tree_pol <- read.tree(file = paste(getwd(), "/Data/4/", i,"/Pol-out-", i, ".nwk", sep=""))
    tree_pol_consensus <- read.tree(file = paste("Data/4/", i,"/Pol-out-", i, "_consensus.nwk", sep=""))
    
    dist = dist.nodes(x = tree_pol)
    bigtree.mppd <- pdist.clusttree(tree = tree_pol, distmat= dist)
    bigtree.support <- as.numeric(tree_pol_consensus$node.label)
    bigtree.support[is.na(bigtree.support)] <- 0 
    clustering <- prosperi.cluster(tree= tree_pol,thresh= 0.06,distvec=bigtree.mppd,
                                   rthresh=0.70,reliabilityvec=bigtree.support)  
    
    data <- data.frame(matrix(ncol = 4, nrow = 999))
    data[,2] <- clustering$membership
    data[,1]<- 1:999
    
    is_tip <- tree_pol$edge[,2] <= length(tree_pol$tip.label)
    ordered_tips <- tree_pol$edge[is_tip, 2]
    tips <- tree_pol$tip.label[ordered_tips]
    
    data[1:500,3]<- tips
    data[,4] <- sapply(1:999, function(x) {unlist(strsplit(data[x,3], split = '[.]'))[2]})
    colnames(data) <- c("Nodes", "Cluster", "Sequence", "Subtype")
    
    upper<-upper.tri(dist[1:500,1:500], diag=F) #turn into a upper triangle
    dist_table<-dist[1:500, 1:500] #take a copy of the original cor-mat
    dist_table[!upper]<-NA#set everything not in upper triangle o NA
    BL_table<-na.omit(melt(dist_table, value.name ="BL")) #use melt to reshape the matrix into triplets, na.omit to get rid of the NA rows
    colnames(BL_table)<-c("Seq1", "Seq2", "Distance")
    
    attach(BL_table)
    BL_table <- BL_table[order(Seq1, Seq2),]
    mege_BL_seq <- merge(BL_table,data, by.x = "Seq1", by.y = "Nodes")
    mege_BL_seq <- merge(mege_BL_seq,data, by.x = "Seq2", by.y = "Nodes")
    attach(mege_BL_seq)
    mege_BL_seq <- mege_BL_seq[order(Seq1, Seq2),]
    mege_BL_seq$Samp <- NA
    mege_BL_seq$Samp <- i
    
    # Count Frequency of each factor level
    frequency <- count(data[1:500,],vars = c("Cluster","Subtype"))
    No_Clusters <- data.frame()
    No_Clusters <- count(unique(frequency[,1]))
    track <- (nrow(No_Clusters)-1)
    
    for( c in 1: track){
      
      clust <- subset(mege_BL_seq, mege_BL_seq[,4]== c & mege_BL_seq[,7]== c)
      div_pair <- data.frame(matrix(ncol = 7, nrow = nrow(clust)))
      div_pair[,1] <- clust[,5]
      div_pair[,2] <- clust[,8]
      div_pair[,3] <- clust[,3]
      div_pair[,6] <- clust[,4]
      div_pair[,7] <- clust[,10]
      
      Ntip <- length(tree_pol$tip.label)
      for (i in 1: nrow(div_pair)){
        
        print(i)
        node_ancestor <-  FindAncestor(c(as.character(div_pair[i,1]), as.character(div_pair[i,2])), tree_pol)
        div_pair[i,4] <- bigtree.support[node_ancestor-Ntip] 
        
        # Annotating seq pairs with bootstrap support >= x% as divergent (*change this)
        if (div_pair[i,4] >= 0.7 && div_pair[i,3] >= 0.01){
          div_pair[i,5] <- 1 
        } else {  div_pair[i,5] <- 0 }
      }
      
      div_pair_final <- subset(div_pair, (div_pair[,5] == 1))
      most_divergent <- div_pair_final[which(div_pair_final[,3] ==  max(div_pair_final[,3])), ]
      divergent <-  rbind(divergent, most_divergent)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}

write.csv(divergent, "Output/4/divergent_pairs_from_clusters.csv")







