

#######################################################################
# Examine distance distribution of pol whole tree

require(ape, quietly=TRUE)

BL_distribution <- data.frame()
colnames(BL_distribution) <- c("Sample", "Median", "Tenth", "Fifteenth", "Twentieth", "TwentyFifth", "Thirtieth")
BL_Dist_All <- data.frame()

for(i in 1:112){
  tryCatch({
    
    BL_distribution[1,1] <- i
    tree_pol <- read.tree(file = paste("Data/", i,"/Pol-out-", i, ".nwk", sep=""))
    BL_distribution[1,2] <- median(dist.nodes(tree_pol)[upper.tri(dist.nodes(tree_pol), diag=FALSE)])
    x <- quantile(dist.nodes(tree_pol),seq(0.1,0.3,by=0.05))
    BL_distribution[1,3] <- x[1]
    BL_distribution[1,4] <- x[2]
    BL_distribution[1,5] <- x[3]
    BL_distribution[1,6] <- x[4]
    BL_distribution[1,7] <- x[5]
    BL_Dist_All <- rbind(BL_Dist_All, BL_distribution)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}

write.table(x = BL_Dist_All, file = "Data/4/Whole_Tree_Pol_Dist.csv", sep = ",")

########################################################################
