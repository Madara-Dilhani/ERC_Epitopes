#####################################################################################################
library("ppcor")
library("stats")
library("dplyr")
#####################################################################################################

# Generate comparison vector
in_lv <- paste0("IN", 1:13)
rt_lv <- paste0("RT", 1:27)
combo1 <- expand.grid(in_lv, rt_lv, stringsAsFactors=FALSE)


f_Corr <- function(dataname, x){
  cor.test(x = dataname[, which(names(dataname) %in% combo1[x,1])], 
            y = dataname[, which(names(dataname) %in% combo1[x,2])], 
            method = "pearson")$estimate
}

f_Corr_PVal <- function(dataname, x){
  cor.test(x = dataname[, which(names(dataname) %in% combo1[x,1])], 
            y = dataname[, which(names(dataname) %in% combo1[x,2])], 
            method = "pearson")$p.value
}


#read branch length values from excel file
filename= "Data/5/DivPair_Distances_filtered_ls_1_Subtype_New.csv"
data_BL = read.csv(filename, header = TRUE, sep = ",")

## Calculate correlation: This will take a while if you run it
pcor.vec <- sapply(1:nrow(combo1), function(x) {f_Corr(data_BL, x)})
pcor.pval.vec <- sapply(1:nrow(combo1), function(x) {f_Corr_PVal(data_BL, x)})
data_corr <- as.data.frame(cbind(combo1, pcor.vec, pcor.pval.vec))
colnames(data_corr) <- c("IN_Label", "RT_Label", "Correlation", "Corr_Pval")

y <- read.csv ("Data/5/Epitope_Label.csv")
mege_Label <- merge(data_corr,y, by.x = "IN_Label", by.y = "Label")
mege_Label_final <- merge(mege_Label,y, by.x = "RT_Label", by.y = "Label")
csvname = "Output/5/DivPair_Corr.csv"
write.table(x = mege_Label_final, file = csvname, append = FALSE, 
            col.names = TRUE, sep = ",", row.names = FALSE) 

#############################################################################################
