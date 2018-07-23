# STATISTICS
# Load data
library(dplyr)
library(rafalib)
cell_type = "mutYmutM" #Enter the cell type for analysis
full_data = read.csv("170421_analysis_04.csv")
cell_data =  filter(full_data, cell == cell_type)

# Population data
cell_All_GtoT <- select(cell_data, T) %>% unlist()

# Separate data based on the starting base
cell_Astart_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),1,1)=="A") %>%
  select(T) %>% unlist()
cell_Cstart_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),1,1)=="C") %>%
  select(T) %>% unlist()
cell_Gstart_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),1,1)=="G") %>%
  select(T) %>% unlist()
cell_Tstart_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),1,1)=="T") %>%
  select(T) %>% unlist()

#Separate data based on the ending base
cell_Aend_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),3,3)=="A") %>%
  select(T) %>% unlist()
cell_Cend_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),3,3)=="C") %>%
  select(T) %>% unlist()
cell_Gend_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),3,3)=="G") %>%
  select(T) %>% unlist()
cell_Tend_GtoT = filter(cell_data, substring(as.character(cell_data$lesion),3,3)=="T") %>%
  select(T) %>% unlist()

# Computing t-test p-values
# mutY
print(c("The average G to T mutation level for all sequence contexts is:", mean(cell_All_GtoT)))

start_stats <- matrix(0, nrow = 4, ncol = 2)
rownames(start_stats) <- c("A", "C", "G", "T")
colnames(start_stats) <- c("mean", "p-value")

start_stats[,1] <- c(mean(cell_Astart_GtoT), mean(cell_Cstart_GtoT), mean(cell_Gstart_GtoT), mean(cell_Tstart_GtoT))
start_stats[,2] <- c(t.test(cell_All_GtoT, cell_Astart_GtoT)$p.value,
                     t.test(cell_All_GtoT, cell_Cstart_GtoT)$p.value,
                     t.test(cell_All_GtoT, cell_Gstart_GtoT)$p.value,
                     t.test(cell_All_GtoT, cell_Tstart_GtoT)$p.value)

end_stats <- matrix(0, nrow = 4, ncol = 2)
rownames(end_stats) <- c("A", "C", "G", "T")
colnames(end_stats) <- c("mean", "p-value")

end_stats[,1] <- c(mean(cell_Aend_GtoT), mean(cell_Cend_GtoT), mean(cell_Gend_GtoT), mean(cell_Tend_GtoT))
end_stats[,2] <- c(t.test(cell_All_GtoT, cell_Aend_GtoT)$p.value,
                   t.test(cell_All_GtoT, cell_Cend_GtoT)$p.value,
                   t.test(cell_All_GtoT, cell_Gend_GtoT)$p.value,
                   t.test(cell_All_GtoT, cell_Tend_GtoT)$p.value)

# PLOTTING
cell_type = "mutYmutM" #Enter the cell type for analysis
full_data = read.csv("170421_analysis_04.csv")
cell_data =  filter(full_data, cell == cell_type)

cell_All_GtoT <- select(cell_data, c(lesion, T))

# Compute the mutation statistics for all sequence contexts
contexts = factor(c("AXA", "AXC", "AXG", "AXT", "CXA", "CXC", "CXG", "CXT", 
                    "GXA", "GXC", "GXG", "GXT", "TXA", "TXC", "TXG", "TXT"))
cell_All_GtoT_stats <- data.frame(matrix(0, nrow = 16, ncol = 3))
colnames(cell_All_GtoT_stats) <- c("lesion", "mean", "SD")
rownames(cell_All_GtoT_stats) <- contexts
for(i in contexts){
  values = filter(cell_data, lesion==i) %>% select(T) %>% unlist
  cell_All_GtoT_stats[i, "lesion"] = i
  cell_All_GtoT_stats[i, "mean"] = mean(values)
  cell_All_GtoT_stats[i, "SD"] = sd(values)
}

# Plot based on the 5' base
mypar(1,2)
boxplot(cell_All_GtoT$T ~ as.factor(substr(as.character(cell_All_GtoT$lesion),1,1)), xlab="5' base to 8-oxoG", ylab="G to T mutation level (%)", cex.lab=0.8, cex.axis=0.7) #boxplot
stripchart(cell_All_GtoT_stats$mean ~ as.factor(substr(as.character(cell_All_GtoT_stats$lesion),1,1)), vertical=TRUE, pch=16, add=TRUE) #add points
text(as.factor(substr(as.character(cell_All_GtoT_stats$lesion),1,1)), cell_All_GtoT_stats$mean, label=cell_All_GtoT_stats$lesion, pos=4, cex=0.6)

# Plot based on the 3' base
boxplot(cell_All_GtoT$T ~ as.factor(substr(as.character(cell_All_GtoT$lesion),3,3)), xlab="3' base to 8-oxoG", ylab="G to T mutation level (%)", cex.lab=0.8, cex.axis=0.7) #boxplot
stripchart(cell_All_GtoT_stats$mean ~ as.factor(substr(as.character(cell_All_GtoT_stats$lesion),3,3)), vertical=TRUE, pch=16, add=TRUE) #add points
text(as.factor(substr(as.character(cell_All_GtoT_stats$lesion),3,3)), cell_All_GtoT_stats$mean, label=cell_All_GtoT_stats$lesion, pos=4, cex=0.6)