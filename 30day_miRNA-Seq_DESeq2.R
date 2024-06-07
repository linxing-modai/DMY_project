setwd("~/linxing2/workspace/DMY_sci_data/results/DEMs")
getwd()
raw_df <- read.table(file = "ALL_miRNA_count_2.txt",header = T,sep = "\t")
raw_df[,14] <- as.numeric(raw_df[,14])
raw_df[,13] <- as.numeric(raw_df[,13])
raw_df[,12] <- as.numeric(raw_df[,12])
raw_df[,15] <- as.numeric(raw_df[,15])
raw_df[,16] <- as.numeric(raw_df[,16])
raw_df[,17] <- as.numeric(raw_df[,17])
raw_df[,18] <- as.numeric(raw_df[,18])
raw_df[,19] <- as.numeric(raw_df[,19])
raw_df[,11] <- as.numeric(raw_df[,11])




count_df <- raw_df[,c(11:19)]
rownames(count_df) <-as.character(raw_df[,1])

colnames(count_df) <-c("30d-XX-1","30d-XX-2","30d-XX-3",
                       "30d-XYD-1","30d-XYD-2","30d-XYD-3",
                       "30d-XY-1","30d-XY-2","30d-XY-3")

count_df_XXXYD <- count_df[,c(1:6)]
count_df_XYDXY <- count_df[,c(4:9)]
count_df_XXXY <- count_df[,c(1:3,7:9)]


#---------------------------------------------------------------------------------------------------------------------
#XXXYD

#filter
count_df_XXXYD.filter <- count_df_XXXYD[rowSums(count_df_XXXYD) > 1,]


#condition table
sample_df_XXXYD <- data.frame(
  condition = c( rep("XX",3), rep("XYD",3)),
  cell_line = "Ola"
)
rownames(sample_df_XXXYD) <- colnames(count_df_XXXYD.filter)

library(DESeq2)

#directly get test result
deseq2.obj <-DESeqDataSetFromMatrix(countData = count_df_XXXYD.filter,colData = sample_df_XXXYD,design = ~condition)
deseq2.obj


deseq2.obj <-estimateSizeFactors(deseq2.obj)
sizeFactors(deseq2.obj)

#test(RLE)
deseq2.obj <- DESeq(deseq2.obj)  

#get results
deseq2.obj.res <-results(deseq2.obj)
deseq2.obj.res.df <- as.data.frame(deseq2.obj.res)
write.csv(deseq2.obj.res.df,file="30day_XXXYD_deseq2_res.csv")

#MA plot 
DESeq2::plotMA(deseq2.obj.res,alpha=0.001)

#总的差异基因list
colnames(deseq2.obj.res)
select.log2FC <- (deseq2.obj.res$log2FoldChange > 1 | deseq2.obj.res$log2FoldChange < -1)
select.log2FC
table(select.log2FC)

select.adjp <- (deseq2.obj.res$padj <0.05)
table(select.adjp)

select.vec <- (select.log2FC & select.adjp)
table(select.vec)

degs_list <- rownames(deseq2.obj.res)[select.vec]
degs_list <- as.character(degs_list)
write.csv(degs_list,file="30day_XXXYD_deseq2_list.csv")

#---------------------------------------------------------------------------------------------------------------------
#XYDXY

#filter
count_df_XYDXY.filter <- count_df_XYDXY[rowSums(count_df_XYDXY) > 1,]


#condition table
sample_df_XYDXY <- data.frame(
  condition = c( rep("XYD",3), rep("XY",3)),
  cell_line = "Ola"
)
rownames(sample_df_XYDXY) <- colnames(count_df_XYDXY.filter)

library(DESeq2)

#directly get test result
deseq2.obj <-DESeqDataSetFromMatrix(countData = count_df_XYDXY.filter,colData = sample_df_XYDXY,design = ~condition)
deseq2.obj


deseq2.obj <-estimateSizeFactors(deseq2.obj)
sizeFactors(deseq2.obj)

#test(RLE)
deseq2.obj <- DESeq(deseq2.obj)  

#get results
deseq2.obj.res <-results(deseq2.obj)
deseq2.obj.res.df <- as.data.frame(deseq2.obj.res)
write.csv(deseq2.obj.res.df,file="30day_XYDXY_deseq2_res.csv")

#MA plot 
DESeq2::plotMA(deseq2.obj.res,alpha=0.001)

#总的差异基因list
colnames(deseq2.obj.res)
select.log2FC <- (deseq2.obj.res$log2FoldChange > 1 | deseq2.obj.res$log2FoldChange < -1)
select.log2FC
table(select.log2FC)

select.adjp <- (deseq2.obj.res$padj <0.05)
table(select.adjp)

select.vec <- (select.log2FC & select.adjp)
table(select.vec)

degs_list <- rownames(deseq2.obj.res)[select.vec]
degs_list <- as.character(degs_list)
write.csv(degs_list,file="30day_XYDXY_deseq2_list.csv")


#---------------------------------------------------------------------------------------------------------------------
#XXXY

#filter
count_df_XXXY.filter <- count_df_XXXY[rowSums(count_df_XXXY) > 1,]


#condition table
sample_df_XXXY <- data.frame(
  condition = c( rep("XX",3), rep("XY",3)),
  cell_line = "Ola"
)
rownames(sample_df_XXXY) <- colnames(count_df_XXXY.filter)

library(DESeq2)

#directly get test result
deseq2.obj <-DESeqDataSetFromMatrix(countData = count_df_XXXY.filter,colData = sample_df_XXXY,design = ~condition)
deseq2.obj


deseq2.obj <-estimateSizeFactors(deseq2.obj)
sizeFactors(deseq2.obj)

#test(RLE)
deseq2.obj <- DESeq(deseq2.obj)  

#get results
deseq2.obj.res <-results(deseq2.obj)
deseq2.obj.res.df <- as.data.frame(deseq2.obj.res)
write.csv(deseq2.obj.res.df,file="30day_XXXY_deseq2_res.csv")

#MA plot 
DESeq2::plotMA(deseq2.obj.res,alpha=0.001)

#总的差异基因list
colnames(deseq2.obj.res)
select.log2FC <- (deseq2.obj.res$log2FoldChange > 1 | deseq2.obj.res$log2FoldChange < -1)
select.log2FC
table(select.log2FC)

select.adjp <- (deseq2.obj.res$padj <0.05)
table(select.adjp)

select.vec <- (select.log2FC & select.adjp)
table(select.vec)

degs_list <- rownames(deseq2.obj.res)[select.vec]
degs_list <- as.character(degs_list)
write.csv(degs_list,file="30day_XXXY_deseq2_list.csv")
