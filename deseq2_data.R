args <- commandArgs(T)
ctrl_1 <- args[1]
ctrl_2 <- args[2]
treat_1 <- args[3]
treat_2 <- args[4]
out_dir <- args[5]
fold_change <- as.numeric(args[6])
padj_value <- as.numeric(args[7])

# setwd("D:/MyDoc/OneDrive/Doc/Lab/analysis/yxh_ce")
# sample <- "ce1_ce2_ce3_ce4"
# ctrl_1 <- "result_treat/ce1/ce1_featureCounts.txt"
# ctrl_2 <- "result_treat/ce2/ce2_featureCounts.txt"
# treat_1 <- "result_treat/ce3/ce3_featureCounts.txt"
# treat_2 <- "result_treat/ce4/ce4_featureCounts.txt"
# out_dir <- "degenes/ce1_ce2_ce3_ce4"
# fold_change <- 0.585
# padj_value <- 0.05

library("DESeq2")
data <- c(ctrl_1, ctrl_2, treat_1, treat_2)

get_head_name <- function(name)
{
  head_name <- strsplit(name, '/')
  xlen <- length(head_name[[1]])
  head_name <- head_name[[1]][xlen - 1]
  return(head_name)
}

#read all featurecounts results
allgene <- read.table(data[1], header=F)[-1,c(1,7)]
head_name <- get_head_name(data[1])
names(allgene) <- c("gene", head_name)
data2 <- tail(data, -1)

for (name in data2)
{
  fc_table <- read.table(name, header=F)[-1,c(1,7)]
  head_name <- get_head_name(name)
  names(fc_table) <- c("gene", head_name)
  allgene <- merge(allgene, fc_table, by="gene", all=T)
}

#treat the feature count data
us_count <- allgene[,-1]
us_count <- as.data.frame(lapply(us_count, as.numeric))
us_count <- round(us_count, digits=0) #get interger
rownames(us_count) <- allgene[,1]

#prepare the data
us_count <- as.matrix(us_count)
condition <- factor(c("ctrl", "ctrl", "treat", "treat"),
                    levels=c("ctrl", "treat"))
coldata <- data.frame(row.names=colnames(us_count), condition)

#deseq2 analyze
dds <- DESeqDataSetFromMatrix(us_count, coldata, design=~condition)
dds <- DESeq(dds)             #standardization dds
res <- results(dds)           #get result

#get up/down/all gene from result, set the cut off
res <- res[order(res$padj),]
resdata <-merge(as.data.frame(res),
                as.data.frame(counts(dds, normalize=TRUE)),
                by="row.names",sort=FALSE)
deseq_res <- data.frame(resdata)
up_diff <- subset(deseq_res,(padj < padj_value) & (log2FoldChange > fold_change))
down_diff <- subset(deseq_res,(padj < padj_value) & (log2FoldChange < -fold_change))

# write result list
if (!file.exists(out_dir))
  {dir.create(out_dir, recursive=TRUE)}
write.table(up_diff$Row.names, 
            paste(out_dir, "/upgenelist.txt",sep=""),
            quote=F, row.names=F, col.names=F, sep="\t")
write.table(down_diff$Row.names, 
            paste(out_dir, "/downgenelist.txt",sep=""),
            quote=F, row.names=F, col.names=F, sep="\t")

