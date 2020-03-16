##load packages

library(edgeR)
library(ggplot2)
require(grDevices)

##Load data

setwd('/Users/sebastianhutchinson/Documents/RNAseq/Lucy_VSG_seq/F18FTSEUHT1406_TRYzdpR/Clean/Publication_analysis/')
raw_data <- read.table(file = '',
                       header = TRUE,
                       row.names = 'Geneid')

genomeBED <- as.data.frame(read.table('All_genes.bed', 
                                      header = FALSE, 
                                      sep = "\t",
                                      as.is = "V1"))

chromosomes <- as.data.frame(read.table('Chromosome_lengths.txt', 
                                        header = TRUE, 
                                        as.is = "Contig"))

VSG_ES <- as.data.frame(read.table('BES_chlengs.txt', 
                                   header = TRUE, 
                                   as.is = "Contig"))

VSG_ES_genes <- as.data.frame(read.table('VSG_ES_Genes.bed', 
                                         header = FALSE, 
                                         as.is = "V1", 
                                         sep = "\t"))
VSGs_BED <- read.table('VariantSurfaceGlycoprotein.bed', 
                       header = FALSE, 
                       row.names = 'V4', 
                       sep = '\t')

# combine data sets into a matrix
# edit the names of samples for your data
geneCounts <- raw_data[ , 6:13]
samples <- c('X70aU', 'X70bU', 
             'X70ai', 'X70bi', 
             'RAD50.3U', 'RAD50.4U', 
             'RAD50.3i', 'RAD50.4i')
head(geneCounts)

#Perform differential expression 
group <- factor(c(1,1,2,2,3,3,4,4))
y <- DGEList(counts=geneCounts, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
Cond1_UvsI <- exactTest(y, pair = c(1,2))
Cond2_UvsI <- exactTest(y, pair = c(3,4))

#Merge tables and generage significant lists -Cond1
Cond1_BED <- as.data.frame(merge(VSGs_BED,Cond1_UvsI$table, all.x = TRUE, by = 0, sort = TRUE))
Cond1_bool <- Cond1_BED$logFC > 2 & Cond1_BED$PValue < 0.05
Cond1_BED_sig <- Cond1_BED[Cond1_bool,]
#Save tables and reinitialise with as.is function to avoid factors (important for the genome plot)
write.table(x = Cond1_BED_sig, file = "Cond1_BED_sig.bed", quote = FALSE, sep = "\t", row.names = FALSE)
Cond1_BED_sig <- as.data.frame(read.table('Cond1_BED_sig.bed', header = TRUE, sep = "\t",as.is = "V1"))
  
#Merge tables and generage significant lists -Cond2
Cond2_BED <- as.data.frame(merge(VSGs_BED, Cond2_UvsI$table, all.x = TRUE, by = 0, sort = TRUE))
Cond2_bool <- Cond2_BED$logFC > 2 & Cond2_BED$PValue < 0.05
Cond2_BED_sig <- Cond2_BED[Cond2_bool,]
#Save tables and reinitialise with as.is function to avoid factors (important for the genome plot)
write.table(x = Cond2_BED_sig, file = "Cond2_BED_sig.bed", quote = FALSE, sep = "\t", row.names = FALSE)
Cond2_BED_sig <- as.data.frame(read.table('Cond2_BED_sig.bed', header = TRUE, sep = "\t",as.is = "V1"))

merged_Cond1_Cond2 <- merge(x = Cond2_BED_sig, y= Cond1_BED_sig, all = TRUE, by = "Row.names", sort = TRUE)
write.table(x = merged_Cond1_Cond2, file = "MergedCond1_Cond2.txt", quote = FALSE, sep = "\t")

png("MDS_Plot.png", width = 500,height = 500, units = "px")
plotMDS(y)


#Draw Volcano plots
Cond1_volcano <- ggplot(Cond1_BED , 
       aes(x = logFC, 
           y = -log10(PValue)))+ 
  xlim(c(-10,15)) +
  ylim(c(0,25)) + 
  geom_point(shape = 21, 
             alpha = 0.5, 
             size = 2, 
             aes(fill = Cond1_bool, colour = Cond1_bool)) +
  scale_fill_manual(values = c("black", "red")) +
  scale_colour_manual(values = c("black", "red"))
ggsave(filename = "Cond1_volcano.png", plot = Cond1_volcano, device = "png", width = 10, height = 15, dpi = "print")  

Cond2_volcano <- ggplot(Cond2_BED , 
       aes(x = logFC, 
           y = -log10(PValue))) + 
  xlim(c(-10,15)) +
  ylim(c(0,25)) +
  geom_point(shape = 21, 
             alpha = 0.5, 
             size = 2, 
             aes(fill = Cond2_bool, colour = Cond2_bool)) +
  scale_fill_manual(values = c("black", "red")) +
  scale_colour_manual(values = c("black", "red"))
ggsave(filename = "Cond2_volcano.png", plot = Cond2_volcano, device = "png", width = 10, height = 15, dpi = "print")

#Draw genome plot

pdf("genome_plot.pdf", paper = "a4")
## set up the plot region:
op <- par(bg = "White", bty = 'n')
plot(x= c(0,7000000), y=c(0,120000),type = "n", xlab = "", ylab = "", axes = FALSE)

# Draw on the gene models from Nicolai's paper
for(i in 1:nrow(genomeBED)) {
  for (j in 1:nrow(chromosomes)) {
    if (genomeBED$V1[i] == chromosomes$Contig[j]){
      rect(xleft = (chromosomes$xStart[j] + genomeBED$V2[i]),
           xright = (chromosomes$xStart[j] +  genomeBED$V3[i]),
           ybottom = chromosomes$ystart[j],
           ytop = chromosomes$yend[j] + if(all(genomeBED$V6[i] == "+")){750}else{-750},
           border = 'gray')
    }
  }
}

#draw on the Cond1 significant diff genes
for(i in 1:nrow(Cond1_BED_sig)) {
  for (j in 1:nrow(chromosomes)) {
    if (Cond1_BED_sig$V1[i] == chromosomes$Contig[j]){
      if(!(grepl(pattern = "BES",Cond1_BED_sig$V1[i]))){
      rect(xleft = (chromosomes$xStart[j] + Cond1_BED_sig$V2[i]),
           xright = (chromosomes$xStart[j] +  Cond1_BED_sig$V3[i]),
           ybottom = chromosomes$ystart[j] + 1000,
           ytop = chromosomes$yend[j] + 1750,
           border = 'blue')
      }
    }
  }
}
#draw on the Cond2 significant diff genes
for(i in 1:nrow(Cond2_BED_sig)) {
  for (j in 1:nrow(chromosomes)) {
    if (Cond2_BED_sig$V1[i] == chromosomes$Contig[j]){
      if(!(grepl(pattern = "BES",Cond2_BED_sig$V1[i]))){
      rect(xleft = (chromosomes$xStart[j] + Cond2_BED_sig$V2[i]),
           xright = (chromosomes$xStart[j] +  Cond2_BED_sig$V3[i]),
           ybottom = chromosomes$ystart[j] + 2000,
           ytop = chromosomes$yend[j] + 2750,
           border = 'red')
      }
    }
  }
}

##Draw chromosomes
for(i in 1:nrow(chromosomes)) {
  rect(xleft = chromosomes$xStart[i],
       xright = chromosomes$xStart[i] + chromosomes$Length[i],
       ybottom = chromosomes$ystart[i],
       ytop = chromosomes$yend[i],
       border = 'black')
}


##Draw VSG-ES genes
for(i in 1:nrow(VSG_ES_genes)) {
  for (j in 1:nrow(VSG_ES)) {
        if (VSG_ES_genes$V1[i] == VSG_ES$Contig[j]){
      rect(xleft = (VSG_ES$xStart[j] + (VSG_ES_genes$V2[i]*10)),
           xright = (VSG_ES$xStart[j]  +  (VSG_ES_genes$V3[i]*10)),
           ybottom = VSG_ES$ystart[j],
           ytop = VSG_ES$yend[j] + 1000,
           border = 'gray',
           col = 'gray')
    }
  }
}


##Draw VSG-ES contigs
for(i in 1:nrow(VSG_ES)) {
  rect(xleft = VSG_ES$xStart[i],
       xright = VSG_ES$xStart[i] + (VSG_ES$Length[i] * 10),
       ybottom = VSG_ES$ystart[i],
       ytop = VSG_ES$yend[i],
       border = 'black')
}

#draw on the Cond1 significant diff genes for BES
for(i in 1:nrow(Cond1_BED_sig)) {
  for (j in 1:nrow(VSG_ES)) {
    if (Cond1_BED_sig$V1[i] == VSG_ES$Contig[j]){
      rect(xleft = (VSG_ES$xStart[j] + (Cond1_BED_sig$V2[i]*10)),
           xright = (VSG_ES$xStart[j] +  (Cond1_BED_sig$V3[i]*10)),
           ybottom = VSG_ES$ystart[j] + 1250,
           ytop = VSG_ES$yend[j] + 2000,
           border = 'blue',
           col = 'blue')
    }
  }
}
#draw on the Cond2 significant diff genes for BES
for(i in 1:nrow(Cond2_BED_sig)) {
  for (j in 1:nrow(VSG_ES)) {
    if (Cond2_BED_sig$V1[i] == VSG_ES$Contig[j]){
      rect(xleft = (VSG_ES$xStart[j] + (Cond2_BED_sig$V2[i]*10)),
           xright = (VSG_ES$xStart[j] +  (Cond2_BED_sig$V3[i]*10)),
           ybottom = VSG_ES$ystart[j] + 2250,
           ytop = VSG_ES$yend[j] + 3000,
           border = 'red',
           col = "red")
    }
  }
}
dev.off() 