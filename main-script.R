## Setup
### Bioconductor and CRAN libraries used
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ggrepel)
library(cowplot)
library(BiocManager)
library(devtools)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationDbi)
library(org.At.tair.db)
library(enrichplot)
library(SPIA)
library(WGCNA)

## Load in data
dataRaw <- read.table("data/countData.txt", header=T, row.names=1) 
meta <- read.table("data/colData.txt", header=T, row.names=1)
data <- round(dataRaw)
entrez <- read.delim("data/entrez.txt", header=T)
entrez <- entrez[,c(1,2,4)]

tx2gene <- read_csv("data/tx2gene.csv")
annot <- tx2gene %>% 
  dplyr::select(ensgene, symbol) %>% 
  dplyr::distinct()

##create data sets
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Type)
ddsB <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Type)
ddsC <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Type)

ddsB$Type <- relevel(ddsB$Type, ref="B")
ddsC$Type <- relevel(ddsC$Type, ref="C")

##PCA and heatmap 
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="Type")
rld_mat <- assay(rld)   
rld_cor <- cor(rld_mat)    ## cor() is a base R function
pheatmap(rld_cor, annotation=meta)


##create all pairwise results tables (lfc = shrunken log fold changes)
dds <- DESeq(dds)
ddsB <- DESeq(ddsB)
ddsC <- DESeq(ddsC)

res_tableBA <- results(dds, contrast=c("Type", "B", "A"), alpha = 0.05, lfcThreshold = 0.41)
res_tableCA <- results(dds, contrast=c("Type", "C", "A"), alpha = 0.05, lfcThreshold = 0.41)
res_tableDA <- results(dds, contrast=c("Type", "D", "A"), alpha = 0.05, lfcThreshold = 0.41)
res_tableBAlfc <- lfcShrink(dds, coef="Type_B_vs_A", res=res_tableBA)
res_tableCAlfc <- lfcShrink(dds, coef="Type_C_vs_A", res=res_tableCA)
res_tableDAlfc <- lfcShrink(dds, coef="Type_D_vs_A", res=res_tableDA)

res_tableCAstrict <- results(dds, contrast=c("Type", "C", "A"), alpha = 0.05, lfcThreshold = 1)
res_tableCAlfcstrict <- lfcShrink(dds, coef="Type_C_vs_A", res=res_tableCAstrict)
res_tableCA_tbstrict <- res_tableCAlfcstrict %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
res_tableCA_tbstrict <- res_tableCA_tbstrict[,c(1:7,9,8)]
sigCAstrict <- res_tableCA_tbstrict %>%
  filter(padj < padj.cutoff)

res_tableBAstrict <- results(dds, contrast=c("Type", "B", "A"), alpha = 0.05, lfcThreshold = 1)
res_tableBAlfcstrict <- lfcShrink(dds, coef="Type_B_vs_A", res=res_tableBAstrict)
res_tableBA_tbstrict <- res_tableBAlfcstrict %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
res_tableBA_tbstrict <- res_tableBA_tbstrict[,c(1:7,9,8)]
sigBAstrict <- res_tableBA_tbstrict %>%
  filter(padj < padj.cutoff)

res_tableDBstrict <- results(ddsB, contrast=c("Type", "D", "B"), alpha = 0.05, lfcThreshold = 1)
res_tableDBlfcstrict <- lfcShrink(ddsB, coef="Type_D_vs_B", res=res_tableDBstrict)
res_tableDB_tbstrict <- res_tableDBlfcstrict %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
res_tableDB_tbstrict <- res_tableDB_tbstrict[,c(1:7,9,8)]
sigDBstrict <- res_tableDB_tbstrict %>%
  filter(padj < padj.cutoff)

res_tableDCstrict <- results(ddsC, contrast=c("Type", "D", "C"), alpha = 0.05, lfcThreshold = 1)
res_tableDClfcstrict <- lfcShrink(ddsC, coef="Type_D_vs_C", res=res_tableDCstrict)
res_tableDC_tbstrict <- res_tableDClfcstrict %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
res_tableDC_tbstrict <- res_tableDC_tbstrict[,c(1:7,9,8)]
sigDCstrict <- res_tableDC_tbstrict %>%
  filter(padj < padj.cutoff)

sigCAstrict_mapman <- sigCAstrict %>% arrange(desc(log2FoldChange))
sigCAstrict_mapman <- sigCAstrict_mapman[,c(1,3)]
write_delim(sigCAstrict_mapman, "mapman/sigCA_lfcthreshhold_1.txt")
sigCA_mapman <- sigCA %>% arrange(desc(log2FoldChange))
sigCA_mapman <- sigCA_mapman[,c(1,3)]
write_delim(sigCA_mapman, "mapman/sigCA_lfcthreshhold_0.41.txt")

sigBAstrict_mapman <- sigBAstrict %>% arrange(desc(log2FoldChange))
sigBAstrict_mapman <- sigBAstrict_mapman[,c(1,3)]
write_delim(sigBAstrict_mapman, "mapman/sigBA_lfcthreshhold_1.txt")
sigBA_mapman <- sigBA %>% arrange(desc(log2FoldChange))
sigBA_mapman <- sigBA_mapman[,c(1,3)]
write_delim(sigBA_mapman, "mapman/sigBA_lfcthreshhold_0.41.txt")

sigDBstrict_mapman <- sigDBstrict %>% arrange(desc(log2FoldChange))
sigDBstrict_mapman <- sigDBstrict_mapman[,c(1,3)]
write_delim(sigDBstrict_mapman, "mapman/sigDB_lfcthreshhold_1.txt")
sigDB_mapman <- sigDB %>% arrange(desc(log2FoldChange))
sigDB_mapman <- sigDB_mapman[,c(1,3)]
write_delim(sigDB_mapman, "mapman/sigDB_lfcthreshhold_0.41.txt")

sigDCstrict_mapman <- sigDCstrict %>% arrange(desc(log2FoldChange))
sigDCstrict_mapman <- sigDCstrict_mapman[,c(1,3)]
write_delim(sigDCstrict_mapman, "mapman/sigDC_lfcthreshhold_1.txt")
sigDC_mapman <- sigDC %>% arrange(desc(log2FoldChange))
sigDC_mapman <- sigDC_mapman[,c(1,3)]
write_delim(sigDC_mapman, "mapman/sigDC_lfcthreshhold_0.41.txt")

res_tableCB <- results(ddsB, contrast=c("Type", "C", "B"), alpha = 0.05, lfcThreshold = 0.41)
res_tableCBlfc <- lfcShrink(ddsB, coef="Type_C_vs_B", res=res_tableCB)
res_tableDB <- results(ddsB, contrast=c("Type", "D", "B"), alpha = 0.05, lfcThreshold = 0.41)
res_tableDBlfc <- lfcShrink(ddsB, coef="Type_D_vs_B", res=res_tableDB)

res_tableDC <- results(ddsC, contrast=c("Type", "D", "C"), alpha = 0.05, lfcThreshold = 0.41)
res_tableDClfc <- lfcShrink(ddsC, coef="Type_D_vs_C", res=res_tableDC)

#make easier to read data tables of the results
res_tableBA_tb <- res_tableBAlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
res_tableCA_tb <- res_tableCAlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
res_tableDA_tb <- res_tableDAlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))

res_tableCB_tb <- res_tableCBlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
res_tableDB_tb <- res_tableDBlfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))

res_tableDC_tb <- res_tableDClfc %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))

res_tableBA_tb <- res_tableBA_tb[,c(1:7,9,8)]
res_tableCA_tb <- res_tableCA_tb[,c(1:7,9,8)]
res_tableDA_tb <- res_tableDA_tb[,c(1:7,9,8)]

res_tableCB_tb <- res_tableCB_tb[,c(1:7,9,8)]
res_tableDB_tb <- res_tableDB_tb[,c(1:7,9,8)]

res_tableDC_tb <- res_tableDC_tb[,c(1:7,9,8)]


#p value cutoff for significant genes
padj.cutoff <- 0.05

sigBA <- res_tableBA_tb %>%
  filter(padj < padj.cutoff)
sigCA <- res_tableCA_tb %>%
  filter(padj < padj.cutoff)
sigDA <- res_tableDA_tb %>%
  filter(padj < padj.cutoff)
sigCB <- res_tableCB_tb %>%
  filter(padj < padj.cutoff)
sigDB <- res_tableDB_tb %>%
  filter(padj < padj.cutoff)
sigDC <- res_tableDC_tb %>%
  filter(padj < padj.cutoff)

write.csv(sigBA, "results/BA/sigBA.csv")
write.csv(sigCA, "results/CA/sigCA.csv")
write.csv(sigDA, "results/DA/sigDA.csv")
write.csv(sigCB, "results/CB/sigCB.csv")
write.csv(sigDB, "results/DB/sigDB.csv")
write.csv(sigDC, "results/DC/sigDC.csv")

#create results tables with entrez ids, because KEGG mapper needs entrez ids

resBA_entrez <- res_tableBA_tb[which(duplicated(res_tableBA_tb$Entrez) == F), ]
resBA_entrez <- dplyr::filter(resBA_entrez, Entrez != "NA")
resCA_entrez <- res_tableCA_tb[which(duplicated(res_tableCA_tb$Entrez) == F), ]
resCA_entrez <- dplyr::filter(resCA_entrez, Entrez != "NA")
resDA_entrez <- res_tableDA_tb[which(duplicated(res_tableDA_tb$Entrez) == F), ]
resDA_entrez <- dplyr::filter(resDA_entrez, Entrez != "NA")

resCB_entrez <- res_tableCB_tb[which(duplicated(res_tableCB_tb$Entrez) == F), ]
resCB_entrez <- dplyr::filter(resCB_entrez, Entrez != "NA")
resDB_entrez <- res_tableDB_tb[which(duplicated(res_tableDB_tb$Entrez) == F), ]
resDB_entrez <- dplyr::filter(resDB_entrez, Entrez != "NA")

resDC_entrez <- res_tableDC_tb[which(duplicated(res_tableDC_tb$Entrez) == F), ]
resDC_entrez <- dplyr::filter(resDC_entrez, Entrez != "NA")

BAentrezlfc <- resBA_entrez$log2FoldChange
names(BAentrezlfc) <- resBA_entrez$Entrez
BAentrezlfc <- sort(BAentrezlfc, decreasing = TRUE)
BAlfcordered <- res_tableBA_tb$log2FoldChange
names(BAlfcordered) <- res_tableBA_tb$gene
BAlfcordered <- sort(BAlfcordered, decreasing=TRUE)

CAentrezlfc <- resCA_entrez$log2FoldChange
names(CAentrezlfc) <- resCA_entrez$Entrez
CAentrezlfc <- sort(CAentrezlfc, decreasing = TRUE)
CAlfcordered <- res_tableCA_tb$log2FoldChange
names(CAlfcordered) <- res_tableCA_tb$gene
CAlfcordered <- sort(CAlfcordered, decreasing=TRUE)

DAentrezlfc <- resDA_entrez$log2FoldChange
names(DAentrezlfc) <- resDA_entrez$Entrez
DAentrezlfc <- sort(DAentrezlfc, decreasing = TRUE)
DAlfcordered <- res_tableDA_tb$log2FoldChange
names(DAlfcordered) <- res_tableDA_tb$gene
DAlfcordered <- sort(DAlfcordered, decreasing=TRUE)

CBentrezlfc <- resCB_entrez$log2FoldChange
names(CBentrezlfc) <- resCB_entrez$Entrez
CBentrezlfc <- sort(CBentrezlfc, decreasing = TRUE)
CBlfcordered <- res_tableCB_tb$log2FoldChange
names(CBlfcordered) <- res_tableCB_tb$gene
CBlfcordered <- sort(CBlfcordered, decreasing=TRUE)

DBentrezlfc <- resDB_entrez$log2FoldChange
names(DBentrezlfc) <- resDB_entrez$Entrez
DBentrezlfc <- sort(DBentrezlfc, decreasing = TRUE)
DBlfcordered <- res_tableDB_tb$log2FoldChange
names(DBlfcordered) <- res_tableDB_tb$gene
DBlfcordered <- sort(DBlfcordered, decreasing=TRUE)

DCentrezlfc <- resDC_entrez$log2FoldChange
names(DCentrezlfc) <- resDC_entrez$Entrez
DCentrezlfc <- sort(DCentrezlfc, decreasing = TRUE)
DClfcordered <- res_tableDC_tb$log2FoldChange
names(DClfcordered) <- res_tableDC_tb$gene
DClfcordered <- sort(DClfcordered, decreasing=TRUE)


## gene set enrichment analysis using gene sets from KEGG pathways
BAgseKEGG <- gseKEGG(geneList = BAlfcordered,
                    organism = "ath",
                    pvalueCutoff = 0.05, 
                    verbose = FALSE)
BAgseKEGG_results <- BAgseKEGG@result
CAgseKEGG <- gseKEGG(geneList = CAlfcordered,
                     organism = "ath",
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)
CAgseKEGG_results <- CAgseKEGG@result
DAgseKEGG <- gseKEGG(geneList = DAlfcordered,
                     organism = "ath",
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)
DAgseKEGG_results <- DAgseKEGG@result

CBgseKEGG <- gseKEGG(geneList = CBlfcordered,
                     organism = "ath",
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)
CBgseKEGG_results <- CBgseKEGG@result
DBgseKEGG <- gseKEGG(geneList = DBlfcordered,
                     organism = "ath",
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)
DBgseKEGG_results <- DBgseKEGG@result

DCgseKEGG <- gseKEGG(geneList = DClfcordered,
                     organism = "ath",
                     pvalueCutoff = 0.05, 
                     verbose = FALSE)
DCgseKEGG_results <- DCgseKEGG@result

write.csv(BAgseKEGG_results, "results/BA/BA_gseKEGG_results.csv")
write.csv(CAgseKEGG_results, "results/CA/CA_gseKEGG_results.csv")
write.csv(DAgseKEGG_results, "results/DA/DA_gseKEGG_results.csv")

write.csv(CBgseKEGG_results, "results/CB/CB_gseKEGG_results.csv")
write.csv(DBgseKEGG_results, "results/DB/DB_gseKEGG_results.csv")

write.csv(DCgseKEGG_results, "results/DC/DC_gseKEGG_results.csv")

detach("package:dplyr", unload=TRUE)
## Output images for all significant KEGG pathways
get_kegg_plotsBA <- function(x) {
  pathview(gene.data = BAentrezlfc, 
           out.suffix = "BA",
           pathway.id = BAgseKEGG_results$ID[x],
           species = "ath",
           low = list(gene = "red", cpd = "yellow"),
           high = list(gene = "green", cpd = "blue"),
           limit = list(gene = 2, # value gives the max/min limit for foldchanges
                        cpd = 1))
}
purrr::map(1:length(BAgseKEGG_results$ID), get_kegg_plotsBA)
get_kegg_plotsCA <- function(x) {
  pathview(gene.data = CAentrezlfc, 
           out.suffix = "CA",
           pathway.id = CAgseKEGG_results$ID[x],
           species = "ath",
           low = list(gene = "red", cpd = "yellow"),
           high = list(gene = "green", cpd = "blue"),
           limit = list(gene = 2, # value gives the max/min limit for foldchanges
                        cpd = 1))
}
purrr::map(1:length(CAgseKEGG_results$ID), get_kegg_plotsCA)
get_kegg_plotsDA <- function(x) {
  pathview(gene.data = DAentrezlfc, 
           out.suffix = "DA",
           pathway.id = DAgseKEGG_results$ID[x],
           species = "ath",
           low = list(gene = "red", cpd = "yellow"),
           high = list(gene = "green", cpd = "blue"),
           limit = list(gene = 2, # value gives the max/min limit for foldchanges
                        cpd = 1))
}
purrr::map(1:length(DAgseKEGG_results$ID), get_kegg_plotsDA)

get_kegg_plotsCB <- function(x) {
  pathview(gene.data = CBentrezlfc, 
           out.suffix = "CB",
           pathway.id = CBgseKEGG_results$ID[x],
           species = "ath",
           low = list(gene = "red", cpd = "yellow"),
           high = list(gene = "green", cpd = "blue"),
           limit = list(gene = 2, # value gives the max/min limit for foldchanges
                        cpd = 1))
}
purrr::map(1:length(CBgseKEGG_results$ID), get_kegg_plotsCB)
get_kegg_plotsDB <- function(x) {
  pathview(gene.data = DBentrezlfc, 
           out.suffix = "DB",
           pathway.id = DBgseKEGG_results$ID[x],
           species = "ath",
           low = list(gene = "red", cpd = "yellow"),
           high = list(gene = "green", cpd = "blue"),
           limit = list(gene = 2, # value gives the max/min limit for foldchanges
                        cpd = 1))
}
purrr::map(1:length(DBgseKEGG_results$ID), get_kegg_plotsDB)

get_kegg_plotsDC <- function(x) {
  pathview(gene.data = DCentrezlfc, 
           out.suffix = "DC",
           pathway.id = DCgseKEGG_results$ID[x],
           species = "ath",
           low = list(gene = "red", cpd = "yellow"),
           high = list(gene = "green", cpd = "blue"),
           limit = list(gene = 2, # value gives the max/min limit for foldchanges
                        cpd = 1))
}
purrr::map(1:length(DCgseKEGG_results$ID), get_kegg_plotsDC)

pathview(gene.data = DCentrezlfc,
         out.suffix = "DC",
         pathway.id = "ath03008",
         species = "ath",
         low = list(gene = "red", cpd = "yellow"),
         high = list(gene = "green", cpd = "blue"),
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))

# ## gene set enrichment analysis using gene sets associated with BP Gene Ontology terms
# CAgseaGO <- gseGO(geneList = CAentrezlfc, 
#                 OrgDb = org.At.tair.db, 
#                 ont = 'BP', 
#                 pvalueCutoff = 0.05,
#                 verbose = FALSE) 
# CAgseaGO_results <- CAgseaGO@result

library(dplyr)

bhlh_meta <- meta %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

normalized_counts <- counts(dds, normalized=T) %>% #create tibbles of normalized counts, annotated by gene symbol
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  left_join(annot, by=c("gene" = "ensgene")) %>%
  left_join(entrez, by=c("gene" = "From"))
write.csv(normalized_counts, "normalized_counts.csv")

########## plot counts for one gene
# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="AT2G42300", intgroup="Type", returnData=TRUE)
# What is the data output of plotCounts()?
d %>% View()

# Plot the normalized counts, using the Types (rownames(d) as labels)
ggplot(d, aes(x = Type, y = count, color = Type)) +
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) +
  theme_bw() +
  ggtitle("BHLH48") +
  theme(plot.title = element_text(hjust = 0.5))

######heatmaps (use normalized counts of significant genes)
total_norm <- normalized_counts[,c(1:13)]
total_norm <- total_norm[,c(1:4,8:10,5:7,11:13)]

allSigs <- full_join(sigBAstrict, sigDBstrict)
allSigs <- full_join(allSigs, sigDCstrict)
allSigs <- full_join(allSigs, sigCAstrict)
allSigs_norm <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% allSigs$gene)  
allSigs_norm <- allSigs_norm[,c(1:4,8:10,5:7,11:13)]

sigCA_norm <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% sigCA$gene)  
sigCA_norm <- sigCA_norm[,c(1:4,8:10,5:7,11:13)]

sigBA_norm <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% sigBA$gene)  
sigBA_norm <- sigBA_norm[,c(1:4,8:10,5:7,11:13)]

sigDA_norm <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% sigDA$gene)  
sigDA_norm <- sigDA_norm[,c(1:4,8:10,5:7,11:13)]

sigDB_norm <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% sigDB$gene)  
sigDB_norm <- sigDB_norm[,c(1:4,8:10,5:7,11:13)]

sigCB_norm <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% sigCB$gene)  
sigCB_norm <- sigCB_norm[,c(1:4,8:10,5:7,11:13)]

sigDC_norm <- normalized_counts[,c(1:13)] %>% 
  filter(gene %in% sigDC$gene)  
sigDC_norm <- sigDC_norm[,c(1:4,8:10,5:7,11:13)]

# CAupregulated <- res_tableCA_tb %>% filter(log2FoldChange>0)
# CAonly_norm <- sigCA_norm %>%
#   filter(!gene %in% sigDB_norm$gene) %>%
#   filter(!gene %in% sigBA_norm$gene) %>%
#   filter(gene %in% CAupregulated$gene) %>% 
#   left_join(annot, by=c("gene" = "ensgene"))
# 
# commonGenes = CAonly_norm %>% filter(gene %in% pinkGenes$gene)

# write.csv(pinkGenes, "heatmapExtracted.csv")
# write.csv(CAonly_norm, "manuallyFiltered.csv")
# write.csv(commonGenes, "common.csv")

heat_colors <- rev(brewer.pal(11, "RdBu"))

# use this code to interact with the heatmap (highlight it and click Run at the top)
# change sigCA_norm to sigBA_norm, etc if needed
ht = pheatmap(allSigs_norm[2:13],
              color = heat_colors,
              cluster_rows = T,
              show_rownames = F,
              annotation = meta,
              border_color = NA,
              fontsize = 10,
              scale = "row",
              fontsize_row = 10,
              cluster_cols=FALSE,
              height = 15)
ht_shiny(ht)

new_total_norm <- total_norm[which(rowSums(total_norm[2:13]) > 0), ]

new2_total_norm <- total_norm[apply(total_norm, 1, function(row) all(row > 0.1  )), ]

nonzerostest <- normalized_counts[rowSums(normalized_counts[,2:13]) > 0,] 

vs CLASS SORTING
E-11 sigBA (glu vs min WT)  
phyper(44, 4736, (26171-4736), 96, lower.tail=FALSE)
E-15 sigDC (Glu vs Min bHLH)
phyper(72, 9625, (26171-9625), 96, lower.tail=FALSE)
E-06 sigCA (bhlh vs WT min) 
phyper(42, 6246, (26171-6246), 96, lower.tail=FALSE)
E-37 sigDB (bhlh vs WT glu) 
phyper(69, 3896, (26171-3896), 96, lower.tail=FALSE)

vs SUPPLEMENTAL
E-10 sigBA (glu vs min WT)
phyper(220, 4736, (26171-4736), 823, lower.tail=FALSE)
E-22 sigDC (Glu vs Min bHLH)
phyper(435, 9625, (26171-9625), 823, lower.tail=FALSE)
E-37 sigCA (bhlh vs WT min) 
phyper(358, 6246, (26171-6246), 823, lower.tail=FALSE)
E-24 sigDB (bhlh vs WT glu)
phyper(232, 3896, (26171-3896), 823, lower.tail=FALSE)

vs NAKABAYASHI
sigBA 1
phyper(2061, 4736, (26171-4736), 12264, lower.tail=FALSE)
sigDC 0.73
phyper(4485, 9625, (26171-9625), 12264, lower.tail=FALSE)
sigCA 1
phyper(2588, 6246, (26171-6246), 12264, lower.tail=FALSE)
sigDB 1
phyper(1587, 3896, (26171-3896), 12264, lower.tail=FALSE)

vs GLC-TOR
sigBAstrict E-32
phyper(288, 1602, (26171-1602), 2357, lower.tail=FALSE)
sigDCstrict E-248 
phyper(1078, 4588, (26171-4588), 2357, lower.tail=FALSE)
sigCAstrict E-35
phyper(385, 2304, (26171-2304), 2357, lower.tail=FALSE)
sigDBstrict E-25
phyper(245, 1408, (26171-1408), 2357, lower.tail=FALSE)
sigBA E-64
phyper(748, 4736, (26171-4736), 2357, lower.tail=FALSE)
sigDC 0
phyper(1790, 9625, (26171-9625), 2357, lower.tail=FALSE)
sigCA E-106
phyper(1021, 6246, (26171-6246), 2357, lower.tail=FALSE)
sigDB E-57
phyper(634, 3896, (26171-3896), 2357, lower.tail=FALSE)

vs torin
sigBA E-10
phyper(666, 4736, (26171-4736), 2974, lower.tail=FALSE)
sigDC E-73
phyper(1548, 9625, (26171-9625), 2974, lower.tail=FALSE)
sigCA E-50
phyper(1047, 6246, (26171-6246), 2974, lower.tail=FALSE)
sigDB E-08
phyper(545, 3896, (26171-3896), 2974, lower.tail=FALSE)

vs abi5-1_LFC_1.0
sigBA 
phyper(553, 4736, (26171-4736), 1919, lower.tail=FALSE)
sigDC 
phyper(968, 9625, (26171-9625), 1919, lower.tail=FALSE)
sigCA 
phyper(653, 6246, (26171-6246), 1919, lower.tail=FALSE)
sigDB 
phyper(541, 3896, (26171-3896), 1919, lower.tail=FALSE)

vs abi5-1_LFC_0.5
sigBA 
phyper(1664, 4736, (26171-4736), 6702, lower.tail=FALSE)
sigDC 
phyper(3069, 9625, (26171-9625), 6702, lower.tail=FALSE)
sigCA 
phyper(2055, 6246, (26171-6246), 6702, lower.tail=FALSE)
sigDB 
phyper(1455, 3896, (26171-3896), 6702, lower.tail=FALSE)

vs aba (7 days old)
sigBA 
phyper(2194, 4736, (26171-4736), 8755, lower.tail=FALSE)
sigDC 
phyper(4231, 9625, (26171-9625), 8755, lower.tail=FALSE)
sigCA 
phyper(3111, 6246, (26171-6246), 8755, lower.tail=FALSE)
sigDB 
phyper(1881, 3896, (26171-3896), 8755, lower.tail=FALSE)

vs aba (12 days old)
sigBA 
phyper(786, 4736, (26171-4736), 3016, lower.tail=FALSE)
sigDC 
phyper(1330, 9625, (26171-9625), 3016, lower.tail=FALSE)
sigCA 
phyper(1084, 6246, (26171-6246), 3016, lower.tail=FALSE)
sigDB 
phyper(591, 3896, (26171-3896), 3016, lower.tail=FALSE)

vs dfpm (12 days old)
sigBA 
phyper(342, 4736, (26171-4736), 1376, lower.tail=FALSE)
sigDC 
phyper(623, 9625, (26171-9625), 1376, lower.tail=FALSE)
sigCA 
phyper(643, 6246, (26171-6246), 1376, lower.tail=FALSE)
sigDB 
phyper(242, 3896, (26171-3896), 1376, lower.tail=FALSE)

sigCA vs sigDB
phyper(1831, 6246, (26171-6246), 3896, lower.tail=FALSE)

abi5_full <- read_csv("abi5_full.csv")
abi5targets <- read_csv("abi5targets.csv") %>%
  left_join(resBA_entrez[, c(1,3,6)], by=c("gene" = "gene")) %>%
  left_join(resDC_entrez[, c(1,3,6)], by=c("gene" = "gene")) %>%
  left_join(resCA_entrez[, c(1,3,6)], by=c("gene" = "gene")) %>%
  left_join(resDB_entrez[, c(1,3,6)], by=c("gene" = "gene")) %>%
  left_join(abaVSwt[, c(1,10)], by=c("gene" = "gene")) %>% 
  left_join(dfpm[,c(9,6,3)], by=c("gene" = "To")) %>%
  left_join(aba_gse28800[,c(9,6,3)], by=c("gene" = "To")) %>% 
  rename(
    LFC_BA = log2FoldChange.x,
    LFC_DC = log2FoldChange.y,
    LFC_CA = log2FoldChange.x.x,
    LFC_DB = log2FoldChange.y.y,
    LFC_ABA_7d = log2FoldChange.x.x.x,
    LFC_DFPM = logFC.x,
    LFC_ABA_12d = logFC.y,
    padj_DFPM = P.Value.x,
    padj_ABA_12d = P.Value.y,
    padj_BA = padj.x,
    padj_DC = padj.y,
    padj_CA = padj.x.x,
    padj_DB = padj.y.y
  )
abi5targets <- abi5targets[,c(1,2,3,12,16,13,14,8,10,4,6,17,15,9,11,5,7)]
write.csv(abi5targets, "abi5targets_output.csv")

aba_gse28800 <- read_tsv("aba_gse28800.tsv") %>%
  left_join(affymetrix_to_agi, by=c("ID" = "From"))
aba_gse28800coex <- aba_gse28800[aba_gse28800$adj.P.Val<0.05,c(9,6)]
sigCAaba_gse28800 <- sigCA[, c(1,3)] %>%
  left_join((aba_gse28800coex), by=c("gene" = "To"))
sigCAaba_gse28800 <- sigCAaba_gse28800 %>% 
  rename(
    LFC_sigCA = log2FoldChange,
    LFC_sigaba_gse28800 = logFC
  )
sigCAaba_gse28800 <- sigCAaba_gse28800[is.finite(sigCAaba_gse28800$LFC_sigaba_gse28800),]
sigCAaba_gse28800$sign <- sigCAaba_gse28800$LFC_sigCA * sigCAaba_gse28800$LFC_sigaba_gse28800
sigCAaba_gse28800$order[sigCAaba_gse28800$sign>0&sigCAaba_gse28800$LFC_sigaba_gse28800>0] <- 1
sigCAaba_gse28800$order[sigCAaba_gse28800$sign>0&sigCAaba_gse28800$LFC_sigaba_gse28800<0] <- 2
sigCAaba_gse28800$order[sigCAaba_gse28800$sign<0&sigCAaba_gse28800$LFC_sigaba_gse28800>0] <- 3
sigCAaba_gse28800$order[sigCAaba_gse28800$sign<0&sigCAaba_gse28800$LFC_sigaba_gse28800<0] <- 4
sigCAaba_gse28800 <- sigCAaba_gse28800 %>% arrange(order, desc(sign))
write.csv(sigCAaba_gse28800, "sigCAaba_gse28800_sheet.csv")
sigBAaba_gse28800 <- sigBA[, c(1,3)] %>%
  left_join((aba_gse28800coex), by=c("gene" = "To"))
sigBAaba_gse28800 <- sigBAaba_gse28800 %>% 
  rename(
    LFC_sigBA = log2FoldChange,
    LFC_sigaba_gse28800 = logFC
  )
sigBAaba_gse28800 <- sigBAaba_gse28800[is.finite(sigBAaba_gse28800$LFC_sigaba_gse28800),]
sigBAaba_gse28800$sign <- sigBAaba_gse28800$LFC_sigBA * sigBAaba_gse28800$LFC_sigaba_gse28800
sigBAaba_gse28800$order[sigBAaba_gse28800$sign>0&sigBAaba_gse28800$LFC_sigaba_gse28800>0] <- 1
sigBAaba_gse28800$order[sigBAaba_gse28800$sign>0&sigBAaba_gse28800$LFC_sigaba_gse28800<0] <- 2
sigBAaba_gse28800$order[sigBAaba_gse28800$sign<0&sigBAaba_gse28800$LFC_sigaba_gse28800>0] <- 3
sigBAaba_gse28800$order[sigBAaba_gse28800$sign<0&sigBAaba_gse28800$LFC_sigaba_gse28800<0] <- 4
sigBAaba_gse28800 <- sigBAaba_gse28800 %>% arrange(order, desc(sign))
write.csv(sigBAaba_gse28800, "sigBAaba_gse28800_sheet.csv")
sigDBaba_gse28800 <- sigDB[, c(1,3)] %>%
  left_join((aba_gse28800coex), by=c("gene" = "To"))
sigDBaba_gse28800 <- sigDBaba_gse28800 %>% 
  rename(
    LFC_sigDB = log2FoldChange,
    LFC_sigaba_gse28800 = logFC
  )
sigDBaba_gse28800 <- sigDBaba_gse28800[is.finite(sigDBaba_gse28800$LFC_sigaba_gse28800),]
sigDBaba_gse28800$sign <- sigDBaba_gse28800$LFC_sigDB * sigDBaba_gse28800$LFC_sigaba_gse28800
sigDBaba_gse28800$order[sigDBaba_gse28800$sign>0&sigDBaba_gse28800$LFC_sigaba_gse28800>0] <- 1
sigDBaba_gse28800$order[sigDBaba_gse28800$sign>0&sigDBaba_gse28800$LFC_sigaba_gse28800<0] <- 2
sigDBaba_gse28800$order[sigDBaba_gse28800$sign<0&sigDBaba_gse28800$LFC_sigaba_gse28800>0] <- 3
sigDBaba_gse28800$order[sigDBaba_gse28800$sign<0&sigDBaba_gse28800$LFC_sigaba_gse28800<0] <- 4
sigDBaba_gse28800 <- sigDBaba_gse28800 %>% arrange(order, desc(sign))
write.csv(sigDBaba_gse28800, "sigDBaba_gse28800_sheet.csv")
sigDCaba_gse28800 <- sigDC[, c(1,3)] %>%
  left_join((aba_gse28800coex), by=c("gene" = "To"))
sigDCaba_gse28800 <- sigDCaba_gse28800 %>% 
  rename(
    LFC_sigDC = log2FoldChange,
    LFC_sigaba_gse28800 = logFC
  )
sigDCaba_gse28800 <- sigDCaba_gse28800[is.finite(sigDCaba_gse28800$LFC_sigaba_gse28800),]
sigDCaba_gse28800$sign <- sigDCaba_gse28800$LFC_sigDC * sigDCaba_gse28800$LFC_sigaba_gse28800
sigDCaba_gse28800$order[sigDCaba_gse28800$sign>0&sigDCaba_gse28800$LFC_sigaba_gse28800>0] <- 1
sigDCaba_gse28800$order[sigDCaba_gse28800$sign>0&sigDCaba_gse28800$LFC_sigaba_gse28800<0] <- 2
sigDCaba_gse28800$order[sigDCaba_gse28800$sign<0&sigDCaba_gse28800$LFC_sigaba_gse28800>0] <- 3
sigDCaba_gse28800$order[sigDCaba_gse28800$sign<0&sigDCaba_gse28800$LFC_sigaba_gse28800<0] <- 4
sigDCaba_gse28800 <- sigDCaba_gse28800 %>% arrange(order, desc(sign))
write.csv(sigDCaba_gse28800, "sigDCaba_gse28800_sheet.csv")


dfpm <- read_tsv("dfpm.tsv") %>%
  left_join(affymetrix_to_agi, by=c("ID" = "From"))
dfpmcoex <- dfpm[dfpm$adj.P.Val<0.05,c(9,6)]
sigCAdfpm <- sigCA[, c(1,3)] %>%
  left_join((dfpmcoex), by=c("gene" = "To"))
sigCAdfpm <- sigCAdfpm %>% 
  rename(
    LFC_sigCA = log2FoldChange,
    LFC_sigdfpm = logFC
  )
sigCAdfpm <- sigCAdfpm[is.finite(sigCAdfpm$LFC_sigdfpm),]
sigCAdfpm$sign <- sigCAdfpm$LFC_sigCA * sigCAdfpm$LFC_sigdfpm
sigCAdfpm$order[sigCAdfpm$sign>0&sigCAdfpm$LFC_sigdfpm>0] <- 1
sigCAdfpm$order[sigCAdfpm$sign>0&sigCAdfpm$LFC_sigdfpm<0] <- 2
sigCAdfpm$order[sigCAdfpm$sign<0&sigCAdfpm$LFC_sigdfpm>0] <- 3
sigCAdfpm$order[sigCAdfpm$sign<0&sigCAdfpm$LFC_sigdfpm<0] <- 4
sigCAdfpm <- sigCAdfpm %>% arrange(order, desc(sign))
write.csv(sigCAdfpm, "sigCAdfpm_sheet.csv")
sigBAdfpm <- sigBA[, c(1,3)] %>%
  left_join((dfpmcoex), by=c("gene" = "To"))
sigBAdfpm <- sigBAdfpm %>% 
  rename(
    LFC_sigBA = log2FoldChange,
    LFC_sigdfpm = logFC
  )
sigBAdfpm <- sigBAdfpm[is.finite(sigBAdfpm$LFC_sigdfpm),]
sigBAdfpm$sign <- sigBAdfpm$LFC_sigBA * sigBAdfpm$LFC_sigdfpm
sigBAdfpm$order[sigBAdfpm$sign>0&sigBAdfpm$LFC_sigdfpm>0] <- 1
sigBAdfpm$order[sigBAdfpm$sign>0&sigBAdfpm$LFC_sigdfpm<0] <- 2
sigBAdfpm$order[sigBAdfpm$sign<0&sigBAdfpm$LFC_sigdfpm>0] <- 3
sigBAdfpm$order[sigBAdfpm$sign<0&sigBAdfpm$LFC_sigdfpm<0] <- 4
sigBAdfpm <- sigBAdfpm %>% arrange(order, desc(sign))
write.csv(sigBAdfpm, "sigBAdfpm_sheet.csv")
sigDBdfpm <- sigDB[, c(1,3)] %>%
  left_join((dfpmcoex), by=c("gene" = "To"))
sigDBdfpm <- sigDBdfpm %>% 
  rename(
    LFC_sigDB = log2FoldChange,
    LFC_sigdfpm = logFC
  )
sigDBdfpm <- sigDBdfpm[is.finite(sigDBdfpm$LFC_sigdfpm),]
sigDBdfpm$sign <- sigDBdfpm$LFC_sigDB * sigDBdfpm$LFC_sigdfpm
sigDBdfpm$order[sigDBdfpm$sign>0&sigDBdfpm$LFC_sigdfpm>0] <- 1
sigDBdfpm$order[sigDBdfpm$sign>0&sigDBdfpm$LFC_sigdfpm<0] <- 2
sigDBdfpm$order[sigDBdfpm$sign<0&sigDBdfpm$LFC_sigdfpm>0] <- 3
sigDBdfpm$order[sigDBdfpm$sign<0&sigDBdfpm$LFC_sigdfpm<0] <- 4
sigDBdfpm <- sigDBdfpm %>% arrange(order, desc(sign))
write.csv(sigDBdfpm, "sigDBdfpm_sheet.csv")
sigDCdfpm <- sigDC[, c(1,3)] %>%
  left_join((dfpmcoex), by=c("gene" = "To"))
sigDCdfpm <- sigDCdfpm %>% 
  rename(
    LFC_sigDC = log2FoldChange,
    LFC_sigdfpm = logFC
  )
sigDCdfpm <- sigDCdfpm[is.finite(sigDCdfpm$LFC_sigdfpm),]
sigDCdfpm$sign <- sigDCdfpm$LFC_sigDC * sigDCdfpm$LFC_sigdfpm
sigDCdfpm$order[sigDCdfpm$sign>0&sigDCdfpm$LFC_sigdfpm>0] <- 1
sigDCdfpm$order[sigDCdfpm$sign>0&sigDCdfpm$LFC_sigdfpm<0] <- 2
sigDCdfpm$order[sigDCdfpm$sign<0&sigDCdfpm$LFC_sigdfpm>0] <- 3
sigDCdfpm$order[sigDCdfpm$sign<0&sigDCdfpm$LFC_sigdfpm<0] <- 4
sigDCdfpm <- sigDCdfpm %>% arrange(order, desc(sign))
write.csv(sigDCdfpm, "sigDCdfpm_sheet.csv")

abaVSwt <- read_tsv("wt.tsv") %>%
  left_join(read_tsv("aba.tsv"), by=c("gene" = "gene"))
abaVSwt$wtavg <- rowMeans(abaVSwt[ , c(2,3,4)], na.rm=TRUE)
abaVSwt$abaavg <- rowMeans(abaVSwt[ , c(5,6,7)], na.rm=TRUE)
abaVSwt$log2FoldChange <- log2(abaVSwt$abaavg/abaVSwt$wtavg)
aba <- abaVSwt[,c(1,10)]
aba$log2FoldChange <- as.numeric(as.character(aba$log2FoldChange))
aba <- aba[is.finite(aba$log2FoldChange),]
aba <- aba[abs(aba$log2FoldChange)>0.5,]
sigCAaba <- sigCA[, c(1,3)] %>%
  left_join((aba), by=c("gene" = "gene"))
sigCAaba$log2FoldChange.y <- as.numeric(as.character(sigCAaba$log2FoldChange.y))
sigCAaba <- sigCAaba %>% 
  rename(
    LFC_sigCA = log2FoldChange.x,
    LFC_sigaba = log2FoldChange.y
  )
sigCAaba$sign <- sigCAaba$LFC_sigCA * sigCAaba$LFC_sigaba
sigCAaba$order[sigCAaba$sign>0&sigCAaba$LFC_sigaba>0] <- 1
sigCAaba$order[sigCAaba$sign>0&sigCAaba$LFC_sigaba<0] <- 2
sigCAaba$order[sigCAaba$sign<0&sigCAaba$LFC_sigaba>0] <- 3
sigCAaba$order[sigCAaba$sign<0&sigCAaba$LFC_sigaba<0] <- 4
sigCAaba <- sigCAaba %>% arrange(order, desc(sign))
write.csv(sigCAaba, "sigCAaba_sheet.csv")
sigBAaba <- sigBA[, c(1,3)] %>%
  left_join((aba), by=c("gene" = "gene"))
sigBAaba$log2FoldChange.y <- as.numeric(as.character(sigBAaba$log2FoldChange.y))
sigBAaba <- sigBAaba %>% 
  rename(
    LFC_sigBA = log2FoldChange.x,
    LFC_sigaba = log2FoldChange.y
  )
sigBAaba$sign <- sigBAaba$LFC_sigBA * sigBAaba$LFC_sigaba
sigBAaba$order[sigBAaba$sign>0&sigBAaba$LFC_sigaba>0] <- 1
sigBAaba$order[sigBAaba$sign>0&sigBAaba$LFC_sigaba<0] <- 2
sigBAaba$order[sigBAaba$sign<0&sigBAaba$LFC_sigaba>0] <- 3
sigBAaba$order[sigBAaba$sign<0&sigBAaba$LFC_sigaba<0] <- 4
sigBAaba <- sigBAaba %>% arrange(order, desc(sign))
write.csv(sigBAaba, "sigBAaba_sheet.csv")
sigDBaba <- sigDB[, c(1,3)] %>%
  left_join((aba), by=c("gene" = "gene"))
sigDBaba$log2FoldChange.y <- as.numeric(as.character(sigDBaba$log2FoldChange.y))
sigDBaba <- sigDBaba %>% 
  rename(
    LFC_sigDB = log2FoldChange.x,
    LFC_sigaba = log2FoldChange.y
  )
sigDBaba$sign <- sigDBaba$LFC_sigDB * sigDBaba$LFC_sigaba
sigDBaba$order[sigDBaba$sign>0&sigDBaba$LFC_sigaba>0] <- 1
sigDBaba$order[sigDBaba$sign>0&sigDBaba$LFC_sigaba<0] <- 2
sigDBaba$order[sigDBaba$sign<0&sigDBaba$LFC_sigaba>0] <- 3
sigDBaba$order[sigDBaba$sign<0&sigDBaba$LFC_sigaba<0] <- 4
sigDBaba <- sigDBaba %>% arrange(order, desc(sign))
write.csv(sigDBaba, "sigDBaba_sheet.csv")
sigDCaba <- sigDC[, c(1,3)] %>%
  left_join((aba), by=c("gene" = "gene"))
sigDCaba$log2FoldChange.y <- as.numeric(as.character(sigDCaba$log2FoldChange.y))
sigDCaba <- sigDCaba %>% 
  rename(
    LFC_sigDC = log2FoldChange.x,
    LFC_sigaba = log2FoldChange.y
  )
sigDCaba$sign <- sigDCaba$LFC_sigDC * sigDCaba$LFC_sigaba
sigDCaba$order[sigDCaba$sign>0&sigDCaba$LFC_sigaba>0] <- 1
sigDCaba$order[sigDCaba$sign>0&sigDCaba$LFC_sigaba<0] <- 2
sigDCaba$order[sigDCaba$sign<0&sigDCaba$LFC_sigaba>0] <- 3
sigDCaba$order[sigDCaba$sign<0&sigDCaba$LFC_sigaba<0] <- 4
sigDCaba <- sigDCaba %>% arrange(order, desc(sign))
write.csv(sigDCaba, "sigDCaba_sheet.csv")

sigDBsigCA <- read_csv("sigDBsigCA.csv") %>%
  left_join(sigCA[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join(sigDB[, c(1,3)], by=c("gene" = "gene"))
sigDBsigCA$log2FoldChange.x <- as.numeric(as.character(sigDBsigCA$log2FoldChange.x))
sigDBsigCA$log2FoldChange.y <- as.numeric(as.character(sigDBsigCA$log2FoldChange.y))
sigDBsigCA <- sigDBsigCA %>% 
  rename(
    LFC_sigCA = log2FoldChange.x,
    LFC_sigDB = log2FoldChange.y
  )
sigDBsigCA$sign <- sigDBsigCA$LFC_sigDB * sigDBsigCA$LFC_sigCA
sigDBsigCA$order[sigDBsigCA$sign>0&sigDBsigCA$LFC_sigCA>0] <- 1
sigDBsigCA$order[sigDBsigCA$sign>0&sigDBsigCA$LFC_sigCA<0] <- 2
sigDBsigCA$order[sigDBsigCA$sign<0&sigDBsigCA$LFC_sigCA>0] <- 3
sigDBsigCA$order[sigDBsigCA$sign<0&sigDBsigCA$LFC_sigCA<0] <- 4
sigDBsigCA <- sigDBsigCA %>% arrange(order, desc(sign))
write.csv(sigDBsigCA, "sigDBsigCA_sheet.csv")

abi5_LFC0.5 <- read_csv("abi5genes.csv")
sigCAabi5_LFC0.5 <- sigCA[, c(1,3)]%>%
  left_join((abi5_LFC0.5), by=c("gene" = "gene"))
sigCAabi5_LFC0.5 <- sigCAabi5_LFC0.5[is.finite(sigCAabi5_LFC0.5$log2FoldChange.y),]
sigCAabi5_LFC0.5$log2FoldChange.x <- as.numeric(as.character(sigCAabi5_LFC0.5$log2FoldChange.x))
sigCAabi5_LFC0.5$log2FoldChange.y <- as.numeric(as.character(sigCAabi5_LFC0.5$log2FoldChange.y))
sigCAabi5_LFC0.5 <- sigCAabi5_LFC0.5 %>% 
  rename(
    LFC_sigCA = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigCAabi5_LFC0.5$sign <- sigCAabi5_LFC0.5$LFC_sigCA * sigCAabi5_LFC0.5$LFC_sigabi5
sigCAabi5_LFC0.5$order[sigCAabi5_LFC0.5$sign>0&sigCAabi5_LFC0.5$LFC_sigabi5>0] <- 1
sigCAabi5_LFC0.5$order[sigCAabi5_LFC0.5$sign>0&sigCAabi5_LFC0.5$LFC_sigabi5<0] <- 2
sigCAabi5_LFC0.5$order[sigCAabi5_LFC0.5$sign<0&sigCAabi5_LFC0.5$LFC_sigabi5>0] <- 3
sigCAabi5_LFC0.5$order[sigCAabi5_LFC0.5$sign<0&sigCAabi5_LFC0.5$LFC_sigabi5<0] <- 4
sigCAabi5_LFC0.5 <- sigCAabi5_LFC0.5 %>% arrange(order, desc(sign))
write.csv(sigCAabi5_LFC0.5, "sigCAabi5_LFC0.5_sheet.csv")
sigBAabi5_LFC0.5 <- sigBA[, c(1,3)]%>%
  left_join((abi5_LFC0.5), by=c("gene" = "gene"))
sigBAabi5_LFC0.5 <- sigBAabi5_LFC0.5[is.finite(sigBAabi5_LFC0.5$log2FoldChange.y),]
sigBAabi5_LFC0.5$log2FoldChange.x <- as.numeric(as.character(sigBAabi5_LFC0.5$log2FoldChange.x))
sigBAabi5_LFC0.5$log2FoldChange.y <- as.numeric(as.character(sigBAabi5_LFC0.5$log2FoldChange.y))
sigBAabi5_LFC0.5 <- sigBAabi5_LFC0.5 %>% 
  rename(
    LFC_sigBA = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigBAabi5_LFC0.5$sign <- sigBAabi5_LFC0.5$LFC_sigBA * sigBAabi5_LFC0.5$LFC_sigabi5
sigBAabi5_LFC0.5$order[sigBAabi5_LFC0.5$sign>0&sigBAabi5_LFC0.5$LFC_sigabi5>0] <- 1
sigBAabi5_LFC0.5$order[sigBAabi5_LFC0.5$sign>0&sigBAabi5_LFC0.5$LFC_sigabi5<0] <- 2
sigBAabi5_LFC0.5$order[sigBAabi5_LFC0.5$sign<0&sigBAabi5_LFC0.5$LFC_sigabi5>0] <- 3
sigBAabi5_LFC0.5$order[sigBAabi5_LFC0.5$sign<0&sigBAabi5_LFC0.5$LFC_sigabi5<0] <- 4
sigBAabi5_LFC0.5 <- sigBAabi5_LFC0.5 %>% arrange(order, desc(sign))
write.csv(sigBAabi5_LFC0.5, "sigBAabi5_LFC0.5_sheet.csv")
sigDBabi5_LFC0.5 <- sigDB[, c(1,3)]%>%
  left_join((abi5_LFC0.5), by=c("gene" = "gene"))
sigDBabi5_LFC0.5 <- sigDBabi5_LFC0.5[is.finite(sigDBabi5_LFC0.5$log2FoldChange.y),]
sigDBabi5_LFC0.5$log2FoldChange.x <- as.numeric(as.character(sigDBabi5_LFC0.5$log2FoldChange.x))
sigDBabi5_LFC0.5$log2FoldChange.y <- as.numeric(as.character(sigDBabi5_LFC0.5$log2FoldChange.y))
sigDBabi5_LFC0.5 <- sigDBabi5_LFC0.5 %>% 
  rename(
    LFC_sigDB = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigDBabi5_LFC0.5$sign <- sigDBabi5_LFC0.5$LFC_sigDB * sigDBabi5_LFC0.5$LFC_sigabi5
sigDBabi5_LFC0.5$order[sigDBabi5_LFC0.5$sign>0&sigDBabi5_LFC0.5$LFC_sigabi5>0] <- 1
sigDBabi5_LFC0.5$order[sigDBabi5_LFC0.5$sign>0&sigDBabi5_LFC0.5$LFC_sigabi5<0] <- 2
sigDBabi5_LFC0.5$order[sigDBabi5_LFC0.5$sign<0&sigDBabi5_LFC0.5$LFC_sigabi5>0] <- 3
sigDBabi5_LFC0.5$order[sigDBabi5_LFC0.5$sign<0&sigDBabi5_LFC0.5$LFC_sigabi5<0] <- 4
sigDBabi5_LFC0.5 <- sigDBabi5_LFC0.5 %>% arrange(order, desc(sign))
write.csv(sigDBabi5_LFC0.5, "sigDBabi5_LFC0.5_sheet.csv")
sigDCabi5_LFC0.5 <- sigDC[, c(1,3)]%>%
  left_join((abi5_LFC0.5), by=c("gene" = "gene"))
sigDCabi5_LFC0.5 <- sigDCabi5_LFC0.5[is.finite(sigDCabi5_LFC0.5$log2FoldChange.y),]
sigDCabi5_LFC0.5$log2FoldChange.x <- as.numeric(as.character(sigDCabi5_LFC0.5$log2FoldChange.x))
sigDCabi5_LFC0.5$log2FoldChange.y <- as.numeric(as.character(sigDCabi5_LFC0.5$log2FoldChange.y))
sigDCabi5_LFC0.5 <- sigDCabi5_LFC0.5 %>% 
  rename(
    LFC_sigDC = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigDCabi5_LFC0.5$sign <- sigDCabi5_LFC0.5$LFC_sigDC * sigDCabi5_LFC0.5$LFC_sigabi5
sigDCabi5_LFC0.5$order[sigDCabi5_LFC0.5$sign>0&sigDCabi5_LFC0.5$LFC_sigabi5>0] <- 1
sigDCabi5_LFC0.5$order[sigDCabi5_LFC0.5$sign>0&sigDCabi5_LFC0.5$LFC_sigabi5<0] <- 2
sigDCabi5_LFC0.5$order[sigDCabi5_LFC0.5$sign<0&sigDCabi5_LFC0.5$LFC_sigabi5>0] <- 3
sigDCabi5_LFC0.5$order[sigDCabi5_LFC0.5$sign<0&sigDCabi5_LFC0.5$LFC_sigabi5<0] <- 4
sigDCabi5_LFC0.5 <- sigDCabi5_LFC0.5 %>% arrange(order, desc(sign))
write.csv(sigDCabi5_LFC0.5, "sigDCabi5_LFC0.5_sheet.csv")


abi5 <- read_csv("abi5genes.csv")
sigCAabi5 <- read_csv("sigCAabi5.csv") %>%
  left_join(sigCA[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((abi5), by=c("gene" = "gene"))
sigCAabi5$log2FoldChange.x <- as.numeric(as.character(sigCAabi5$log2FoldChange.x))
sigCAabi5$log2FoldChange.y <- as.numeric(as.character(sigCAabi5$log2FoldChange.y))
sigCAabi5 <- sigCAabi5 %>% 
  rename(
    LFC_sigCA = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigCAabi5$sign <- sigCAabi5$LFC_sigCA * sigCAabi5$LFC_sigabi5
sigCAabi5$order[sigCAabi5$sign>0&sigCAabi5$LFC_sigabi5>0] <- 1
sigCAabi5$order[sigCAabi5$sign>0&sigCAabi5$LFC_sigabi5<0] <- 2
sigCAabi5$order[sigCAabi5$sign<0&sigCAabi5$LFC_sigabi5>0] <- 3
sigCAabi5$order[sigCAabi5$sign<0&sigCAabi5$LFC_sigabi5<0] <- 4
sigCAabi5 <- sigCAabi5 %>% arrange(order, desc(sign))
write.csv(sigCAabi5, "sigCAabi5_sheet.csv")

sigBAabi5 <- read_csv("sigBAabi5.csv") %>%
  left_join(sigBA[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((abi5), by=c("gene" = "gene"))
sigBAabi5$log2FoldChange.x <- as.numeric(as.character(sigBAabi5$log2FoldChange.x))
sigBAabi5$log2FoldChange.y <- as.numeric(as.character(sigBAabi5$log2FoldChange.y))
sigBAabi5 <- sigBAabi5 %>% 
  rename(
    LFC_sigBA = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigBAabi5$sign <- sigBAabi5$LFC_sigBA * sigBAabi5$LFC_sigabi5
sigBAabi5$order[sigBAabi5$sign>0&sigBAabi5$LFC_sigabi5>0] <- 1
sigBAabi5$order[sigBAabi5$sign>0&sigBAabi5$LFC_sigabi5<0] <- 2
sigBAabi5$order[sigBAabi5$sign<0&sigBAabi5$LFC_sigabi5>0] <- 3
sigBAabi5$order[sigBAabi5$sign<0&sigBAabi5$LFC_sigabi5<0] <- 4
sigBAabi5 <- sigBAabi5 %>% arrange(order, desc(sign))
write.csv(sigBAabi5, "sigBAabi5_sheet.csv")

sigDBabi5 <- read_csv("sigDBabi5.csv") %>%
  left_join(sigDB[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((abi5), by=c("gene" = "gene"))
sigDBabi5$log2FoldChange.x <- as.numeric(as.character(sigDBabi5$log2FoldChange.x))
sigDBabi5$log2FoldChange.y <- as.numeric(as.character(sigDBabi5$log2FoldChange.y))
sigDBabi5 <- sigDBabi5 %>% 
  rename(
    LFC_sigDB = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigDBabi5$sign <- sigDBabi5$LFC_sigDB * sigDBabi5$LFC_sigabi5
sigDBabi5$order[sigDBabi5$sign>0&sigDBabi5$LFC_sigabi5>0] <- 1
sigDBabi5$order[sigDBabi5$sign>0&sigDBabi5$LFC_sigabi5<0] <- 2
sigDBabi5$order[sigDBabi5$sign<0&sigDBabi5$LFC_sigabi5>0] <- 3
sigDBabi5$order[sigDBabi5$sign<0&sigDBabi5$LFC_sigabi5<0] <- 4
sigDBabi5 <- sigDBabi5 %>% arrange(order, desc(sign))
write.csv(sigDBabi5, "sigDBabi5_sheet.csv")

sigDCabi5 <- read_csv("sigDCabi5.csv") %>%
  left_join(sigDC[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((abi5), by=c("gene" = "gene"))
sigDCabi5$log2FoldChange.x <- as.numeric(as.character(sigDCabi5$log2FoldChange.x))
sigDCabi5$log2FoldChange.y <- as.numeric(as.character(sigDCabi5$log2FoldChange.y))
sigDCabi5 <- sigDCabi5 %>% 
  rename(
    LFC_sigDC = log2FoldChange.x,
    LFC_sigabi5 = log2FoldChange.y
  )
sigDCabi5$sign <- sigDCabi5$LFC_sigDC * sigDCabi5$LFC_sigabi5
sigDCabi5$order[sigDCabi5$sign>0&sigDCabi5$LFC_sigabi5>0] <- 1
sigDCabi5$order[sigDCabi5$sign>0&sigDCabi5$LFC_sigabi5<0] <- 2
sigDCabi5$order[sigDCabi5$sign<0&sigDCabi5$LFC_sigabi5>0] <- 3
sigDCabi5$order[sigDCabi5$sign<0&sigDCabi5$LFC_sigabi5<0] <- 4
sigDCabi5 <- sigDCabi5 %>% arrange(order, desc(sign))
write.csv(sigDCabi5, "sigDCabi5_sheet.csv")

supp <- read_csv("supp.csv")
sigCAsupp <- read_csv("sigCAsupp.csv") %>%
  left_join(sigCA[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((supp), by=c("gene" = "gene"))
sigCAsupp$log2FoldChange.x <- as.numeric(as.character(sigCAsupp$log2FoldChange.x))
sigCAsupp$log2FoldChange.y <- as.numeric(as.character(sigCAsupp$log2FoldChange.y))
sigCAsupp <- sigCAsupp %>% 
  rename(
    LFC_sigCA = log2FoldChange.x,
    LFC_sigSupp = log2FoldChange.y
  )
sigCAsupp$sign <- sigCAsupp$LFC_sigCA * sigCAsupp$LFC_sigSupp
sigCAsupp$order[sigCAsupp$sign>0&sigCAsupp$LFC_sigSupp>0] <- 1
sigCAsupp$order[sigCAsupp$sign>0&sigCAsupp$LFC_sigSupp<0] <- 2
sigCAsupp$order[sigCAsupp$sign<0&sigCAsupp$LFC_sigSupp>0] <- 3
sigCAsupp$order[sigCAsupp$sign<0&sigCAsupp$LFC_sigSupp<0] <- 4
sigCAsupp <- sigCAsupp %>% arrange(order, desc(sign))
write.csv(sigCAsupp, "sigCAsupp_sheet.csv")

sigBAsupp <- read_csv("sigBAsupp.csv") %>%
  left_join(sigBA[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((supp), by=c("gene" = "gene"))
sigBAsupp$log2FoldChange.x <- as.numeric(as.character(sigBAsupp$log2FoldChange.x))
sigBAsupp$log2FoldChange.y <- as.numeric(as.character(sigBAsupp$log2FoldChange.y))
sigBAsupp <- sigBAsupp %>% 
  rename(
    LFC_sigBA = log2FoldChange.x,
    LFC_sigSupp = log2FoldChange.y
  )
sigBAsupp$sign <- sigBAsupp$LFC_sigBA * sigBAsupp$LFC_sigSupp
sigBAsupp$order[sigBAsupp$sign>0&sigBAsupp$LFC_sigSupp>0] <- 1
sigBAsupp$order[sigBAsupp$sign>0&sigBAsupp$LFC_sigSupp<0] <- 2
sigBAsupp$order[sigBAsupp$sign<0&sigBAsupp$LFC_sigSupp>0] <- 3
sigBAsupp$order[sigBAsupp$sign<0&sigBAsupp$LFC_sigSupp<0] <- 4
sigBAsupp <- sigBAsupp %>% arrange(order, desc(sign))
write.csv(sigBAsupp, "sigBAsupp_sheet.csv")

sigDBsupp <- read_csv("sigDBsupp.csv") %>%
  left_join(sigDB[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((supp), by=c("gene" = "gene"))
sigDBsupp$log2FoldChange.x <- as.numeric(as.character(sigDBsupp$log2FoldChange.x))
sigDBsupp$log2FoldChange.y <- as.numeric(as.character(sigDBsupp$log2FoldChange.y))
sigDBsupp <- sigDBsupp %>% 
  rename(
    LFC_sigDB = log2FoldChange.x,
    LFC_sigSupp = log2FoldChange.y
  )
sigDBsupp$sign <- sigDBsupp$LFC_sigDB * sigDBsupp$LFC_sigSupp
sigDBsupp$order[sigDBsupp$sign>0&sigDBsupp$LFC_sigSupp>0] <- 1
sigDBsupp$order[sigDBsupp$sign>0&sigDBsupp$LFC_sigSupp<0] <- 2
sigDBsupp$order[sigDBsupp$sign<0&sigDBsupp$LFC_sigSupp>0] <- 3
sigDBsupp$order[sigDBsupp$sign<0&sigDBsupp$LFC_sigSupp<0] <- 4
sigDBsupp <- sigDBsupp %>% arrange(order, desc(sign))
write.csv(sigDBsupp, "sigDBsupp_sheet.csv")

sigDCsupp <- read_csv("sigDCsupp.csv") %>%
  left_join(sigDC[, c(1,3)], by=c("gene" = "gene")) %>%
  left_join((supp), by=c("gene" = "gene"))
sigDCsupp$log2FoldChange.x <- as.numeric(as.character(sigDCsupp$log2FoldChange.x))
sigDCsupp$log2FoldChange.y <- as.numeric(as.character(sigDCsupp$log2FoldChange.y))
sigDCsupp <- sigDCsupp %>% 
  rename(
    LFC_sigDC = log2FoldChange.x,
    LFC_sigSupp = log2FoldChange.y
  )
sigDCsupp$sign <- sigDCsupp$LFC_sigDC * sigDCsupp$LFC_sigSupp
sigDCsupp$order[sigDCsupp$sign>0&sigDCsupp$LFC_sigSupp>0] <- 1
sigDCsupp$order[sigDCsupp$sign>0&sigDCsupp$LFC_sigSupp<0] <- 2
sigDCsupp$order[sigDCsupp$sign<0&sigDCsupp$LFC_sigSupp>0] <- 3
sigDCsupp$order[sigDCsupp$sign<0&sigDCsupp$LFC_sigSupp<0] <- 4
sigDCsupp <- sigDCsupp %>% arrange(order, desc(sign))
write.csv(sigDCsupp, "sigDCsupp_sheet.csv")

sigCA_bins <- read_csv("mapman master.csv") %>%
  inner_join(sigCA_mapman, by=c("gene" = "gene"))
write.csv(sigCA_bins, "sigCA_bins.csv")

torgenesbackup <- torgenes
torgenes <- read_csv("torgenes.csv")
torgenes <- as.data.frame(sapply(torgenes, toupper))

torsigDC <- read_csv("torsigDC.csv") %>%
  left_join(torgenes, by=c("gene" = "gene")) %>%
  left_join(sigDC[, c(1,3)], by=c("gene" = "gene"))
torsigDC <- torsigDC[, c(1,3,2)]
torsigDC$log2FoldChange.x <- as.numeric(as.character(torsigDC$log2FoldChange.x))
torsigDC$log2FoldChange.y <- as.numeric(as.character(torsigDC$log2FoldChange.y))
torsigDC <- torsigDC %>% 
  rename(
    LFC_tor = log2FoldChange.x,
    LFC_sigDC = log2FoldChange.y
  )
torsigDC$LFC_tor <- torsigDC$LFC_tor * -1
torsigDC$sign <- torsigDC$LFC_sigDC * torsigDC$LFC_tor
torsigDC$order[torsigDC$sign>0&torsigDC$LFC_tor>0] <- 1
torsigDC$order[torsigDC$sign>0&torsigDC$LFC_tor<0] <- 2
torsigDC$order[torsigDC$sign<0&torsigDC$LFC_tor>0] <- 3
torsigDC$order[torsigDC$sign<0&torsigDC$LFC_tor<0] <- 4
torsigDC <- torsigDC %>% arrange(order, desc(sign))

torsigBA <- read_csv("torsigBA.csv") %>%
  left_join(torgenes, by=c("gene" = "gene")) %>%
  left_join(sigBA[, c(1,3)], by=c("gene" = "gene"))
torsigBA <- torsigBA[, c(1,3,2)]
torsigBA$log2FoldChange.x <- as.numeric(as.character(torsigBA$log2FoldChange.x))
torsigBA$log2FoldChange.y <- as.numeric(as.character(torsigBA$log2FoldChange.y))
torsigBA <- torsigBA %>% 
  rename(
    LFC_tor = log2FoldChange.x,
    LFC_sigBA = log2FoldChange.y
  )
torsigBA$LFC_tor <- torsigBA$LFC_tor * -1
torsigBA$sign <- torsigBA$LFC_sigBA * torsigBA$LFC_tor
torsigBA$order[torsigBA$sign>0&torsigBA$LFC_tor>0] <- 1
torsigBA$order[torsigBA$sign>0&torsigBA$LFC_tor<0] <- 2
torsigBA$order[torsigBA$sign<0&torsigBA$LFC_tor>0] <- 3
torsigBA$order[torsigBA$sign<0&torsigBA$LFC_tor<0] <- 4
torsigBA <- torsigBA %>% arrange(order, desc(sign))

torsigCA <- read_csv("torsigCA.csv") %>%
  left_join(torgenes, by=c("gene" = "gene")) %>%
  left_join(sigCA[, c(1,3)], by=c("gene" = "gene"))
torsigCA <- torsigCA[, c(1,3,2)]
torsigCA$log2FoldChange.x <- as.numeric(as.character(torsigCA$log2FoldChange.x))
torsigCA$log2FoldChange.y <- as.numeric(as.character(torsigCA$log2FoldChange.y))
torsigCA <- torsigCA %>% 
  rename(
    LFC_tor = log2FoldChange.x,
    LFC_sigCA = log2FoldChange.y
  )
torsigCA$LFC_tor <- torsigCA$LFC_tor * -1
torsigCA$sign <- torsigCA$LFC_sigCA * torsigCA$LFC_tor
torsigCA$order[torsigCA$sign>0&torsigCA$LFC_tor>0] <- 1
torsigCA$order[torsigCA$sign>0&torsigCA$LFC_tor<0] <- 2
torsigCA$order[torsigCA$sign<0&torsigCA$LFC_tor>0] <- 3
torsigCA$order[torsigCA$sign<0&torsigCA$LFC_tor<0] <- 4
torsigCA <- torsigCA %>% arrange(order, desc(sign))

torsigDB <- read_csv("torsigDB.csv") %>%
  left_join(torgenes, by=c("gene" = "gene")) %>%
  left_join(sigDB[, c(1,3)], by=c("gene" = "gene"))
torsigDB <- torsigDB[, c(1,3,2)]
torsigDB$log2FoldChange.x <- as.numeric(as.character(torsigDB$log2FoldChange.x))
torsigDB$log2FoldChange.y <- as.numeric(as.character(torsigDB$log2FoldChange.y))
torsigDB <- torsigDB %>% 
  rename(
    LFC_tor = log2FoldChange.x,
    LFC_sigDB = log2FoldChange.y
  )
torsigDB$LFC_tor <- torsigDB$LFC_tor * -1
torsigDB$sign <- torsigDB$LFC_sigDB * torsigDB$LFC_tor
torsigDB$order[torsigDB$sign>0&torsigDB$LFC_tor>0] <- 1
torsigDB$order[torsigDB$sign>0&torsigDB$LFC_tor<0] <- 2
torsigDB$order[torsigDB$sign<0&torsigDB$LFC_tor>0] <- 3
torsigDB$order[torsigDB$sign<0&torsigDB$LFC_tor<0] <- 4
torsigDB <- torsigDB %>% arrange(order, desc(sign))

write.csv(torsigDC, "torsigDC_sheet.csv")
write.csv(torsigCA, "torsigCA_sheet.csv")
write.csv(torsigDB, "torsigDB_sheet.csv")
write.csv(torsigBA, "torsigBA_sheet.csv")

torsigDClimited <- torsigDC
torsigDClimited[2] <- replace(torsigDClimited[2], torsigDClimited[2] > 5, 5)
torsigDClimited[3] <- replace(torsigDClimited[3], torsigDClimited[3] > 5, 5)
torsigDClimited[2] <- replace(torsigDClimited[2], torsigDClimited[2] < -5, -5)
torsigDClimited[3] <- replace(torsigDClimited[3], torsigDClimited[3] < -5, -5)
torsigDClimited <- torsigDClimited[ , c(1,3,2,4)]
torsigDClimited <- torsigDClimited %>% arrange(desc(LFC_tor))
torsigDClimited_invert <- torsigDClimited
torsigDClimited_invert[3] <- torsigDClimited_invert[3]*-1

pheatmap::pheatmap(torsigDBlimited_invert[2:3],
                   cluster_rows = T,
                   show_rownames = F,
                   color = heat_colors,
                   border_color = NA,
                   filename = "tor_sigDB_invert.png",
                   treeheight_row = 0, 
                   treeheight_col = 0,
                   cluster_cols=FALSE)

pheatmap::pheatmap(new2_total_norm[2:13],
                   color = heat_colors,
                   cluster_rows = T,
                   show_rownames = F,
                   annotation = meta,
                   border_color = NA,
                   fontsize = 10,
                   scale = "row",
                   fontsize_row = 10,
                   filename = "results/total/completeheatmap_NOzeros.png",
                   cluster_cols=FALSE,
                   height = 15)

pheatmap::pheatmap(sigCA_norm[2:13],
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = meta,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         filename = "results/CA/heatmaps/CA_heatmaptotal.png",
         cluster_cols=FALSE,
         height = 15)

pheatmap::pheatmap(sigBA_norm[2:13],
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = meta,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         filename = "results/BA/heatmaps/BA_heatmaptotal.png",
         cluster_cols=FALSE,
         height = 15)

pheatmap::pheatmap(sigDA_norm[2:13],
                   color = heat_colors,
                   cluster_rows = T,
                   show_rownames = F,
                   annotation = meta,
                   border_color = NA,
                   fontsize = 10,
                   scale = "row",
                   fontsize_row = 10,
                   filename = "results/DA/heatmaps/DA_heatmaptotal.png",
                   cluster_cols=FALSE,
                   height = 15)

pheatmap::pheatmap(sigDB_norm[2:13],
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = meta,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         filename = "results/DB/heatmaps/DB_heatmaptotal.png",
         cluster_cols=FALSE,
         height = 15)

pheatmap::pheatmap(sigCB_norm[2:13],
                   color = heat_colors,
                   cluster_rows = T,
                   show_rownames = F,
                   annotation = meta,
                   border_color = NA,
                   fontsize = 10,
                   scale = "row",
                   fontsize_row = 10,
                   filename = "results/CB/heatmaps/CB_heatmaptotal.png",
                   cluster_cols=FALSE,
                   height = 15)

pheatmap::pheatmap(sigDC_norm[2:13],
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = meta,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         filename = "results/DC/heatmaps/DC_heatmaptotal.png",
         cluster_cols=FALSE,
         height = 15)

########plot volcanos
res_tableCA_tb <- res_tableCA_tb %>% 
  mutate(threshold_CA = padj < 0.05)
res_tableCA_tb <- res_tableCA_tb %>% mutate(genelabels = "")
res_tableCA_tb <- res_tableCA_tb %>% arrange(padj)
res_tableCA_tb$genelabels[1:10] <- as.character(res_tableCA_tb$symbol[1:10])
ggplot(res_tableCA_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_CA)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("C vs A") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave("results/CA/CA_volcano.png")

res_tableBA_tb <- res_tableBA_tb %>% 
  mutate(threshold_BA = padj < 0.05)
res_tableBA_tb <- res_tableBA_tb %>% mutate(genelabels = "")
res_tableBA_tb <- res_tableBA_tb %>% arrange(padj)
res_tableBA_tb$genelabels[1:10] <- as.character(res_tableBA_tb$symbol[1:10])
ggplot(res_tableBA_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_BA)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("B vs A") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave("results/BA/BA_volcano.png")

res_tableDA_tb <- res_tableDA_tb %>% 
  mutate(threshold_DA = padj < 0.05)
res_tableDA_tb <- res_tableDA_tb %>% mutate(genelabels = "")
res_tableDA_tb <- res_tableDA_tb %>% arrange(padj)
res_tableDA_tb$genelabels[1:10] <- as.character(res_tableDA_tb$symbol[1:10])
ggplot(res_tableDA_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_DA)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("D vs A") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave("results/DA/DA_volcano.png")

res_tableDB_tb <- res_tableDB_tb %>% 
  mutate(threshold_DB = padj < 0.05)
res_tableDB_tb <- res_tableDB_tb %>% mutate(genelabels = "")
res_tableDB_tb <- res_tableDB_tb %>% arrange(padj)
res_tableDB_tb$genelabels[1:10] <- as.character(res_tableDB_tb$symbol[1:10])
ggplot(res_tableDB_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_DB)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("D vs B") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave("results/DB/DB_volcano.png")

res_tableCB_tb <- res_tableCB_tb %>% 
  mutate(threshold_CB = padj < 0.05)
res_tableCB_tb <- res_tableCB_tb %>% mutate(genelabels = "")
res_tableCB_tb <- res_tableCB_tb %>% arrange(padj)
res_tableCB_tb$genelabels[1:10] <- as.character(res_tableCB_tb$symbol[1:10])
ggplot(res_tableCB_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_CB)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("C vs B") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave("results/CB/CB_volcano.png")

res_tableDC_tb <- res_tableDC_tb %>% 
  mutate(threshold_DC = padj < 0.05)
res_tableDC_tb <- res_tableDC_tb %>% mutate(genelabels = "")
res_tableDC_tb <- res_tableDC_tb %>% arrange(padj)
res_tableDC_tb$genelabels[1:10] <- as.character(res_tableDC_tb$symbol[1:10])
ggplot(res_tableDC_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_DC)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("D vs C") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
ggsave("results/DC/DC_volcano.png")


# ########likelihood ratio test 
# dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
# res_LRT <- results(dds_lrt)
# 
# # Create a tibble for LRT results
# res_LRT_tb <- res_LRT %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble()
# 
# # Subset to return genes with padj <
# sigLRT_genes <- res_LRT_tb %>% 
#   filter(padj < 0.001)
# 
# # Subset results for faster cluster finding 
# clustering_sig_genes <- sigLRT_genes %>%
#   arrange(padj) %>%
#   head(n=10000)
# 
# # Obtain rlog values for those significant genes
# cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
# 
# # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
# meta2 <- read.table("data/testMeta.txt", header=T, row.names=1)
# clusters <- degPatterns(cluster_rlog, metadata = meta2, time = "Group", col="Type")

# # Extract the Group 1 genes
# cluster_groups <- clusters$df
# group1 <- clusters$df %>%
#   filter(cluster == 1)

all_genes <- as.character(res_tableCA_tb$gene) #background dataset for hypergeometric testing
sigCA_genes <- as.character(sigCA_norm$gene) #significant genes dataset
sigBA_genes <- as.character(sigBA_norm$gene)
sigDA_genes <- as.character(sigDA_norm$gene)
sigCB_genes <- as.character(sigCB_norm$gene)
sigDB_genes <- as.character(sigDB_norm$gene)
sigDC_genes <- as.character(sigDC_norm$gene)

BA1_DC1_merged_genes <- as.character(BA1_DC1_merged$gene)
CA2_CA3_DB9_merged_genes <- as.character(CA2_CA3_DB9_merged$gene)
DC2_CA10_merged_genes <- as.character(DC2_CA10_merged$gene)
DC4_DB6_merged_genes <- as.character(DC4_DB6_merged$gene)

## Run GO enrichment analysis 
egoDC4_DB6_merged <- enrichGO(gene = DC4_DB6_merged_genes, 
                  universe = all_genes,
                  keyType = "TAIR",
                  OrgDb = org.At.tair.db, 
                  ont = "ALL", 
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  readable = TRUE)
cluster_summary <- data.frame(egoDC4_DB6_merged)
cluster_summaryBP <- cluster_summary %>% filter(ONTOLOGY == "BP")
cluster_summaryMF <- cluster_summary %>% filter(ONTOLOGY == "MF")
write.csv(cluster_summaryBP, "bhlh48_heatmap_analysis/merged_sets/DC4_DB6_merged_BP.csv")
write.csv(cluster_summaryMF, "bhlh48_heatmap_analysis/merged_sets/DC4_DB6_merged_MF.csv")

egoCA <- enrichGO(gene = sigCA_genes, 
                universe = all_genes,
                keyType = "TAIR",
                OrgDb = org.At.tair.db, 
                ont = "ALL", 
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                pvalueCutoff = 1,
                readable = TRUE)
cluster_summary <- data.frame(egoCA)
cluster_summaryBP <- cluster_summary %>% filter(ONTOLOGY == "BP")
cluster_summaryMF <- cluster_summary %>% filter(ONTOLOGY == "MF")
write.csv(cluster_summaryBP, "results/CA/clusterProfiler_CA_BP.csv")
write.csv(cluster_summaryMF, "results/CA/clusterProfiler_CA_MF.csv")

egoBA <- enrichGO(gene = sigBA_genes, 
                  universe = all_genes,
                  keyType = "TAIR",
                  OrgDb = org.At.tair.db, 
                  ont = "ALL", 
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  readable = TRUE)
cluster_summary <- data.frame(egoBA)
cluster_summaryBP <- cluster_summary %>% filter(ONTOLOGY == "BP")
cluster_summaryMF <- cluster_summary %>% filter(ONTOLOGY == "MF")
write.csv(cluster_summaryBP, "results/BA/clusterProfiler_BA_BP.csv")
write.csv(cluster_summaryMF, "results/BA/clusterProfiler_BA_MF.csv")

egoDA <- enrichGO(gene = sigDA_genes, 
                  universe = all_genes,
                  keyType = "TAIR",
                  OrgDb = org.At.tair.db, 
                  ont = "ALL", 
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  readable = TRUE)
cluster_summary <- data.frame(egoDA)
cluster_summaryBP <- cluster_summary %>% filter(ONTOLOGY == "BP")
cluster_summaryMF <- cluster_summary %>% filter(ONTOLOGY == "MF")
write.csv(cluster_summaryBP, "results/DA/clusterProfiler_DA_BP.csv")
write.csv(cluster_summaryMF, "results/DA/clusterProfiler_DA_MF.csv")

egoCB <- enrichGO(gene = sigCB_genes, 
                  universe = all_genes,
                  keyType = "TAIR",
                  OrgDb = org.At.tair.db, 
                  ont = "ALL", 
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  readable = TRUE)
cluster_summary <- data.frame(egoCB)
cluster_summaryBP <- cluster_summary %>% filter(ONTOLOGY == "BP")
cluster_summaryMF <- cluster_summary %>% filter(ONTOLOGY == "MF")
write.csv(cluster_summaryBP, "results/CB/clusterProfiler_CB_BP.csv")
write.csv(cluster_summaryMF, "results/CB/clusterProfiler_CB_MF.csv")

egoDB <- enrichGO(gene = sigDB_genes, 
                  universe = all_genes,
                  keyType = "TAIR",
                  OrgDb = org.At.tair.db, 
                  ont = "ALL", 
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  readable = TRUE)
cluster_summary <- data.frame(egoDB)
cluster_summaryBP <- cluster_summary %>% filter(ONTOLOGY == "BP")
cluster_summaryMF <- cluster_summary %>% filter(ONTOLOGY == "MF")
write.csv(cluster_summaryBP, "results/DB/clusterProfiler_DB_BP.csv")
write.csv(cluster_summaryMF, "results/DB/clusterProfiler_DB_MF.csv")

egoDC <- enrichGO(gene = sigDC_genes, 
                  universe = all_genes,
                  keyType = "TAIR",
                  OrgDb = org.At.tair.db, 
                  ont = "ALL", 
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  pvalueCutoff = 1,
                  readable = TRUE)
cluster_summary <- data.frame(egoDC)
cluster_summaryBP <- cluster_summary %>% filter(ONTOLOGY == "BP")
cluster_summaryMF <- cluster_summary %>% filter(ONTOLOGY == "MF")
egoDC@result$geneID 
write.csv(cluster_summaryBP, "results/DC/clusterProfiler_DC_BP.csv")
write.csv(cluster_summaryMF, "results/DC/clusterProfiler_DC_MF.csv")

selectedCategories <- c("regulation of abscisic acid biosynthetic process", 
                        "abscisic acid metabolic process", 
                        "cellular response to abscisic acid stimulus",
                        "abscisic acid-activated signaling pathway",
                        "abscisic acid biosynthetic process",
                        "regulation of DNA-binding transcription factor activity")

# ## plots
# dotplot(ego, showCategory=selectedCategories)
# ego2<-pairwise_termsim(ego)
# emapplot(ego2, showCategory=selectedCategories)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
CA_foldchanges <- sigCA$log2FoldChange
names(CA_foldchanges) <- sigCA$gene
CA_foldchanges <- ifelse(CA_foldchanges > 2, 2, CA_foldchanges)
CA_foldchanges <- ifelse(CA_foldchanges < -2, -2, CA_foldchanges)
pdf("results/CA/CA_cnetplot.pdf", width = 20, height = 20)
plot <- cnetplot(egoCA, 
         categorySize="pvalue", 
         showCategory = selectedCategories, 
         foldChange=CA_foldchanges)
plot
dev.off()

BA_foldchanges <- sigBA$log2FoldChange
names(BA_foldchanges) <- sigBA$gene
BA_foldchanges <- ifelse(BA_foldchanges > 2, 2, BA_foldchanges)
BA_foldchanges <- ifelse(BA_foldchanges < -2, -2, BA_foldchanges)
pdf("results/BA/BA_cnetplot.pdf", width = 20, height = 20)
plot <- cnetplot(egoBA, 
                 categorySize="pvalue", 
                 showCategory = selectedCategories, 
                 foldChange=BA_foldchanges)
plot
dev.off()

DA_foldchanges <- sigDA$log2FoldChange
names(DA_foldchanges) <- sigDA$gene
DA_foldchanges <- ifelse(DA_foldchanges > 2, 2, DA_foldchanges)
DA_foldchanges <- ifelse(DA_foldchanges < -2, -2, DA_foldchanges)
pdf("results/DA/DA_cnetplot.pdf", width = 20, height = 20)
plot <- cnetplot(egoDA, 
                 categorySize="pvalue", 
                 showCategory = selectedCategories, 
                 foldChange=DA_foldchanges)
plot
dev.off()

CB_foldchanges <- sigCB$log2FoldChange
names(CB_foldchanges) <- sigCB$gene
CB_foldchanges <- ifelse(CB_foldchanges > 2, 2, CB_foldchanges)
CB_foldchanges <- ifelse(CB_foldchanges < -2, -2, CB_foldchanges)
pdf("results/CB/CB_cnetplot.pdf", width = 20, height = 20)
plot <- cnetplot(egoCB, 
                 categorySize="pvalue", 
                 showCategory = selectedCategories, 
                 foldChange=CB_foldchanges)
plot
dev.off()

DB_foldchanges <- sigDB$log2FoldChange
names(DB_foldchanges) <- sigDB$gene
DB_foldchanges <- ifelse(DB_foldchanges > 2, 2, DB_foldchanges)
DB_foldchanges <- ifelse(DB_foldchanges < -2, -2, DB_foldchanges)
pdf("results/DB/DB_cnetplot.pdf", width = 20, height = 20)
plot <- cnetplot(egoDB, 
                 categorySize="pvalue", 
                 showCategory = selectedCategories, 
                 foldChange=DB_foldchanges)
plot
dev.off()

DC_foldchanges <- sigDC$log2FoldChange
names(DC_foldchanges) <- sigDC$gene
DC_foldchanges <- ifelse(DC_foldchanges > 2, 2, DC_foldchanges)
DC_foldchanges <- ifelse(DC_foldchanges < -2, -2, DC_foldchanges)
pdf("results/DC/DC_cnetplot.pdf", width = 20, height = 20)
plot <- cnetplot(egoDC, 
                 categorySize="pvalue", 
                 showCategory = selectedCategories, 
                 foldChange=DC_foldchanges)
plot
dev.off()

write.csv(sigBA1, "bhlh48_heatmap_analysis/BA/sigBA1.csv")
write.csv(sigBA2, "bhlh48_heatmap_analysis/BA/sigBA2.csv")
write.csv(sigBA3, "bhlh48_heatmap_analysis/BA/sigBA3.csv")
write.csv(sigBA4, "bhlh48_heatmap_analysis/BA/sigBA4.csv")
write.csv(sigBA5, "bhlh48_heatmap_analysis/BA/sigBA5.csv")
write.csv(sigBA6, "bhlh48_heatmap_analysis/BA/sigBA6.csv")
write.csv(sigBA7, "bhlh48_heatmap_analysis/BA/sigBA7.csv")
write.csv(sigBA8, "bhlh48_heatmap_analysis/BA/sigBA8.csv")
write.csv(sigBA9, "bhlh48_heatmap_analysis/BA/sigBA9.csv")
write.csv(sigBA10, "bhlh48_heatmap_analysis/BA/sigBA10.csv")

write.csv(sigCA1, "bhlh48_heatmap_analysis/CA/sigCA1.csv")
write.csv(sigCA2, "bhlh48_heatmap_analysis/CA/sigCA2.csv")
write.csv(sigCA3, "bhlh48_heatmap_analysis/CA/sigCA3.csv")
write.csv(sigCA4, "bhlh48_heatmap_analysis/CA/sigCA4.csv")
write.csv(sigCA5, "bhlh48_heatmap_analysis/CA/sigCA5.csv")
write.csv(sigCA6, "bhlh48_heatmap_analysis/CA/sigCA6.csv")
write.csv(sigCA7, "bhlh48_heatmap_analysis/CA/sigCA7.csv")
write.csv(sigCA8, "bhlh48_heatmap_analysis/CA/sigCA8.csv")
write.csv(sigCA9, "bhlh48_heatmap_analysis/CA/sigCA9.csv")
write.csv(sigCA10, "bhlh48_heatmap_analysis/CA/sigCA10.csv")

write.csv(sigDB1, "bhlh48_heatmap_analysis/DB/sigDB1.csv")
write.csv(sigDB2, "bhlh48_heatmap_analysis/DB/sigDB2.csv")
write.csv(sigDB3, "bhlh48_heatmap_analysis/DB/sigDB3.csv")
write.csv(sigDB4, "bhlh48_heatmap_analysis/DB/sigDB4.csv")
write.csv(sigDB5, "bhlh48_heatmap_analysis/DB/sigDB5.csv")
write.csv(sigDB6, "bhlh48_heatmap_analysis/DB/sigDB6.csv")
write.csv(sigDB7, "bhlh48_heatmap_analysis/DB/sigDB7.csv")
write.csv(sigDB8, "bhlh48_heatmap_analysis/DB/sigDB8.csv")
write.csv(sigDB9, "bhlh48_heatmap_analysis/DB/sigDB9.csv")

write.csv(sigDC1, "bhlh48_heatmap_analysis/DC/sigDC1.csv")
write.csv(sigDC2, "bhlh48_heatmap_analysis/DC/sigDC2.csv")
write.csv(sigDC3, "bhlh48_heatmap_analysis/DC/sigDC3.csv")
write.csv(sigDC4, "bhlh48_heatmap_analysis/DC/sigDC4.csv")

BA1_DC1 <- inner_join(sigBA1, sigDC1)
BA1_DC1_merged <- full_join(sigBA1, sigDC1)

CA2_CA3 <- full_join(sigCA2, sigCA3)
CA2_CA3_DB9 <- inner_join(CA2_CA3, sigDB9)
CA2_CA3_DB9_merged <- full_join(CA2_CA3, sigDB9)

DC2_CA10 <- inner_join(sigDC2, sigCA10)
DC2_CA10_merged <- full_join(sigDC2, sigCA10)

DC4_DB6 <- inner_join(sigDC4, sigDB6)
DC4_DB6_merged <- full_join(sigDC4, sigDB6)

write.csv(BA1_DC1_merged, "bhlh48_heatmap_analysis/merged_sets/BA1_DC1_merged.csv")
write.csv(CA2_CA3_DB9_merged, "bhlh48_heatmap_analysis/merged_sets/CA2_CA3_DB9_merged.csv")
write.csv(DC2_CA10_merged, "bhlh48_heatmap_analysis/merged_sets/DC2_CA10_merged.csv")
write.csv(DC4_DB6_merged, "bhlh48_heatmap_analysis/merged_sets/DC4_DB6_merged.csv")

DC3_CA2 <- inner_join(sigDC3, sigCA2)
write.csv(DC3_CA2, "DC3_CA2.csv")

input <- read_csv("DC3_CA2_proteinmodification.csv") %>%
  left_join(sigDC, by=c("gene" = "gene")) %>%
  left_join(normalized_counts, by=c("gene" = "gene"))
input$A_mean = (input$RF01 + input$RF02 + input$RF03) / 3
input$B_mean = (input$RF04 + input$RF05 + input$RF06) / 3
input$C_mean = (input$RF07 + input$RF08 + input$RF09) / 3
input$D_mean = (input$RF10 + input$RF11 + input$RF12) / 3

DC3_CA2_proteinmodification <- input[, c(1,7,6,25,27,26,28,8)]
write.csv(DC3_CA2_proteinmodification, "DC3_CA2_proteinmodification.csv")

sigBAmerge <- sigBA[ , c(1,6)] %>% rename(padj.BA = padj)
sigCAmerge <- sigCA[ , c(1,6)] %>% rename(padj.CA = padj)
sigDBmerge <- sigDB[ , c(1,6)] %>% rename(padj.DB = padj)
sigDCmerge <- sigDC[ , c(1,6)] %>% rename(padj.DC = padj)
countmeans <- normalized_counts[ , c(1,14,16)]

countmeans$A_WTmin = (normalized_counts$RF01 + normalized_counts$RF02 + normalized_counts$RF03) / 3
countmeans$C_bHLHmin = (normalized_counts$RF07 + normalized_counts$RF08 + normalized_counts$RF09) / 3
countmeans$B_WTglu = (normalized_counts$RF04 + normalized_counts$RF05 + normalized_counts$RF06) / 3
countmeans$D_bHLHglu = (normalized_counts$RF10 + normalized_counts$RF11 + normalized_counts$RF12) / 3

mergedlist <- read_csv("input.csv") %>%
  left_join(sigBAmerge, by=c("gene" = "gene")) %>%
  left_join(sigCAmerge, by=c("gene" = "gene")) %>%
  left_join(sigDBmerge, by=c("gene" = "gene")) %>%
  left_join(sigDCmerge, by=c("gene" = "gene")) %>%
  left_join(countmeans, by=c("gene" = "gene"))
mergedlist <- mergedlist[ , c(1,6,2,3,4,5,8,9,10,11,7)]
write.csv(mergedlist, "input.csv")

sigBAaba_gse28800 <- left_join(sigBAaba_gse28800, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigBAaba_gse28800, "sigBAaba_gse28800.csv")
sigBAaba <- left_join(sigBAaba, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigBAaba, "sigBAaba.csv")
sigBAabi5 <- left_join(sigBAabi5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigBAabi5, "sigBAabi5_LFC_1.csv")
sigBAabi5_LFC0.5 <- left_join(sigBAabi5_LFC0.5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigBAabi5_LFC0.5, "sigBAabi5_LFC_0.5.csv")
sigBAdfpm <- left_join(sigBAdfpm, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigBAdfpm, "sigBAdfpm.csv")
torsigBA <- left_join(torsigBA, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(torsigBA, "torsigBA.csv")

sigCAaba_gse28800 <- left_join(sigCAaba_gse28800, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigCAaba_gse28800, "sigCAaba_gse28800.csv")
sigCAaba <- left_join(sigCAaba, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigCAaba, "sigCAaba.csv")
sigCAabi5 <- left_join(sigCAabi5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigCAabi5, "sigCAabi5_LFC_1.csv")
sigCAabi5_LFC0.5 <- left_join(sigCAabi5_LFC0.5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigCAabi5_LFC0.5, "sigCAabi5_LFC_0.5.csv")
sigCAdfpm <- left_join(sigCAdfpm, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigCAdfpm, "sigCAdfpm.csv")
torsigCA <- left_join(torsigCA, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(torsigCA, "torsigCA.csv")

sigDBaba_gse28800 <- left_join(sigDBaba_gse28800, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDBaba_gse28800, "sigDBaba_gse28800.csv")
sigDBaba <- left_join(sigDBaba, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDBaba, "sigDBaba.csv")
sigDBabi5 <- left_join(sigDBabi5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDBabi5, "sigDBabi5_LFC_1.csv")
sigDBabi5_LFC0.5 <- left_join(sigDBabi5_LFC0.5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDBabi5_LFC0.5, "sigDBabi5_LFC_0.5.csv")
sigDBdfpm <- left_join(sigDBdfpm, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDBdfpm, "sigDBdfpm.csv")
torsigDB <- left_join(torsigDB, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(torsigDB, "torsigDB.csv")

sigDCaba_gse28800 <- left_join(sigDCaba_gse28800, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDCaba_gse28800, "sigDCaba_gse28800.csv")
sigDCaba <- left_join(sigDCaba, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDCaba, "sigDCaba.csv")
sigDCabi5 <- left_join(sigDCabi5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDCabi5, "sigDCabi5_LFC_1.csv")
sigDCabi5_LFC0.5 <- left_join(sigDCabi5_LFC0.5, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDCabi5_LFC0.5, "sigDCabi5_LFC_0.5.csv")
sigDCdfpm <- left_join(sigDCdfpm, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(sigDCdfpm, "sigDCdfpm.csv")
torsigDC <- left_join(torsigDC, resCA_entrez[,c(1,2)], by=c("gene" = "gene"))
write.csv(torsigDC, "torsigDC.csv")

all_LFC <- resBA_entrez[,c(1,7,2,3)] %>%
  left_join(resDC_entrez[,c(1,3)], by=c("gene" = "gene")) %>%
  left_join(resCA_entrez[,c(1,3)], by=c("gene" = "gene")) %>%
  left_join(resDB_entrez[,c(1,3)], by=c("gene" = "gene")) %>%
  left_join(abaVSwt[, c(1,10)], by=c("gene" = "gene")) %>% 
  left_join(dfpm[,c(9,6)], by=c("gene" = "To")) %>%
  left_join(aba_gse28800[,c(9,6)], by=c("gene" = "To")) %>%
  left_join(resCA_entrez[,c(1,6)], by=c("gene" = "gene")) %>%
  left_join(resDC_entrez[,c(1,6)], by=c("gene" = "gene")) %>%
  left_join(resCA_entrez[,c(1,6)], by=c("gene" = "gene")) %>%
  left_join(resDB_entrez[,c(1,6)], by=c("gene" = "gene")) %>%
  left_join(dfpm[,c(9,3)], by=c("gene" = "To")) %>%
  left_join(aba_gse28800[,c(9,3)], by=c("gene" = "To")) %>% 
  left_join(resBA_entrez[,c(1,8)], by=c("gene" = "gene")) %>%
  rename(
    LFC_BA = log2FoldChange.x,
    LFC_DC = log2FoldChange.y,
    LFC_CA = log2FoldChange.x.x,
    LFC_DB = log2FoldChange.y.y,
    padj_BA = padj.x,
    padj_DC = padj.y,
    padj_CA = padj.x.x,
    padj_DB = padj.y.y,
    LFC_ABA_7d = log2FoldChange,
    LFC_DFPM = logFC.x,
    LFC_ABA_12d = logFC.y,
    padj_DFPM = P.Value.x,
    padj_ABA_12d = P.Value.y)
  
