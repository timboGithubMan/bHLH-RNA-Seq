# bHLH-RNA-Seq
RNA-Seq data and analysis of a mutant overexpressing BHLH48, a transcription factor that interacts with ABI4 and is suspected to play a role in regulating glucose response.

First, Salmon was used to calculate count data from RNA seq fastQ files.
This unnormalized count data was then imported to R.

Then, DESeqDataSetFromMatrix files were created specifying countData, colData, and design type
	dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Type)
	colData tells DESeq2 information about your samples (so that you can place them into groups)
	and design tells DESeq2 how to compare your data (in this case, across Type, which represents conditions)

Then, we performed DESeq2 on the comparisons
	dds <- DESeq(dds)

Then created results tables. Testing for the likelihood of log fold change to be at least greater than 0.41
    with a false discovery rate of 0.05.
	res_tableBA <- results(dds, contrast=c("Type", "B", "A"), alpha = 0.05, lfcThreshold = 0.41)

Then, we shrunk the log fold changes
	res_tableCAlfc <- lfcShrink(dds, coef="Type_C_vs_A", res=res_tableCA)
	This filters out the noise from genes with low counts or inconsistent counts (bad data).

Then, we created easy to read tables of the results
	res_tableCA_tb <- res_tableCAlfc .........
	
Finally, we made subsets containing only genes with p<0.05 across different comparisons
	sigCA <- res_tableCA_tb %>%
	filter(padj < padj.cutoff)
	
And created a table of the counts, normalized by DESeq2
	normalized_counts
	
For K-means clustering:
First, a table of medians was created called count_medians
Then, this table was filtered for genes with p<0.05 across one comparison of interest: count_medians_filtered
Then, the rows were scaled by z-score: scaled_countsMEDIANS_filtered (or scaledata)

Summary of files and folders:

"data" 
- gene counts, information about the samples, and gene ID conversions

"mapman_analysis" 
- this contains text files formatted for use with MapMan, 
  a tool that organizes genes into different GO terms to create visual maps
  
"coexpression_analysis" 
- this contains files related to coexpression analysis,
  which compares effects of different experiments/treatments to look for overlap
  
"new_bhlh48_heatmap_analysis" 
- this is the previous kind of heatmap cluster analysis we were doing.
  in this kind of analysis, we made a heatmap for each pairwise comparison
  and then manually pulled out clusters. the heatmaps were made using all the default 
  settings (using pearson correlation, and scaled by z-score)
  
"results"
- there are a couple of files in this folder:
*all_LFC - all log fold changes, for easy comparison, including some other similar experiments/treatments
*PCA - a PCA plot showing that replicates of conditions group together (data is consistent across replicates)
*Correlation Heatmap - a map showing correlation between samples (also shows that data is consistent across replicates)
*complete_heatmap - a heatmap of all samples 
*normalized_counts.csv - gene counts normalized by DESeq2

-then there is a folder for each comparison: "CA", "DB", etc
within each folder, there are the files:
*sigCA.csv - a spreadsheet of all genes significantly different between C and A
*clusterProfiler_CA_MF (and _BP) clusterProfiler GO analysis of sigCA genes. 
 the clusterProfiler analysis is kind of outdated, there are better GO tools to use
*CA_volcano - a volcano plot of sigCA genes
*CA_cnetplot - created using the clusterProfiler results. a visualization of a couple of GO terms
"KEGG pathways" - a visualization of the regulation of different KEGG pathways across the pairwise comparison



