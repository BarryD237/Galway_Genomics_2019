Heatmaps were generated using the following code: 

rld <- rlog(dds, blind = FALSE) 

^ This function transforms the count data to the log2 scale

resOrdered <- subset(res, padj < 0.05)
resOrdered <- resOrdered[order(resOrdered$log2FoldChange),]

^ Captures all statistically significant genes (padj < 0.05) i.e the differentially expressed genes
  Ordered by log2 fold change 

select_genes <- rownames(resOrdered)
df = as.data.frame(colData(rld)[,c("celltype")])
colnames(df) = "celltype"
rownames(df) = colnames(rld)

^ select_genes extracts Ensembl gene ID's from resOrdered
  df defines the comparison (we are comparing cell types)

mat  <- assay(rld)[ select_genes, ]
mat  <- mat - rowMeans(mat)
gns <- getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name"),
             filters = c("ensembl_gene_id"),
             values = row.names(mat),
             mart = mart)
row.names(mat)[match(gns[,1], row.names(mat))] <- gns[,2]

^ Match the corresponding gene symbols to the Ensembl gene ID's using BioMart. 

outfile="Heatmap_Aria_vs_TYTO_S2_CTRL.pdf"
pdf(file=outfile)
pheatmap(mat[,-c(1,3,5,7,9,11,13,15)], border_color = NA, fontsize = 10, scale = "row", fontsize_row = 2, height = 20,  annotation_col = df, cluster_cols = TRUE, cluster_rows = TRUE,
         main="Aria_vs_TYTO_S2 (Control)")
dev.off()

^ Edit as required for the comparisons being made
  mat[,-c(n,n,n,n,n,)], selects the Aria and TYTO_S2 columns for the heatmap.
  Without this all cell types will be present in the heatmap.
  By selecting Aria and TYTO_S2, we tell R we are only interested in this groups differentially expressed genes
  Dense heatmaps = change fontsize_row = 1. 
