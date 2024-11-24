# GO and KEGG pathway enrichment analysis using DGE list as background/universe
# Date: 05-06-2023

# Load libraries
library(clusterProfiler)
library(enrichplot)
library(org.Cf.eg.db)

## Over-representation analysis

# Define the universe gene list
gene_universe <- rownames(results)

# Run enrichGO
## upregulated-genes
upregulated_enrichgo_out <- enrichGO(gene = upregulated_gene_list,
                                     universe = gene_universe,
                                     OrgDb = org.Cf.eg.db,
                                     ont = "ALL",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.05,
                                     readable = TRUE)

## down-regulated genes
downregulated_enrichgo_out <- enrichGO(gene = downregulated_gene_list,
                                     universe = gene_universe,
                                     OrgDb = org.Cf.eg.db,
                                     ont = "ALL",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.05,
                                     readable = TRUE)


## KEGG Enrichment Pathway
## up-regulated genes
upregulated_kegg <- enrichKEGG(gene = upregulated_gene_list,
                                   universe = gene_universe,
                                   organism = "cfa",
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.05,
                                   minGSSize = 10,
                                   maxGSSize = 500)

## down-regulated genes
downregulated_kegg <- enrichKEGG(gene = downregulated_gene_list,
                                 universe = gene_universe,
                                 organism = "cfa",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.05,
                                 minGSSize = 10,
                                 maxGSSize = 500)
