# import DESeq2 data from slamdunk count
DESeq2_Gad2 <- read.delim(*path to Gad2 DESeq2 file*)
DESeq2_vGluT2 <- read.delim(*path to vGluT2 DESeq2 file*)


# find genes that express in both conditions
merged <- merge(DESeq2_Gad2, DESeq2_vGluT2, by.x = "gene_name", by.y = "gene_name")

# Plot heatmap for all gene
pheatmap::pheatmap(merged[,c(2,5)], cluster_cols = F, main = "Heapmap of DEGs in both cells", angle_col = 0,
                   show_rownames = F, labels_col = c("Gad2","vGluT2"))

# Find oppositely deferentially expressed genes
gene_list1 <- c()

for (i in 1:nrow(merged)){
  if (merged$log2FC_deseq2.x[i] > 0) {
    if (merged$log2FC_deseq2.y[i] < 0) {
      gene_list1 <- append(gene_list1, merged$gene_name[i])
    }
  } else {
    if (merged$log2FC_deseq2.y[i] > 0) {
      gene_list1 <- append(gene_list1, merged$gene_name[i])
    } 
  }
}

# Filter
gene_list1 <- as.data.frame(gene_list1)
merged1 <- merge(merged, gene_list1, by.x = "gene_name", by.y = "gene_list1")

# Plot heatmap for filtered DEGs
pheatmap::pheatmap(merged1[, c(2,5)], cluster_cols = F, main = "Heapmap of filtered DEGs in both cells", 
                   angle_col = 0, show_rownames = F, labels_col = c("Gad2","vGluT2"))

# Get heatmap genes order
out <- pheatmap::pheatmap(merged1[, c(2,5)], cluster_cols = F, main = "Heapmap of filtered DEGs in both cells", 
                          angle_col = 0, show_rownames = F, labels_col = c("Gad2","vGluT2"))

heatmap_genes <- merged1[rownames(merged1[out$tree_row[["order"]],]),]

#Enrichment analysis
# convert ensembl to entrez
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=DESeq2_Gad2$gene_name,
  mart=mart)
genes1 <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=DESeq2_vGluT2$gene_name,
  mart=mart)

# remove NA
genes <- na.omit(genes)
genes1 <- na.omit(genes1)

# remove vague translation
f <- plyr::count(genes, 2)
f <- f[f$freq == 1,]
genes <- merge(f, genes, by.x="entrezgene_id", by.y = "entrezgene_id")
genes <- merge(genes, DESeq2_Gad2, by.x = "ensembl_gene_id", by.y = "gene_name")

f <- plyr::count(genes1, 2)
f <- f[f$freq == 1,]
genes1 <- merge(f, genes1, by.x="entrezgene_id", by.y = "entrezgene_id")
genes1 <- merge(genes1, DESeq2_vGluT2, by.x = "ensembl_gene_id", by.y = "gene_name")

#Enrichment
library(clusterProfiler)
# Get background genes
uni <- read.delim(*path to any Slamdunk count file*, comment.char="#")
uni <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=uni$gene_name,
  mart=mart)
f <- plyr::count(uni, 2)
f <- f[f$freq == 1,]
uni <- f$entrezgene_id

# Anlysis
KEGG_Gad2 <- enrichKEGG(gene = genes$entrezgene_id, organism = "mmu", universe = as.character(uni))
MKEGG_Gad2 <- enrichMKEGG(gene = genes$entrezgene_id, organism = "mmu")
CC_Gad2 <- enrichGO(genes$entrezgene_id, OrgDb = "org.Mm.eg.db", ont = "CC")
MF_Gad2 <- enrichGO(genes$entrezgene_id, OrgDb = "org.Mm.eg.db", ont = "MF")
BP_Gad2 <- enrichGO(genes$entrezgene_id, OrgDb = "org.Mm.eg.db", ont = "BP")

KEGG_vGluT2 <- enrichKEGG(gene = genes1$entrezgene_id, organism = "mmu", universe = as.character(uni))
MKEGG_vGluT2 <- enrichMKEGG(gene = genes1$entrezgene_id, organism = "mmu")
CC_vGluT2 <- enrichGO(genes1$entrezgene_id, OrgDb = "org.Mm.eg.db", ont = "CC")
MF_vGluT2 <- enrichGO(genes1$entrezgene_id, OrgDb = "org.Mm.eg.db", ont = "MF")
BP_vGluT2 <- enrichGO(genes1$entrezgene_id, OrgDb = "org.Mm.eg.db", ont = "BP")

#Dot plot
dotplot(KEGG_Gad2, title = "Gad2 KEGG enrichment")
dotplot(MKEGG_Gad2, title = "Gad2 MKEGG enrichment")
dotplot(CC_Gad2, title = "Gad2 CC GO enrichment")
dotplot(MF_Gad2, title = "Gad2 MF GO enrichment")
dotplot(BP_Gad2, title = "Gad2 BP GO enrichment")

dotplot(KEGG_vGluT2, title = "vGluT2 KEGG enrichment")
dotplot(MKEGG_vGluT2, title = "vGluT2 MKEGG enrichment")
dotplot(CC_vGluT2, title = "vGluT2 GO enrichment")
dotplot(MF_vGluT2, title = "vGluT2 GO enrichment")
dotplot(BP_vGluT2, title = "vGluT2 GO enrichment")

# Enrichment map
emapplot(KEGG_Gad2)
emapplot(MKEGG_Gad2)
emapplot(CC_Gad2)
emapplot(MF_Gad2)
emapplot(BP_Gad2)

emapplot(KEGG_vGluT2)
emapplot(MKEGG_vGluT2)
emapplot(CC_vGluT2)
emapplot(MF_vGluT2)
emapplot(BP_vGluT2)

#Pathview
library(pathview)

#Gad2
logFC <- genes$log2FC_deseq2
names(logFC) <- genes$entrezgene_id
pathview(pathway.id = "mmu05022", gene.data = logFC, low = "red", mid = "grey", high = "green",
         species = "mmu",
         limit = c(-10, 10))

# vGluT2
logFC <- genes1$log2FC_deseq2
names(logFC) <- genes1$entrezgene_id
pathview(pathway.id = "mmu05022", gene.data = logFC, low = "red", mid = "grey", high = "green",
         species = "mmu",
         limit = c(-10, 10))


# MODifieR
# Create count matrix
RT_Gad2_1 <- read.delim(*path to Gad2 cell RiboTag condition 1st sample*, comment.char="#")
RT_Gad2_2 <- read.delim(*path to Gad2 cell RiboTag condition 2nd sample*, comment.char="#")
RT_Gad2_3 <- read.delim(*path to Gad2 cell RiboTag condition 3rd  sample*, comment.char="#")
S1_Gad2_1 <- read.delim(*path to Gad2 cell supernatant condition 1st sample*, comment.char="#")
S1_Gad2_2 <- read.delim(*path to Gad2 cell supernatant condition 2nd sample*, comment.char="#")
S1_Gad2_3 <- read.delim(*path to Gad2 cell supernatant condition 3rd sample*, comment.char="#")
Gad2 <- data.frame(RT_Gad2_1$gene_name ,RT_Gad2_1$tcReadCount,RT_Gad2_2$tcReadCount,RT_Gad2_3$tcReadCount,S1_Gad2_1$tcReadCount,S1_Gad2_2$tcReadCount,S1_Gad2_3$tcReadCount)
colnames(Gad2) <- c("ensembl_id","RT_Gad2_1", "RT_Gad2_2", "RT_Gad2_3", "S1_Gad2_1", "S1_Gad2_2", "S1_Gad2_3")

# Convert enseble id to entrezid
Gad2_genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=Gad2$ensembl_id,
  mart=mart)
f <- plyr::count(Gad2_genes, 2)
f <- f[f$freq == 1,]
Gad2_genes <- merge(f, Gad2_genes, by.x="entrezgene_id", by.y = "entrezgene_id")
Gad2_genes <- merge(Gad2_genes, Gad2, by.x = "ensembl_gene_id", by.y = "ensembl_id")
Gad2 <- Gad2_genes[,c(4,5,6,7,8,9)]
rownames(Gad2) <- Gad2_genes$entrezgene_id

RT_vGluT2_1 <- read.delim(*path to vGluT2 cell RiboTag condition 1st sample*, comment.char="#")
RT_vGluT2_2 <- read.delim(*path to vGluT2 cell RiboTag condition 2nd sample*, comment.char="#")
S1_vGluT2_1 <- read.delim(*path to vGluT2 cell supernatant condition 1st sample*, comment.char="#")
S1_vGluT2_2 <- read.delim(*path to vGluT2 cell supernatant condition 2nd sample*, comment.char="#")
S1_vGluT2_3 <- read.delim(*path to vGluT2 cell supernatant condition 3rd sample*, comment.char="#")
vGluT2 <- data.frame(RT_vGluT2_1$gene_name, RT_vGluT2_1$tcReadCount,RT_vGluT2_2$tcReadCount,S1_vGluT2_1$tcReadCount,S1_vGluT2_2$tcReadCount,S1_vGluT2_3$tcReadCount)
colnames(vGluT2) <- c("ensembl_id", "RT_vGluT2_1", "RT_vGluT2_2", "S1_vGluT2_1", "S1_vGluT2_2", "S1_vGluT2_3")

# Convert enseble id to entrezid
vGluT2_genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=vGluT2$ensembl_id,
  mart=mart)
f <- plyr::count(vGluT2_genes, 2)
f <- f[f$freq == 1,]
vGluT2_genes <- merge(f, vGluT2_genes, by.x="entrezgene_id", by.y = "entrezgene_id")
vGluT2_genes <- merge(vGluT2_genes, vGluT2, by.x = "ensembl_gene_id", by.y = "ensembl_id")
vGluT2 <- vGluT2_genes[,c(4,5,6,7,8)]
rownames(vGluT2) <- vGluT2_genes$entrezgene_id

# Create input
library(MODifieR)
Gad2_input <- create_input_rnaseq(count_matrix = as.matrix(Gad2), 
                    group1_indici = c(1,2,3), 
                    group2_indici = c(4,5,6), 
                    group1_label = "RT", 
                    group2_label = "S1", 
                    use_adjusted = F, 
                    normalize_quantiles = F)

vGluT2_input <- create_input_rnaseq(count_matrix = as.matrix(vGluT2), 
                                  group1_indici = c(1,2), 
                                  group2_indici = c(3,4,5), 
                                  group1_label = "RT", 
                                  group2_label = "S1", 
                                  use_adjusted = F, 
                                  normalize_quantiles = F)


# Add ppi network
ppi_network <- read.delim(*path to ppi network*, comment.char="#")
# Build database of cliques
build_clique_db(ppi_network = ppi_network,
                db_folder = *output folder*,
                db_name = *name for the db*)

# Analysis
#Gad2
MODA_GAD2 <- moda(MODifieR_input = Gad2_input, group_of_interest = 1)
MCODE_GAD2 <- mcode(MODifieR_input = Gad2_input, ppi_network = ppi_network)
WGCNA_Gad2 <- wgcna(MODifieR_input = Gad2_input, group_of_interest = 1)
DIFFCOEX_GAD2 <- diffcoex(MODifieR_input = Gad2_input)
CCLIQUE_GAD2 <- correlation_clique(MODifieR_input = Gad2_input, ppi_network = ppi_network)
MD_GAD2 <- modulediscoverer(MODifieR_input = Gad2_input, ppi_network = ppi_network, permutation = 1)
CS_GAD2 <- clique_sum_exact(Gad2_input, db = *path to the sqlite db*)

#vGluT2
MODA_vGluT2 <- moda(MODifieR_input = vGluT2_input, group_of_interest = 1)
MCODE_vGluT2 <- mcode(MODifieR_input = vGluT2_input, ppi_network = ppi_network)
WGCNA_vGluT2 <- wgcna(MODifieR_input = vGluT2_input, group_of_interest = 1)
DIFFCOEX_vGluT2 <- diffcoex(MODifieR_input = vGluT2_input)
CCLIQUE_vGluT2 <- correlation_clique(MODifieR_input = vGluT2_input, ppi_network = ppi_network)
MD_vGluT2 <- modulediscoverer(MODifieR_input = vGluT2_input, ppi_network = ppi_network, permutation = 1)
CS_vGluT2 <- clique_sum_exact(vGluT2_input, db = *path to the sqlite db*)

# Venn diagram
set1 <- CS_GAD2$module_genes
set2 <- WGCNA_Gad2$module_genes
set3 <- MCODE_GAD2$module_genes
library(VennDiagram)
draw.triple.venn(area1 = length(set1),
                 area2 = length(set2),
                 area3 = length(set3),
                 n12 = length(intersect(set1, set2)),
                 n23 = length(intersect(set2, set3)),
                 n13 = length(intersect(set1, set3)),
                 n123 = length(intersect(intersect(set1, set2), set3)),
                 category = c("Clique Sum", "WGCNA", "MCODE"),
                 lty = "blank",
                 fill = c("steelblue1", "yellowgreen", "indianred1") ,
                 cex = 2, cat.cex = 2, cat.fontfamily = rep("serif", 3))

# Enrichment analysis
Gad2_genes_mod <- intersect(intersect(set1, set2), set3)
KEGG_Gad2_mod <- enrichKEGG(gene = Gad2_genes_mod, organism = "mmu", universe = as.character(uni))
MKEGG_Gad2_mod <- enrichMKEGG(gene = Gad2_genes_mod, organism = "mmu")
CC_Gad2_mod <- enrichGO(Gad2_genes_mod, OrgDb = "org.Mm.eg.db", ont = "CC")
MF_Gad2_mod <- enrichGO(Gad2_genes_mod, OrgDb = "org.Mm.eg.db", ont = "MF")
BP_Gad2_mod <- enrichGO(Gad2_genes_mod, OrgDb = "org.Mm.eg.db", ont = "BP")

vGluT2_genes_mod <- intersect(intersect(set1, set2), set3)
KEGG_vGluT2_mod <- enrichKEGG(gene = vGluT2_genes_mod, organism = "mmu", universe = as.character(uni))
MKEGG_vGluT2_mod <- enrichMKEGG(gene = vGluT2_genes_mod, organism = "mmu")
CC_vGluT2_mod <- enrichGO(vGluT2_genes_mod, OrgDb = "org.Mm.eg.db", ont = "CC")
MF_vGluT2_mod <- enrichGO(vGluT2_genes_mod, OrgDb = "org.Mm.eg.db", ont = "MF")
BP_vGluT2_mod <- enrichGO(vGluT2_genes_mod, OrgDb = "org.Mm.eg.db", ont = "BP")

#Dot plot
dotplot(KEGG_Gad2_mod, title = "Gad2 KEGG enrichment for MODifieR gene module")
dotplot(MKEGG_Gad2_mod, title = "Gad2 MKEGG enrichment for MODifieR gene module")
dotplot(CC_Gad2_mod, title = "Gad2 CC GO enrichment for MODifieR gene module")
dotplot(MF_Gad2_mod, title = "Gad2 MF GO enrichment for MODifieR gene module")
dotplot(BP_Gad2_mod, title = "Gad2 BP GO enrichment for MODifieR gene module")

dotplot(KEGG_vGluT2_mod, title = "vGluT2 KEGG enrichment for MODifieR gene module")
dotplot(MKEGG_vGluT2_mod, title = "vGluT2 MKEGG enrichment for MODifieR gene module")
dotplot(CC_vGluT2_mod, title = "vGluT2 GO enrichment for MODifieR gene module")
dotplot(MF_vGluT2_mod, title = "vGluT2 GO enrichment for MODifieR gene module")
dotplot(BP_vGluT2_mod, title = "vGluT2 GO enrichment for MODifieR gene module")

# Combine dot plot
df <- data.frame(CC_Gad2_mod[1:5,])
df <- rbind(df, MF_Gad2_mod[1,])
df <- rbind(df, BP_Gad2_mod[1:5,])
Type <- c("CC","CC","CC","CC","CC","MF", "BP","BP","BP","BP","BP")
df <- cbind(df, Type)
library(ggplot2)
ggplot(df,
       aes(x = Type, y = sort(Description))) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_colour_gradient(limits=c(0, 0.10), low="red", high = "blue") +
  ylab(NULL) + xlab("Functional group") +
  ggtitle("GO pathway enrichment")

# Enrichment map
emapplot(KEGG_Gad2_mod)
emapplot(MKEGG_Gad2_mod)
emapplot(CC_Gad2_mod)
emapplot(MF_Gad2_mod)
emapplot(BP_Gad2_mod)

emapplot(KEGG_vGluT2_mod)
emapplot(MKEGG_vGluT2_mod)
emapplot(CC_vGluT2_mod)
emapplot(MF_vGluT2_mod)
emapplot(BP_vGluT2_mod)
