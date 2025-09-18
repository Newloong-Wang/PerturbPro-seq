library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(paletteer)
library(viridis)
library(ggridges)



### K562 seurat base analysis
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/seurat/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_after/crispr_obj.rds")

crispr_obj[["percent.mt"]] <- PercentageFeatureSet(crispr_obj, pattern = "^MT-")

VlnPlot(crispr_obj, features = "nFeature_RNA", pt.size = 0.1)
ggsave("violin_nFeature_RNA_qc_pre.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "nCount_RNA", pt.size = 0.1)
ggsave("violin_nCount_RNA_qc_pre.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "percent.mt", pt.size = 0.1)
ggsave("violin_percent_mt_qc_pre.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)


FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave("pct_counts_mt_qc_pre.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)

FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("n_genes_by_counts_qc_pre.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)

crispr_obj <- subset(crispr_obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 30)



##### Assign the sgRNA identity to each cell and store it in meta.data, then remove cells without an sgRNA identity
sgTypeCounts_df_filter2 <- read.csv("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_after/cell_sgRNA_identity.csv", row.names = 1,  check.names = FALSE, stringsAsFactors = FALSE)
sgRNA_values <- as.character(sgTypeCounts_df_filter2["sgRNA_identity", ])
barcodes <- colnames(sgTypeCounts_df_filter2)
sgRNA_identity_vector <- rep(NA, ncol(crispr_obj))
matching_indices <- match(barcodes, colnames(crispr_obj))
valid_indices <- matching_indices[!is.na(matching_indices)]
valid_sgRNA_values <- sgRNA_values[!is.na(matching_indices)]
sgRNA_identity_vector[valid_indices] <- valid_sgRNA_values
crispr_obj@meta.data$sgRNA_identity <- sgRNA_identity_vector

crispr_obj
na_cells <- is.na(crispr_obj@meta.data$sgRNA_identity)
crispr_obj <- subset(crispr_obj, cells = colnames(crispr_obj)[!na_cells])
crispr_obj

table(crispr_obj@meta.data[["sgRNA_identity"]])

VlnPlot(crispr_obj, features = "nFeature_RNA", pt.size = 0.1)
ggsave("violin_nFeature_RNA_qc.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "nCount_RNA", pt.size = 0.1)
ggsave("violin_nCount_RNA_qc.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "percent.mt", pt.size = 0.1)
ggsave("violin_percent_mt_qc.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)

FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave("pct_counts_mt_qc.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)

FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("n_genes_by_counts_qc.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)


crispr_obj <- NormalizeData(crispr_obj, 
                            normalization.method = "LogNormalize",
                            scale.factor = 1e4)
crispr_obj <- FindVariableFeatures(crispr_obj)

top10 <- head(VariableFeatures(crispr_obj), 10)
plot1 <- VariableFeaturePlot(crispr_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
combined_plot <- plot1 + plot2
ggsave("highly_variable_genes.pdf", plot = last_plot(), width = 12, height = 5, units = "in", dpi = 300)

all.genes <- rownames(crispr_obj)
crispr_obj <- ScaleData(crispr_obj, features = all.genes)

crispr_obj <- RunPCA(crispr_obj)
ElbowPlot(crispr_obj)
ggsave("ElbowPlot.pdf", plot = last_plot(), width = 4, height = 3, units = "in", dpi = 300)

crispr_obj <- FindNeighbors(crispr_obj, dims = 1:15)
crispr_obj <- RunUMAP(crispr_obj, dims = 1:15)


DimPlot(crispr_obj, reduction = "umap") + theme(legend.position = "bottom")
ggsave("umap.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)

crispr_obj$sgRNA_identity <- factor(crispr_obj$sgRNA_identity, 
                                    levels = c("NT1", "NT2", "AAVS", 
                                               "ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3",
                                               "RPS6-sg1", "RPS6-sg2", "RPS6-sg3",
                                               "GATA1-sg1", "GATA1-sg2", "GATA1-sg3"))
DimPlot(crispr_obj, reduction = "umap", group.by = "sgRNA_identity") + theme(legend.position = "bottom")
ggsave("umap_groupby_sgRNA.pdf", plot = last_plot(), width = 6, height = 7, units = "in", dpi = 300)
DimPlot(crispr_obj, reduction = "umap", split.by = "sgRNA_identity", ncol = 3) +
  theme(legend.position = "bottom")
ggsave("umap_splitby_sgRNA.pdf", plot = last_plot(), width = 10, height = 12, units = "in", dpi = 300)

saveRDS(crispr_obj, file = "crispr_obj_seurat.rds")



##### sgRNA
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/seurat/crispr_obj_seurat.rds")

table(crispr_obj@meta.data[["sgRNA_identity"]])


# Set the assay to sgRNA and perform normalization and scaling
DefaultAssay(crispr_obj) <- "sgRNA"
crispr_obj <- NormalizeData(crispr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(crispr_obj)
crispr_obj <- ScaleData(crispr_obj, features = all.genes)

crispr_obj <- RunPCA(crispr_obj, features = all.genes, reduction.name = "sgRNA_PCA", reduction.key = "sgRNAPCA_")
print(crispr_obj[["sgRNA_PCA"]], dims = 1:5, nfeatures = 10)  
ElbowPlot(crispr_obj, reduction = "sgRNA_PCA") 

crispr_obj <- FindNeighbors(crispr_obj, dims = 1:11, reduction = "sgRNA_PCA")  
crispr_obj <- FindClusters(crispr_obj, resolution = 0.5)  
table(crispr_obj@meta.data[["seurat_clusters"]])  

crispr_obj <- RunUMAP(crispr_obj, dims = 1:11, min.dist = 0.5, reduction = "sgRNA_PCA", reduction.name = "sgRNA_UMAP", reduction.key = "sgRNAUMAP_")

saveRDS(crispr_obj, file = "crispr_obj_seurat_sgRNA.rds")





##### Figure1
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/figure/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")
crispr_obj_seurat <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/seurat/crispr_obj_seurat.rds")
K562<- Read10X("/public/home/wangxinlong/Project/crispr1/data/20250102ABFC20221434-76_K562/sc20241231K562_addsequence/Ab-1_L2_S2401B2401_DBEC_MolsPerCell_MEX")
K562_seurat<-CreateSeuratObject(counts =K562$`Gene Expression`)
pal <-c("#E3A8A2",'#75A3D1','#E9AD95', "#65B8AD")
pal <-c('#75A3D1','#E3A8A2')
pal <-c( '#65B8AD','#DC7C56')



##### Figure1 d
#nCount_sgRNA
# # Set the assay to sgRNA and perform normalization and scaling
DefaultAssay(crispr_obj) <- "sgRNA"
crispr_obj <- NormalizeData(crispr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(crispr_obj)
crispr_obj <- ScaleData(crispr_obj, features = all.genes)
crispr_obj <- RunPCA(crispr_obj, features = all.genes, reduction.name = "sgRNA_PCA", reduction.key = "sgRNAPCA_")
print(crispr_obj[["sgRNA_PCA"]], dims = 1:5, nfeatures = 10)  
ElbowPlot(crispr_obj, reduction = "sgRNA_PCA") 

crispr_obj <- FindNeighbors(crispr_obj, dims = 1:11, reduction = "sgRNA_PCA")  
crispr_obj <- FindClusters(crispr_obj, resolution = 0.5)  
table(crispr_obj@meta.data[["seurat_clusters"]])  

crispr_obj <- RunUMAP(crispr_obj, dims = 1:11, min.dist = 0.5,reduction = "sgRNA_PCA", reduction.name = "sgRNA_UMAP", reduction.key = "sgRNAUMAP_")

mypal <- c("#AED6F1", '#6BB3E6', "#75A3D1", '#1E75B3', '#E3A8A2', '#CF5B5B', '#CC7177', '#B23E4D', '#51132F', '#BDE0B6', '#77BB88', '#55AA6A', "#65B8AD", '#4AA195', '#326C64', '#E69D70', '#DA712F', "#B65A20", '#F78426', "#E26A08", '#B15306', '#7F3C1A')

sgRNA_colors <- c(
  "ITGB1-sg1" = "#AED6F1",  
  "ITGB1-sg2" = '#6BB3E6',  
  "ITGB1-sg3" = "#4F8AC4",  
  "RPS6-sg1"  = "#EB9593",
  "RPS6-sg2"  = "#C26171",
  "RPS6-sg3"  = '#9A3C4C',  
  "GATA1-sg1" =  '#65B8AD',  
  "GATA1-sg2" = "#4AA195",  
  "GATA1-sg3" = "#326C64",  
  "NT1"      = "#FCB99C" ,
  "NT2" = "#FCB99C" ,
  "AAVS" ="#FCB99C" 
)

DimPlot(
  crispr_obj, 
  reduction = "sgRNA_UMAP", 
  group.by = "sgRNA_identity", 
  cols = sgRNA_colors,
  label = TRUE, 
  repel = TRUE
) + 
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5)
  ) 
ggsave("sgRNA_UMAP.pdf", plot = last_plot(), width = 6.5, height = 5, units = "in", dpi = 300)



##### Figure1 e
# No threshold applied ----
VlnPlot(K562_seurat,features=c("nFeature_RNA"),ncol=1, pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")+ 
  scale_fill_manual(values ='#75A3D1' ) + 
  NoLegend()
median(K562_seurat$nFeature_RNA)
ggsave("Nothresholdapplied_allCell_3942.nFeature_RNA.pdf", plot = last_plot(), width = 4, height = 5, units = "in", dpi = 300)



##### Figure1 g
# Set the protein as the default assay and normalize it
DefaultAssay(crispr_obj_seurat) <- "Protein"
crispr_obj_seurat <- NormalizeData(crispr_obj_seurat, normalization.method = "CLR", margin = 2)

# CD29_protein
cd29_protein_data <- FetchData(crispr_obj_seurat, vars = c("CD29-pAbO", "sgRNA_identity"))
cd29_protein_data$sgRNA_identity <- as.character(cd29_protein_data$sgRNA_identity)
cd29_protein_data$sgRNA_identity[cd29_protein_data$sgRNA_identity %in% c("NT1", "NT2", "AAVS")] <- "NT"
cd29_protein_data <- cd29_protein_data %>% filter(sgRNA_identity %in% c("ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3", "NT"))
cd29_protein_data$sgRNA_identity <- factor(cd29_protein_data$sgRNA_identity, levels = c("ITGB1-sg3", "ITGB1-sg2", "ITGB1-sg1", "NT"))

ggplot(cd29_protein_data, aes(x = `CD29-pAbO`, y = sgRNA_identity, fill = sgRNA_identity, color = sgRNA_identity)) +
  geom_density_ridges(scale = 1.5, alpha = 0.7, bandwidth = 0.2, rel_min_height = 0, quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6")) +
  scale_color_manual(values = darken(c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6"), amount = 0.3)) +
  labs(x = "CD29 Protein Expression", y = "sgRNA Identity") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray", size = 0.5),
        panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_line(color = "black", size = 0.5)) +
  scale_y_discrete(expand = expansion(mult = c(0, 0.5)))
ggsave("density_Protein_expression_CD29.pdf", width = 5.5, height = 5, units = "in", dpi = 300)

# ITGB1_RNA
cd29_rna_data <- FetchData(crispr_obj_seurat, vars = c("ITGB1", "sgRNA_identity"))
cd29_rna_data$sgRNA_identity <- as.character(cd29_rna_data$sgRNA_identity)
cd29_rna_data$sgRNA_identity[cd29_rna_data$sgRNA_identity %in% c("NT1", "NT2", "AAVS")] <- "NT"
cd29_rna_data <- cd29_rna_data %>% filter(sgRNA_identity %in% c("ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3", "NT"))
cd29_rna_data$sgRNA_identity <- factor(cd29_rna_data$sgRNA_identity, levels = c("ITGB1-sg3", "ITGB1-sg2", "ITGB1-sg1", "NT"))

ggplot(cd29_rna_data, aes(x = rna_ITGB1, y = sgRNA_identity, fill = sgRNA_identity, color = sgRNA_identity)) +
  geom_density_ridges(scale = 1.5, alpha = 0.7, bandwidth = 0.4, rel_min_height = 0, quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6")) +
  scale_color_manual(values = darken(c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6"), amount = 0.3)) +
  labs(x = "ITGB1 RNA Expression", y = "sgRNA Identity") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray", size = 0.5),
        panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_line(color = "black", size = 0.5)) +
  scale_y_discrete(expand = expansion(mult = c(0, 0.51)))
ggsave("density_RNA_expression_ITGB1.pdf", width = 5.5, height = 5, units = "in", dpi = 300)



##### Figure1 h
# GATA1_protein
gata1_protein_data <- FetchData(crispr_obj_seurat, vars = c("GATA1-pAbO", "sgRNA_identity"))
gata1_protein_data$sgRNA_identity <- as.character(gata1_protein_data$sgRNA_identity)
gata1_protein_data$sgRNA_identity[gata1_protein_data$sgRNA_identity %in% c("NT1", "NT2", "AAVS")] <- "NT"
gata1_protein_data <- gata1_protein_data %>% filter(sgRNA_identity %in% c("GATA1-sg1", "GATA1-sg2", "GATA1-sg3", "NT"))
gata1_protein_data$sgRNA_identity <- factor(gata1_protein_data$sgRNA_identity, levels = c("GATA1-sg3", "GATA1-sg2", "GATA1-sg1", "NT"))

ggplot(gata1_protein_data, aes(x = `GATA1-pAbO`, y = sgRNA_identity, fill = sgRNA_identity, color = sgRNA_identity)) +
  geom_density_ridges(scale = 1.5, alpha = 0.7, bandwidth = 0.3, rel_min_height = 0, quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6")) +
  scale_color_manual(values = darken(c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6"), amount = 0.3)) +
  labs(x = "GATA1 Protein Expression", y = "sgRNA Identity") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray", size = 0.5),
        panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_line(color = "black", size = 0.5)) +
  scale_y_discrete(expand = expansion(mult = c(0, 0.5)))
ggsave("density_Protein_expression_GATA1.pdf", width = 5.5, height = 5, units = "in", dpi = 300)

# GATA1_RNA
gata1_rna_data <- FetchData(crispr_obj_seurat, vars = c("GATA1", "sgRNA_identity"))
gata1_rna_data$sgRNA_identity <- as.character(gata1_rna_data$sgRNA_identity)
gata1_rna_data$sgRNA_identity[gata1_rna_data$sgRNA_identity %in% c("NT1", "NT2", "AAVS")] <- "NT"
gata1_rna_data <- gata1_rna_data %>% filter(sgRNA_identity %in% c("GATA1-sg1", "GATA1-sg2", "GATA1-sg3", "NT"))
gata1_rna_data$sgRNA_identity <- factor(gata1_rna_data$sgRNA_identity, levels = c("GATA1-sg3", "GATA1-sg2", "GATA1-sg1", "NT"))

ggplot(gata1_rna_data, aes(x = rna_GATA1, y = sgRNA_identity, fill = sgRNA_identity, color = sgRNA_identity)) +
  geom_density_ridges(scale = 1.5, alpha = 0.7, bandwidth = 0.4, rel_min_height = 0, quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_manual(values = c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6")) +
  scale_color_manual(values = darken(c("#F8CCCC", "#CEE3F6", "#FDE5CE", "#D6D6D6"), amount = 0.3)) +
  labs(x = "GATA1 RNA Expression", y = "sgRNA Identity") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray", size = 0.5),
        panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.ticks.x = element_line(color = "black", size = 0.5)) +
  scale_y_discrete(expand = expansion(mult = c(0, 0.51)))
ggsave("density_RNA_expression_GATA1.pdf", width = 5.5, height = 5, units = "in", dpi = 300)



##### Figure1 g
# CD29 mRNA Protein ----
CD29_target_sgRNA <- c("NT1",'NT2','AAVS', "ITGB1-sg1", "ITGB1-sg2",'ITGB1-sg3')
CD29_crispr_subset <- subset(crispr_obj_seurat, subset = sgRNA_identity %in% CD29_target_sgRNA)
CD29_crispr_subset$sgRNA_identity <- as.character(CD29_crispr_subset$sgRNA_identity)

# Merge NT1, NT2, and AAVS into “NT”
CD29_crispr_subset$sgRNA_identity[
  CD29_crispr_subset$sgRNA_identity %in% c("NT1", "NT2", "AAVS")
] <- "NT"

CD29_expression_data <- FetchData(
  CD29_crispr_subset,
  vars = c("ITGB1", "CD29-pAbO", "sgRNA_identity"),
  assays = c("RNA", "ADT"),  
  slot = "data"              
)
colnames(CD29_expression_data) <- c("ITGB1_mRNA", "CD29_protein", "sgRNA_identity")

CD29_expression_long <- CD29_expression_data %>%
  pivot_longer(
    cols = c(ITGB1_mRNA, CD29_protein),
    names_to = "measurement",
    values_to = "value"
  ) %>%
  mutate(
    sgRNA_identity = factor(
      sgRNA_identity,
      levels = c("NT", "ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3")
    ),
    measurement = factor(
      measurement,
      levels = c("ITGB1_mRNA", "CD29_protein"),
      labels = c("ITGB1 mRNA Expression", "CD29-pAbO Protein Expression")
    )
  )

comparisons <- list(
  c("NT", "ITGB1-sg1"),
  c("NT", "ITGB1-sg2"),
  c("NT", "ITGB1-sg3")
)


fill_colors <- c("#fbac91", "#AED6F1", '#6BB3E6', "#75A3D1")
border_colors <- darken(fill_colors, amount = 0.3)

ggplot(CD29_expression_long, aes(x = sgRNA_identity, y = value, fill = sgRNA_identity)) + 
  geom_violin(aes(color = sgRNA_identity), alpha = 0.8, width = 1, trim = TRUE) + 
  geom_boxplot(aes(color = sgRNA_identity), width = 0.2, fill = "white", outlier.shape = NA, alpha = 0.5) + 
  facet_wrap(
    ~ measurement, 
    scales = "free_y", 
    strip.position = "left", 
    nrow = 1
  ) + 
  
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", 
    label = "p.signif", 
    step.increase = 0.15,  
    size = 4, 
    bracket.size = 0.5,    
    tip.length = 0.01,
    vjust = 0.5, 
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) + 
  scale_fill_manual(values = fill_colors) + 
  scale_color_manual(values = border_colors) +
  labs(
    x = "sgRNA Group", 
    y = NULL,
    title = "CD29 mRNA and CD29 Protein Expression by sgRNA Group"
  ) + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )
ggsave("vln_CD29_expression.pdf", plot = last_plot(), width = 10, height = 7, units = "in", dpi = 300)



##### Figure1 h
# GATA1 mRNA Protein ----
target_sgRNA <- c("NT1",'NT2','AAVS', "GATA1-sg1", "GATA1-sg2",'GATA1-sg3')
crispr_subset <- subset(crispr_obj_seurat, subset = sgRNA_identity %in% target_sgRNA)
crispr_subset$sgRNA_identity <- as.character(crispr_subset$sgRNA_identity)

# Merge NT1, NT2, and AAVS into “NT”
crispr_subset$sgRNA_identity[
  crispr_subset$sgRNA_identity %in% c("NT1", "NT2", "AAVS")
] <- "NT"

expression_data <- FetchData(
  crispr_subset,
  vars = c("GATA1", "GATA1-pAbO", "sgRNA_identity"),
  assays = c("RNA", "ADT"), 
  slot = "data"
)
colnames(expression_data) <- c("GATA1_mRNA", "GATA1_protein", "sgRNA_identity")

expression_long <- expression_data %>%
  pivot_longer(
    cols = c(GATA1_mRNA, GATA1_protein),
    names_to = "measurement",
    values_to = "value"
  ) %>%
  mutate(
    sgRNA_identity = factor(
      sgRNA_identity,
      levels = c("NT", "GATA1-sg1", "GATA1-sg2", "GATA1-sg3")
    ),
    measurement = factor(
      measurement,
      levels = c("GATA1_mRNA", "GATA1_protein"),
      labels = c("GATA1 mRNA Expression", "GATA1-pAbO Protein Expression")
    )
  )

comparisons <- list(
  c("NT", "GATA1-sg1"),
  c("NT", "GATA1-sg2"),
  c("NT", "GATA1-sg3"))


fill_colors <- c("#fbac91", "#88C8C1", "#4AA195", "#326C64")
border_colors <- darken(fill_colors, amount = 0.3)

ggplot(expression_long, aes(x = sgRNA_identity, y = value, fill = sgRNA_identity)) + 
  geom_violin(aes(color = sgRNA_identity), alpha = 0.7, width = 1, trim = TRUE) + 
  geom_boxplot(aes(color = sgRNA_identity), width = 0.2, fill = "white", outlier.shape = NA, alpha = 0.5) + 
  facet_wrap(
    ~ measurement, 
    scales = "free_y", 
    strip.position = "left", 
    nrow = 1
  ) + 
  stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", 
    label = "p.signif", 
    step.increase = 0.15,
    size = 4, 
    bracket.size = 0.5,
    tip.length = 0.01,
    vjust = 0.5, 
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) + 
  scale_fill_manual(values = fill_colors) + 
  scale_color_manual(values = border_colors) +
  labs(
    x = "sgRNA Group", 
    y = NULL,
    title = "GATA1 mRNA and GATA1 Protein Expression by sgRNA Group"
  ) + 
  theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )
ggsave("vln_GATA1_expression.pdf", plot = last_plot(), width = 10, height = 7, units = "in", dpi = 300)





##### Figure1 f
library(dplyr)
library(ggplot2)

setwd("/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/figure_2/")

file_path <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/merge/combined_merged_all.tsv"

combined_merged <- read.delim(file_path, stringsAsFactors = FALSE)

combined_merged <- combined_merged[combined_merged$group != "GSM4367984", ]
combined_merged$group <- factor(combined_merged$group, levels = c("Chromium", "ECCITE", "PerturbPro", "Perturb"))
combined_merged <- combined_merged %>% arrange(group)

custom_colors <- c(
  "Chromium" = "#F4A88E",
  "PerturbPro" = "#5E89BA",
  "Perturb" = "#9594D1",
  "ECCITE" = "#85C1BB"
)

# boxplot genes
ggplot(combined_merged, aes(x = group, y = nFeature_RNA, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 1) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 0)), 
    vjust = -0.5,                 
    size = 4,
    color = "black"
  ) +
  labs(
    x = NULL,
    y = "Number of genes"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave("boxplot_genes.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)





##### Figure1 i
### GATA1_protein_density
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GATA1_protein_density/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GATA1_protein_density/")


# Load object.
crispr_obj <- readRDS("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")

# Save the knocked-out gene to the gene column in meta.data
sgRNA_identity <- crispr_obj@meta.data$sgRNA_identity
gene <- rep(NA, length(sgRNA_identity))
gene[sgRNA_identity == "NT1"] <- "NT"
gene[sgRNA_identity == "NT2"] <- "NT"
gene[sgRNA_identity == "AAVS"] <- "NT"
gene[sgRNA_identity %in% c("ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3")] <- "ITGB1"
gene[sgRNA_identity %in% c("RPS6-sg1", "RPS6-sg2", "RPS6-sg3")] <- "RPS6"
gene[sgRNA_identity %in% c("GATA1-sg1", "GATA1-sg2", "GATA1-sg3")] <- "GATA1"
gene <- factor(gene, levels = c("NT", "ITGB1", "RPS6", "GATA1"))
crispr_obj@meta.data$gene <- gene

table(crispr_obj@meta.data$gene)


DefaultAssay(crispr_obj) <- "Protein"

subset_obj <- subset(crispr_obj, subset = gene %in% c("GATA1"))

### KDE
gata1_expression <- FetchData(subset_obj, vars = "GATA1-pAbO")$`GATA1-pAbO`

q20 <- quantile(gata1_expression, 0.2)
q30 <- quantile(gata1_expression, 0.3)
q70 <- quantile(gata1_expression, 0.7)
q80 <- quantile(gata1_expression, 0.8)

dens <- density(gata1_expression, adjust = 2, n = 1024)
dens_df <- data.frame(x = dens$x, y = dens$y)

dens_df <- dens_df %>%
  mutate(region = case_when(
    x <= q30 ~ "Bottom 30%",
    x >= q70 ~ "Top 30%",
    TRUE ~ "Middle 40%"
  ))

y_max <- max(dens_df$y)

ggplot(dens_df, aes(x = x, y = y, fill = region)) +
  geom_area(data = subset(dens_df, region == "Bottom 30%"), fill = "#A4D2E0", alpha = 0.4) +
  geom_area(data = subset(dens_df, region == "Middle 40%"), fill = "#EEEEEE", alpha = 0.5) +
  geom_area(data = subset(dens_df, region == "Top 30%"), fill = "#DF9796", alpha = 0.4) +
  geom_line(color = "black", size = 0.8) +
  geom_vline(xintercept = q20, linetype = "dashed", color = "#A4D2E0", size = 1) +
  geom_vline(xintercept = q80, linetype = "dashed", color = "#DF9796", size = 1) +
  geom_vline(xintercept = q30, linetype = "dashed", color = "#A4D2E0", size = 1) +
  geom_vline(xintercept = q70, linetype = "dashed", color = "#DF9796", size = 1) +
  annotate("text", x = q20, y = y_max * 0.9, label = "Bottom 20%", color = "#A4D2E0", hjust = 1.1, size = 4) +
  annotate("text", x = q80, y = y_max * 0.9, label = "Top 20%", color = "#DF9796", hjust = -0.1, size = 4) +
  annotate("text", x = q30, y = y_max * 0.3, label = "Bottom 30%", color = "#A4D2E0", hjust = -0.05, size = 4) +
  annotate("text", x = q70, y = y_max * 0.2, label = "Top 30%", color = "#DF9796", hjust = 1.1, size = 4) +
  labs(
    title = "GATA1 protein expression distribution",
    x = "Expression Level",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("KDE_GATA1_protein.pdf", plot = last_plot(), width = 10, height = 7, units = "in", dpi = 300)





##### Figure1 j
##### mixscale GATA1 3 group(all GATA1 cell / remove top_20 / remove top_30)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(mixtools)
library(Mixscale)
library(ggpubr)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/mixscale_score_3group/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/mixscale_score_3group/")


# Load object.
crispr_obj <- readRDS("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")

# Save the knocked-out gene to the gene column in meta.data
sgRNA_identity <- crispr_obj@meta.data$sgRNA_identity
gene <- rep(NA, length(sgRNA_identity))
gene[sgRNA_identity == "NT1"] <- "NT"
gene[sgRNA_identity == "NT2"] <- "NT"
gene[sgRNA_identity == "AAVS"] <- "NT"
gene[sgRNA_identity %in% c("ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3")] <- "ITGB1"
gene[sgRNA_identity %in% c("RPS6-sg1", "RPS6-sg2", "RPS6-sg3")] <- "RPS6"
gene[sgRNA_identity %in% c("GATA1-sg1", "GATA1-sg2", "GATA1-sg3")] <- "GATA1"
gene <- factor(gene, levels = c("NT", "ITGB1", "RPS6", "GATA1"))
crispr_obj@meta.data$gene <- gene

table(crispr_obj@meta.data$gene)


DefaultAssay(crispr_obj) <- "Protein"

subset_obj <- subset(crispr_obj, subset = gene %in% c("NT", "GATA1"))


# Divide GATA1 sgRNA cells based on GATA1 protein expression levels
gata1_expression <- FetchData(subset_obj, vars = c("GATA1-pAbO", "gene")) %>%
  dplyr::filter(gene == "GATA1")
sorted_cells <- gata1_expression[order(gata1_expression$`GATA1-pAbO`, decreasing = TRUE), , drop = FALSE]

num_cells <- nrow(sorted_cells)
top_30 <- ceiling(0.3 * num_cells)
bottom_30 <- floor(0.3 * num_cells)
subset_obj$protein_30 <- "middle"
top_cells <- rownames(sorted_cells)[1:top_30]
subset_obj$protein_30[top_cells] <- "top_30"
bottom_cells <- rownames(sorted_cells)[(num_cells - bottom_30 + 1):num_cells]
subset_obj$protein_30[bottom_cells] <- "bottom_30"

top_20 <- ceiling(0.2 * num_cells)
top20_cells <- rownames(sorted_cells)[1:top_20]
subset_obj$protein_30[top20_cells] <- "top_20"

subset_obj$protein_30[subset_obj$gene == "NT"] <- "NT"


### all GATA1 cell
subset_obj
DefaultAssay(subset_obj) <- "RNA"

# standard pre-processing
subset_obj = NormalizeData(subset_obj)
subset_obj = FindVariableFeatures(subset_obj)
subset_obj = ScaleData(subset_obj)
subset_obj = RunPCA(subset_obj)

# calculate Perturbation signatures 
subset_obj <- CalcPerturbSig(
  object = subset_obj, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB", 
  split.by = NULL)

# Mixscale
subset_obj = RunMixscale(
  object = subset_obj, 
  assay = "PRTB", 
  slot = "data", 
  labels = "gene", 
  nt.class.name = "NT", 
  min.de.genes = 5, 
  logfc.threshold = 0.25,
  de.assay = "RNA",
  max.de.genes = 100, 
  new.class.name = "mixscale_score", 
  fine.mode = F, 
  verbose = F, 
  split.by = NULL)

meta_df <- subset_obj@meta.data
df_all <- meta_df %>%
  filter(protein_30 %in% c("top_20", "top_30", "middle", "bottom_30")) %>%
  mutate(Group = "All")


### remove top_20
subset_obj <- subset(subset_obj, subset = protein_30 != "top_20")

subset_obj
DefaultAssay(subset_obj) <- "RNA"

# standard pre-processing
subset_obj = NormalizeData(subset_obj)
subset_obj = FindVariableFeatures(subset_obj)
subset_obj = ScaleData(subset_obj)
subset_obj = RunPCA(subset_obj)

# calculate Perturbation signatures 
subset_obj <- CalcPerturbSig(
  object = subset_obj, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB", 
  split.by = NULL)

# Mixscale
subset_obj = RunMixscale(
  object = subset_obj, 
  assay = "PRTB", 
  slot = "data", 
  labels = "gene", 
  nt.class.name = "NT", 
  min.de.genes = 5, 
  logfc.threshold = 0.25,
  de.assay = "RNA",
  max.de.genes = 100, 
  new.class.name = "mixscale_score", 
  fine.mode = F, 
  verbose = F, 
  split.by = NULL)

meta_df <- subset_obj@meta.data
df_remove_top20 <- meta_df %>%
  filter(protein_30 %in% c("top_30", "middle", "bottom_30")) %>%
  mutate(Group = "remove_top20")


### remove top_30
subset_obj <- subset(subset_obj, subset = protein_30 != "top_30")

subset_obj
DefaultAssay(subset_obj) <- "RNA"

# standard pre-processing
subset_obj = NormalizeData(subset_obj)
subset_obj = FindVariableFeatures(subset_obj)
subset_obj = ScaleData(subset_obj)
subset_obj = RunPCA(subset_obj)

# calculate Perturbation signatures 
subset_obj <- CalcPerturbSig(
  object = subset_obj, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB", 
  split.by = NULL)

# Mixscale
subset_obj = RunMixscale(
  object = subset_obj, 
  assay = "PRTB", 
  slot = "data", 
  labels = "gene", 
  nt.class.name = "NT", 
  min.de.genes = 5, 
  logfc.threshold = 0.25,
  de.assay = "RNA",
  max.de.genes = 100, 
  new.class.name = "mixscale_score", 
  fine.mode = F, 
  verbose = F, 
  split.by = NULL)

meta_df <- subset_obj@meta.data
df_remove_top30 <- meta_df %>%
  filter(protein_30 %in% c("middle", "bottom_30")) %>%
  mutate(Group = "remove_top30")


# Merge the three data
plot_df <- bind_rows(df_all, df_remove_top20, df_remove_top30) %>%
  mutate(highlight = case_when(
    Group == "All" & protein_30 == "top_20" ~ "highlight_top20",
    Group == "All" & protein_30 == "top_30" ~ "highlight_top30",
    Group == "remove_top20" & protein_30 == "top_30" ~ "highlight",
    TRUE ~ "normal"
  ))

plot_df$Group <- factor(plot_df$Group, levels = c("All", "remove_top20", "remove_top30"))

comparisons <- list(
  c("All", "remove_top20"),
  c("All", "remove_top30"),
  c("remove_top20", "remove_top30")
)

ggplot(plot_df, aes(x = Group, y = mixscale_score, fill = Group)) +
  geom_violin(width = 0.9, color = NA, alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.6) +
  geom_jitter(data = filter(plot_df, highlight == "normal"),
              color = "#F0F0F0", width = 0.15, size = 1) +
  geom_jitter(data = filter(plot_df, highlight == "highlight_top20"),
              color = "#C2574A", width = 0.15, size = 1.5) +
  geom_jitter(data = filter(plot_df, highlight == "highlight_top30"),
              color = "#80539F", width = 0.15, size = 1.5) +
  geom_jitter(data = filter(plot_df, highlight == "highlight"),
              color = "#80539F", width = 0.15, size = 1.5) +
  stat_summary(fun = median, geom = "text",
               aes(label = round(..y.., 2)),
               vjust = -0.6, size = 4, fontface = "bold", color = "black") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c(
    "All" = "#8DB7D8",
    "remove_top20" = "#679ACD",
    "remove_top30" = "#2F5EAB"
  )) +
  labs(y = "Mixscale score") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none"
  )
ggsave(filename = "mixscale_score_3group.pdf", plot = last_plot(), width = 8, height = 6, dpi = 300)





##### Figure1 k
##### mixscale all cell
library(Seurat)
library(ggridges)
library(ggplot2)
library(Mixscale)
library(stringr)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/")

### load data
seurat_obj <- readRDS("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")

DefaultAssay(seurat_obj) <- "RNA"


# Save the knocked-out gene to the gene column in meta.data
sgRNA_identity <- seurat_obj@meta.data$sgRNA_identity
gene <- rep(NA, length(sgRNA_identity))
gene[sgRNA_identity == "NT1"] <- "NT"
gene[sgRNA_identity == "NT2"] <- "NT"
gene[sgRNA_identity == "AAVS"] <- "NT"
gene[sgRNA_identity %in% c("ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3")] <- "ITGB1"
gene[sgRNA_identity %in% c("RPS6-sg1", "RPS6-sg2", "RPS6-sg3")] <- "RPS6"
gene[sgRNA_identity %in% c("GATA1-sg1", "GATA1-sg2", "GATA1-sg3")] <- "GATA1"
gene <- factor(gene, levels = c("NT", "ITGB1", "RPS6", "GATA1"))
seurat_obj@meta.data$gene <- gene

table(seurat_obj@meta.data$gene)


### Pre-processing and calculating the Mixscale score
seurat_obj

# calculate Perturbation signatures 
seurat_obj <- CalcPerturbSig(
  object = seurat_obj, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB", 
  split.by = NULL)

# Mixscale
seurat_obj = RunMixscale(
  object = seurat_obj, 
  assay = "PRTB", 
  slot = "data", 
  labels = "gene", 
  nt.class.name = "NT", 
  min.de.genes = 5, 
  logfc.threshold = 0.25,
  de.assay = "RNA",
  max.de.genes = 100, 
  new.class.name = "mixscale_score", 
  fine.mode = F, 
  verbose = F, 
  split.by = NULL)


# mixscalescore_proteinGATA1expression
cells_with_GATA1 <- subset(seurat_obj, subset = gene == "GATA1")
GATA1_df <- FetchData(cells_with_GATA1, vars = c("GATA1-pAbO", "mixscale_score"))
GATA1_df <- GATA1_df[order(GATA1_df$mixscale_score), ]
colnames(GATA1_df)[1] <- "protein_GATA1_pAbO"

cor_res <- cor.test(GATA1_df$mixscale_score, GATA1_df$protein_GATA1_pAbO, method = "spearman")
rho <- round(cor_res$estimate, 3)
p_val <- signif(cor_res$p.value, 3)

ggplot(GATA1_df, aes(x = mixscale_score, y = protein_GATA1_pAbO)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  labs(x = "mixscale_score", y = "Protein_GATA1") +
  theme_bw()
ggsave(filename = "mixscalescore_proteinGATA1expression.pdf", plot = last_plot(), width = 8, height = 6, dpi = 300)


### Differential expression (DE) analysis
de_res = Run_wmvRegDE(object = seurat_obj, assay = "RNA", slot = "counts",
                      labels = "gene", nt.class.name = "NT",
                      logfc.threshold = 0,
                      split.by = NULL,
                      full.results = F)

# Since the software did not strictly filter by logfc.threshold, apply logFC filtering manually
logfc_threshold <- 0
de_res <- lapply(de_res, function(df) {
  df[df$log2FC >= logfc_threshold | df$log2FC <= -logfc_threshold, ]
})

for (gene in names(de_res)) {
  file_name <- paste0(gene, "_DE_results_logfc0.tsv")
  write.table(de_res[[gene]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


# Since the software did not strictly filter by logfc.threshold, apply logFC filtering manually
logfc_threshold <- 0.5
de_res <- lapply(de_res, function(df) {
  df[df$log2FC >= logfc_threshold | df$log2FC <= -logfc_threshold, ]
})

for (gene in names(de_res)) {
  file_name <- paste0(gene, "_DE_results_logfc05.tsv")
  write.table(de_res[[gene]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


### GATA1_DE_results_logfc0  Calculate FDR
file_path <- "GATA1_DE_results_logfc0.tsv"
deg <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deg$FDR <- p.adjust(deg$p_weight, method = "BH")
write.table(deg, file = "GATA1_DE_results_logfc0_padjust.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


seurat_obj@meta.data$gene <- as.character(seurat_obj@meta.data$gene)



##### Divide GATA1 sgRNA cells into top 20, middle, and bottom 20 based on GATA1 protein expression levels
##### mixscale_GATA1_top20middlebottoom20
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(mixtools)
library(Mixscale)
library(ggpubr)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/")


# Load object.
crispr_obj <- readRDS("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")

# Save the knocked-out gene to the gene column in meta.data
sgRNA_identity <- crispr_obj@meta.data$sgRNA_identity
gene <- rep(NA, length(sgRNA_identity))
gene[sgRNA_identity == "NT1"] <- "NT"
gene[sgRNA_identity == "NT2"] <- "NT"
gene[sgRNA_identity == "AAVS"] <- "NT"
gene[sgRNA_identity %in% c("ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3")] <- "ITGB1"
gene[sgRNA_identity %in% c("RPS6-sg1", "RPS6-sg2", "RPS6-sg3")] <- "RPS6"
gene[sgRNA_identity %in% c("GATA1-sg1", "GATA1-sg2", "GATA1-sg3")] <- "GATA1"
gene <- factor(gene, levels = c("NT", "ITGB1", "RPS6", "GATA1"))
crispr_obj@meta.data$gene <- gene

table(crispr_obj@meta.data$gene)


DefaultAssay(crispr_obj) <- "Protein"

subset_obj <- subset(crispr_obj, subset = gene %in% c("NT", "GATA1"))


# Divide GATA1 sgRNA cells into top 20, middle, and bottom 20 based on GATA1 protein expression levels
gata1_expression <- FetchData(subset_obj, vars = c("GATA1-pAbO", "gene")) %>%
  dplyr::filter(gene == "GATA1")
sorted_cells <- gata1_expression[order(gata1_expression$`GATA1-pAbO`, decreasing = TRUE), , drop = FALSE]
num_cells <- nrow(sorted_cells)
top_20 <- ceiling(0.2 * num_cells)
bottom_20 <- floor(0.2 * num_cells)
subset_obj$protein_20 <- "middle"
top_cells <- rownames(sorted_cells)[1:top_20]
subset_obj$protein_20[top_cells] <- "top_20"
bottom_cells <- rownames(sorted_cells)[(num_cells - bottom_20 + 1):num_cells]
subset_obj$protein_20[bottom_cells] <- "bottom_20"
subset_obj$protein_20[subset_obj$gene == "NT"] <- "NT"


# Delete cells where protein_20 is “top_20”
subset_obj <- subset(subset_obj, subset = protein_20 != "top_20")

# standard pre-processing
DefaultAssay(subset_obj) <- "RNA"
subset_obj = NormalizeData(subset_obj)
subset_obj = FindVariableFeatures(subset_obj)
subset_obj = ScaleData(subset_obj)
subset_obj = RunPCA(subset_obj)

# calculate Perturbation signatures 
subset_obj <- CalcPerturbSig(
  object = subset_obj, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB", 
  split.by = NULL)

# Mixscale
subset_obj = RunMixscale(
  object = subset_obj, 
  assay = "PRTB", 
  slot = "data", 
  labels = "gene", 
  nt.class.name = "NT", 
  min.de.genes = 5, 
  logfc.threshold = 0.25,
  de.assay = "RNA",
  max.de.genes = 100, 
  new.class.name = "mixscale_score", 
  fine.mode = F, 
  verbose = F, 
  split.by = NULL)


### Differential expression (DE) analysis
de_res = Run_wmvRegDE(object = subset_obj, assay = "RNA", slot = "counts",
                      labels = "gene", nt.class.name = "NT",
                      logfc.threshold = 0,
                      split.by = NULL,
                      full.results = F)

# Since the software did not strictly filter by logfc.threshold, apply logFC filtering manually
logfc_threshold <- 0
de_res <- lapply(de_res, function(df) {
  df[df$log2FC >= logfc_threshold | df$log2FC <= -logfc_threshold, ]
})

for (gene in names(de_res)) {
  file_name <- paste0(gene, "_DE_results_logfc0.tsv")
  write.table(de_res[[gene]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


### GATA1_DE_results_logfc0 Calculate FDR
file_path <- "GATA1_DE_results_logfc0.tsv"
deg <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deg$FDR <- p.adjust(deg$p_weight, method = "BH")
write.table(deg, file = "GATA1_DE_results_logfc0_padjust.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



##### Divide GATA1 sgRNA cells into top 30, middle, and bottom 30 based on GATA1 protein expression levels
##### mixscale_GATA1_top30middlebottoom30
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(mixtools)
library(Mixscale)
library(ggpubr)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/")


# Load object.
crispr_obj <- readRDS("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")

# Save the knocked-out gene to the gene column in meta.data
sgRNA_identity <- crispr_obj@meta.data$sgRNA_identity
gene <- rep(NA, length(sgRNA_identity))
gene[sgRNA_identity == "NT1"] <- "NT"
gene[sgRNA_identity == "NT2"] <- "NT"
gene[sgRNA_identity == "AAVS"] <- "NT"
gene[sgRNA_identity %in% c("ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3")] <- "ITGB1"
gene[sgRNA_identity %in% c("RPS6-sg1", "RPS6-sg2", "RPS6-sg3")] <- "RPS6"
gene[sgRNA_identity %in% c("GATA1-sg1", "GATA1-sg2", "GATA1-sg3")] <- "GATA1"
gene <- factor(gene, levels = c("NT", "ITGB1", "RPS6", "GATA1"))
crispr_obj@meta.data$gene <- gene

table(crispr_obj@meta.data$gene)


DefaultAssay(crispr_obj) <- "Protein"

subset_obj <- subset(crispr_obj, subset = gene %in% c("NT", "GATA1"))


# Divide GATA1 sgRNA cells into top 30, middle, and bottom 30 based on GATA1 protein expression levels
gata1_expression <- FetchData(subset_obj, vars = c("GATA1-pAbO", "gene")) %>%
  dplyr::filter(gene == "GATA1")
sorted_cells <- gata1_expression[order(gata1_expression$`GATA1-pAbO`, decreasing = TRUE), , drop = FALSE]
num_cells <- nrow(sorted_cells)
top_30 <- ceiling(0.3 * num_cells)
bottom_30 <- floor(0.3 * num_cells)
subset_obj$protein_30 <- "middle"
top_cells <- rownames(sorted_cells)[1:top_30]
subset_obj$protein_30[top_cells] <- "top_30"
bottom_cells <- rownames(sorted_cells)[(num_cells - bottom_30 + 1):num_cells]
subset_obj$protein_30[bottom_cells] <- "bottom_30"
subset_obj$protein_30[subset_obj$gene == "NT"] <- "NT"


# Delete cells where protein_30 is “top_30”
subset_obj <- subset(subset_obj, subset = protein_30 != "top_30")

# standard pre-processing
DefaultAssay(subset_obj) <- "RNA"
subset_obj = NormalizeData(subset_obj)
subset_obj = FindVariableFeatures(subset_obj)
subset_obj = ScaleData(subset_obj)
subset_obj = RunPCA(subset_obj)

# calculate Perturbation signatures 
subset_obj <- CalcPerturbSig(
  object = subset_obj, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="gene", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = "PRTB", 
  split.by = NULL)

# Mixscale
subset_obj = RunMixscale(
  object = subset_obj, 
  assay = "PRTB", 
  slot = "data", 
  labels = "gene", 
  nt.class.name = "NT", 
  min.de.genes = 5, 
  logfc.threshold = 0.25,
  de.assay = "RNA",
  max.de.genes = 100, 
  new.class.name = "mixscale_score", 
  fine.mode = F, 
  verbose = F, 
  split.by = NULL)


### Differential expression (DE) analysis
de_res = Run_wmvRegDE(object = subset_obj, assay = "RNA", slot = "counts",
                      labels = "gene", nt.class.name = "NT",
                      logfc.threshold = 0,
                      split.by = NULL,
                      full.results = F)

# Since the software did not strictly filter by logfc.threshold, apply logFC filtering manually
logfc_threshold <- 0
de_res <- lapply(de_res, function(df) {
  df[df$log2FC >= logfc_threshold | df$log2FC <= -logfc_threshold, ]
})

for (gene in names(de_res)) {
  file_name <- paste0(gene, "_DE_results_logfc0.tsv")
  write.table(de_res[[gene]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}


### GATA1_DE_results_logfc0 Calculate FDR
file_path <- "GATA1_DE_results_logfc0.tsv"
deg <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deg$FDR <- p.adjust(deg$p_weight, method = "BH")
write.table(deg, file = "GATA1_DE_results_logfc0_padjust.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



##### go_kegg_mixscale_GATA1cell
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(RColorBrewer)

### GATA1 cell / NT go_kegg
###DEGs
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/GATA1_NT_go_kegg")
DEGs <- read.table("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/GATA1_DE_results_logfc0_padjust.tsv", header = TRUE, sep = "\t")

#  Retain genes with p_weight < 0.01 and absolute log2FC > 0.5
DEGs <- DEGs %>% filter(p_weight < 0.01, abs(log2FC) > 0.5)

DEGs_all <- DEGs  
DEGs_positive <- subset(DEGs, log2FC > 0)  # log2FC > 0 
DEGs_negative <- subset(DEGs, log2FC < 0)  # log2FC < 0 

convert_to_entrez <- function(data) {
  gene_symbols <- data$gene_ID
  conversion <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(conversion)
}

DEGs_all_Entrez <- convert_to_entrez(DEGs_all)
DEGs_positive_Entrez <- convert_to_entrez(DEGs_positive)
DEGs_negative_Entrez <- convert_to_entrez(DEGs_negative)

DEGs_positive
write.table(DEGs_positive,file="GATA1_NT_DEGs_positive.tsv",sep="\t",quote=F,row.names = F)
DEGs_negative
write.table(DEGs_negative,file="GATA1_NT_DEGs_negative.tsv",sep="\t",quote=F,row.names = F)

group_Entrezid <- list(
  all = DEGs_all_Entrez$ENTREZID,
  positive = DEGs_positive_Entrez$ENTREZID,
  negative = DEGs_negative_Entrez$ENTREZID)

#BP
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="GATA1_NT_compareCluster_BP.tsv",sep="\t",quote=F,row.names = F)

#MF
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="GATA1_NT_compareCluster_MF.tsv",sep="\t",quote=F,row.names = F)

#CC
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="GATA1_NT_compareCluster_CC.tsv",sep="\t",quote=F,row.names = F)

#KEGG
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichKEGG",
                         organism="hsa", 
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500)
write.table(result,file="GATA1_NT_compareCluster_KEGG.tsv",sep="\t",quote=F,row.names = F)



##### go_kegg_mixscale_GATA1cell_removetop20
### GATA1cell_removeGATA1proteintop20 / NT DEG
###DEGs
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/go_kegg_mixscale_GATA1_top20middlebottoom20")
DEGs <- read.table("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/GATA1_DE_results_logfc0_padjust.tsv", header = TRUE, sep = "\t")

#  Retain genes with p_weight < 0.01 and absolute log2FC > 0.5
DEGs <- DEGs %>% filter(p_weight < 0.01, abs(log2FC) > 0.5)

DEGs_all <- DEGs
DEGs_positive <- subset(DEGs, log2FC > 0)  # log2FC > 0
DEGs_negative <- subset(DEGs, log2FC < 0)  # log2FC < 0


convert_to_entrez <- function(data) {
  gene_symbols <- data$gene_ID  
  conversion <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(conversion)
}


DEGs_all_Entrez <- convert_to_entrez(DEGs_all)
DEGs_positive_Entrez <- convert_to_entrez(DEGs_positive)
DEGs_negative_Entrez <- convert_to_entrez(DEGs_negative)

DEGs_positive
write.table(DEGs_positive,file="GATA1_NT_DEGs_positive.tsv",sep="\t",quote=F,row.names = F)
DEGs_negative
write.table(DEGs_negative,file="GATA1_NT_DEGs_negative.tsv",sep="\t",quote=F,row.names = F)

group_Entrezid <- list(
  all = DEGs_all_Entrez$ENTREZID,
  positive = DEGs_positive_Entrez$ENTREZID,
  negative = DEGs_negative_Entrez$ENTREZID)

#BP
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="mixscale_GATA1_compareCluster_BP.tsv",sep="\t",quote=F,row.names = F)

#MF
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="mixscale_GATA1_compareCluster_MF.tsv",sep="\t",quote=F,row.names = F)

#CC
result <- compareCluster(geneClusters = group_Entrezid,
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="mixscale_GATA1_compareCluster_CC.tsv",sep="\t",quote=F,row.names = F)

#KEGG
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichKEGG",
                         organism="hsa", 
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500)
write.table(result,file="mixscale_GATA1_compareCluster_KEGG.tsv",sep="\t",quote=F,row.names = F)



##### go_kegg_mixscale_GATA1cell_removetop30
### GATA1cell_removeGATA1proteintop30 / NT DEG
###DEGs
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/go_kegg_mixscale_GATA1_top30middlebottoom30")
DEGs <- read.table("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/GATA1_DE_results_logfc0_padjust.tsv", header = TRUE, sep = "\t")

#  Retain genes with p_weight < 0.01 and absolute log2FC > 0.5
DEGs <- DEGs %>% filter(p_weight < 0.01, abs(log2FC) > 0.5)

DEGs_all <- DEGs 
DEGs_positive <- subset(DEGs, log2FC > 0)  # log2FC > 0 
DEGs_negative <- subset(DEGs, log2FC < 0)  # log2FC < 0 

convert_to_entrez <- function(data) {
  gene_symbols <- data$gene_ID  # 选择基因列
  conversion <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(conversion)
}

DEGs_all_Entrez <- convert_to_entrez(DEGs_all)
DEGs_positive_Entrez <- convert_to_entrez(DEGs_positive)
DEGs_negative_Entrez <- convert_to_entrez(DEGs_negative)

DEGs_positive
write.table(DEGs_positive,file="GATA1_NT_DEGs_positive.tsv",sep="\t",quote=F,row.names = F)
DEGs_negative
write.table(DEGs_negative,file="GATA1_NT_DEGs_negative.tsv",sep="\t",quote=F,row.names = F)

group_Entrezid <- list(
  all = DEGs_all_Entrez$ENTREZID,
  positive = DEGs_positive_Entrez$ENTREZID,
  negative = DEGs_negative_Entrez$ENTREZID)

#BP
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="mixscale_GATA1_compareCluster_BP.tsv",sep="\t",quote=F,row.names = F)

#MF
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="mixscale_GATA1_compareCluster_MF.tsv",sep="\t",quote=F,row.names = F)

#CC
result <- compareCluster(geneClusters = group_Entrezid,
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db, 
                         ont = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = T)
write.table(result,file="mixscale_GATA1_compareCluster_CC.tsv",sep="\t",quote=F,row.names = F)

#KEGG
result <- compareCluster(geneClusters = group_Entrezid, 
                         fun = "enrichKEGG",
                         organism="hsa", 
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         minGSSize = 10,
                         maxGSSize = 500)
write.table(result,file="mixscale_GATA1_compareCluster_KEGG.tsv",sep="\t",quote=F,row.names = F)



### go_plot
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GO/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GO/")

### negative

# All
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/GATA1_NT_go_kegg/GATA1_NT_compareCluster_CC.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "hemoglobin complex"
)
item <- go_data %>% filter(Cluster == "negative" & Description %in% target_descriptions)
item$group <- "All"
go_negative <- item

# remove_top20
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/go_kegg_mixscale_GATA1_top20middlebottoom20/mixscale_GATA1_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "ribosome biogenesis",
  "hemoglobin metabolic process"
)
item <- go_data %>% filter(Cluster == "negative" & Description %in% target_descriptions)
item$group <- "remove_top20"
go_negative <- rbind(go_negative, item)

file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/go_kegg_mixscale_GATA1_top20middlebottoom20/mixscale_GATA1_compareCluster_CC.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "hemoglobin complex"
)
item <- go_data %>% filter(Cluster == "negative" & Description %in% target_descriptions)
item$group <- "remove_top20"
go_negative <- rbind(go_negative, item)

# remove_top30
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/go_kegg_mixscale_GATA1_top30middlebottoom30/mixscale_GATA1_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "ribosome biogenesis",
  "hemoglobin metabolic process",
  "heme metabolic process",
  "heme biosynthetic process",
  "hemoglobin biosynthetic process"
)
item <- go_data %>% filter(Cluster == "negative" & Description %in% target_descriptions)
item$group <- "remove_top30"
go_negative <- rbind(go_negative, item)

file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/go_kegg_mixscale_GATA1_top30middlebottoom30/mixscale_GATA1_compareCluster_CC.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "hemoglobin complex"
)
item <- go_data %>% filter(Cluster == "negative" & Description %in% target_descriptions)
item$group <- "remove_top30"
go_negative <- rbind(go_negative, item)

write.table(go_negative,
            file = "go_negative_3group.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

##### Figure1 k
# DownregulatedGOtermsinGATA1sg_3group
go_negative$group <- factor(go_negative$group, levels = c("All", "remove_top20", "remove_top30"))
go_negative <- go_negative[order(go_negative$Description, go_negative$group), ]
go_negative <- rbind(go_negative[4:nrow(go_negative), ], go_negative[1:3, ])

go_negative$y_label <- factor(seq_len(nrow(go_negative)), levels = rev(seq_len(nrow(go_negative))))

group_colors <- c(
  "All" = "#A6CEE3",         
  "remove_top20" = "#5FA6CE", 
  "remove_top30" = "#1F78B4"
)

go_negative$label_combined <- paste0(
  go_negative$Description,
  str_dup(" ", 5),
  "p.adjust=", signif(go_negative$p.adjust, 2)
)

ggplot(go_negative, aes(x = Count, y = y_label, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label_combined),
            x = 0.1, hjust = 0, vjust = 0.5, size = 4, color = "black") +
  geom_text(aes(label = geneID),
            x = 0.1, hjust = 0, vjust = 4, size = 3.5,
            color = "#75B9D3") +
  scale_fill_manual(values = group_colors) +
  scale_y_discrete(labels = NULL) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  ) +
  labs(title = "Downregulated GO terms in GATA1sg", x = "Count", y = "") +
  coord_cartesian(clip = "off", expand = TRUE)
ggsave(filename = "DownregulatedGOtermsinGATA1sg_3group.pdf", plot = last_plot(), width = 10, height = 9, dpi = 300)
