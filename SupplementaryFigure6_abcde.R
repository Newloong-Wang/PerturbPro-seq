##### SupplementaryFigure6
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(paletteer)
library(writexl)
setwd('/public/home/wangxinlong/Project/crispr1/result/lung/')

###seurat
lung<-Read10X("/public/home/wangxinlong/Project/crispr1/data/20241217_sclung/Ab-5_L4_Q0028W0072_DBEC_MolsPerCell_MEX/")
lung_obj<-CreateSeuratObject(counts =lung$`Gene Expression`)
lung_obj[['Protein']] = CreateAssayObject(counts = lung$`Antibody Capture`, project="Lung")

Assays(lung_obj)
rownames(lung_obj[["Protein"]])
lung_obj[["percent.mt"]] <- PercentageFeatureSet(lung_obj, pattern = "^mt")



##### SupplementaryFigure6 a
VlnPlot(lung_obj,features=c("nFeature_RNA"),ncol=1, pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")+ 
  scale_fill_manual(values = "#75A3D1") + 
  NoLegend() +
  scale_y_continuous(limits = y_lims)
y_lims <- c(0,6000)
median(lung_obj$nFeature_RNA)
ggsave("allCell_nFeature_RNA_median1561.pdf", plot = last_plot(), width = 8, height = 8, units = "in", dpi = 300)
##### SupplementaryFigure6 a
VlnPlot(lung_obj,features=c("nCount_RNA"),ncol=1, pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")+ 
  scale_fill_manual(values = "#75A3D1") + 
  NoLegend()+
  scale_y_continuous(limits = y_lims)
y_lims <- c(0,10000)
median(lung_obj$nCount_RNA)
ggsave("allCell_nCount_RNA_median3334.5.pdf", plot = last_plot(), width = 8, height = 8, units = "in", dpi = 300)

plot1 <- FeatureScatter(lung_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

lung_obj_fliter <- subset(lung_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 20)

VlnPlot(lung_obj_fliter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

lung_obj_fliter<-NormalizeData(lung_obj_fliter,normalization.method = "LogNormalize",scale.factor = 10000)
lung_obj_fliter<-FindVariableFeatures(lung_obj_fliter,selection.method = "vst",nfeatures=700)

marker.list<-c("Vmf","Il1b","Aplnr","Gpihbp1","Cd209a","Flt3","Tbx2","Ednrb","Apln","Lamp3","Scgb1a1","Scgb3a2","Bpifb1","Sftpd","Abca3","Foxj1","Ccdc113","Sftpc","Car4","Vwf","Slc6a2","Thy1","Pdpn","Prox1","Flt4","Mmrn1","Fxyd6","Rtkn2","Ager","Hopx","Pdpn","Retnlg","Ear2","Marco","Abcg1","Dcn","Col14a1","Msln","Gucy1a1","Notch3","Enpp2","Pdgfrb","Pi16","Cygb","Acta2","Tagln","Npnt","Col13a1","Itga8","Cxcl14","Ms4a2","Plac8","Ly6c2","Ptprc","Epcam","Scgb1a1","Mki67","C1qa","C1qb","C1qc","Cd3g","Cd3d","Cd3e","Ms4a1","Bank1","Nkg7","Klrd1",'Hopx','Cd4','Cd19','Cd79a','Itgam','Cd7','Flt3','Hhip')
current_var_features <- VariableFeatures(lung_obj_fliter)
updated_var_features <- unique(c(current_var_features, marker.list))
VariableFeatures(lung_obj_fliter) <- updated_var_features
VariableFeatures(lung_obj_fliter)

top10<-head(VariableFeatures(lung_obj_fliter),10)
plot1<-VariableFeaturePlot(lung_obj_fliter)
plot2<-LabelPoints(plot=plot1,points = top10)
plot1+plot2

all.genes<-rownames(lung_obj_fliter)
lung_obj_fliter<-ScaleData(lung_obj_fliter,features=all.genes)
lung_obj_fliter<-RunPCA(lung_obj_fliter,features=VariableFeatures(object = lung_obj_fliter),npcs =100)
VizDimLoadings(lung_obj_fliter,dims=1:2,reduction = "pca")
ElbowPlot(lung_obj_fliter,ndims=80)

lung_obj_fliter<-FindNeighbors(lung_obj_fliter,dims=1:26)
lung_obj_fliter<-FindClusters(lung_obj_fliter,resolution = 1.5)
lung_obj_fliter<-RunUMAP(lung_obj_fliter,dims = 1:26)
DimPlot(lung_obj_fliter,reduction = "umap",label =TRUE)
table(lung_obj_fliter@meta.data[["RNA_snn_res.1.5"]])

FeaturePlot(lung_obj_fliter, features = c("nFeature_RNA")) #aM 17



##### SupplementaryFigure6 b
#Cell type annotation
new.cluster.ids <- c(
  "AT2",          # 0
  "gCap",         # 1
  "Col13a1",      # 2
  "Col13a1",      # 3
  "Pericyte",     # 4
  "aCap",         # 5
  "Club",         # 6
  "AT2",          # 7
  "T",            # 8
  "Col14a1",      # 9
  "Col13a1",      # 10
  "Monocyte",     # 11
  "V_SMC",        # 12
  "Ciliated",     # 13
  "AT2",          # 14
  "EC_Vein",      # 15
  "Pericyte",     # 16
  "aM",           # 17
  "NK",           # 18
  "EC_Artery",    # 19
  "Col13a1",      # 20
  "AT1",          # 21
  "Mesothelial",  # 22
  "B",            # 23
  "Innate_Lymphoid", # 24
  "Pericyte",     # 25
  "A_SMC",        # 26
  "Neutrophil",  # 27
  "EC_Lymph"      # 28
)
names(new.cluster.ids) <- levels(lung_obj_fliter)
lung_obj_fliter <- RenameIdents(lung_obj_fliter, new.cluster.ids)
lung_obj_fliter$celltype <- Idents(lung_obj_fliter)
DimPlot(lung_obj_fliter,reduction = "umap",label = TRUE, group.by = "celltype")

ordered.season<-c("AT1","AT2",'Ciliated',"Club","aCap","gCap","EC_Artery","EC_Vein",'EC_Lymph',"Pericyte",'Col13a1','Col14a1','V_SMC','A_SMC',"Mesothelial","T",'B','Neutrophil','aM','NK','Innate_Lymphoid','Monocyte')
Idents(lung_obj_fliter) <- factor(Idents(lung_obj_fliter), levels= ordered.season)

mypal <- c("#AED6F1", "#6BB3E6", "#75A3D1", "#1E75B3", 
           "#B23E4D", "#E3A8A2", "#CC7177", "#CF5B5B", "#93344F", 
           "#BDE0B6", "#77BB88", "#448855", "#65B8AD", "#3A7E75", "#2A5B54", 
           "#f8b595", "#B65A20", "#643112", "#B5420D","#F86930", "#CE4A0D", 
           '#F59061')
DimPlot(
  lung_obj_fliter,
  reduction = "umap",
  label = TRUE, 
  pt.size = 1, 
  cols = mypal, 
  label.size = 4, 
  label.box = FALSE, 
  repel = TRUE
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")
  ) +
  guides(col = guide_legend(ncol = 1))
ggsave("UMAP1.pdf", plot = last_plot(), width = 11, height = 9, units = "in", dpi = 300)



# celltype_cellbarcode
cellbarcode <- rownames(lung_obj_fliter@meta.data)
celltype <- lung_obj_fliter@meta.data$celltype
cell_info <- data.frame(cellbarcode = cellbarcode,
                        celltype = celltype,
                        stringsAsFactors = FALSE)
cell_summary <- aggregate(cellbarcode ~ celltype, data = cell_info,
                          FUN = function(x) paste(x, collapse = ","))
cell_summary$celltype <- factor(cell_summary$celltype, levels = ordered.season)
cell_summary <- cell_summary[order(cell_summary$celltype), ]
write.table(cell_summary, file = "celltype_cellbarcode.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write_xlsx(cell_summary, path = "celltype_cellbarcode.xlsx")

DefaultAssay(lung_obj_fliter) <- "RNA"
marker.genes<-FindAllMarkers(lung_obj_fliter, group.by = "celltype", only.pos = T,min.pct = 0.3,logfc.threshold = 1)
marker.genes <- marker.genes %>% filter(p_val_adj < 0.05)


                          
##### SupplementaryFigure6 c
# Count the marker genes for each cell type
marker_counts <- marker.genes %>%
  group_by(cluster) %>%
  summarise(marker_number = n()) %>%
  as.data.frame()

ggplot(marker_counts, aes(x = cluster, y = marker_number, fill = cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mypal) +
  labs(x = "Cell Type", y = "Number of Marker Genes", title = "Marker Gene Count per Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

marker_counts$cluster <- factor(marker_counts$cluster, levels = rev(ordered.season))

ggplot(marker_counts, aes(y = cluster, x = marker_number, color = cluster)) +
  geom_segment(aes(x = 0, xend = marker_number, y = cluster, yend = cluster), 
               color = "gray60", linewidth = 1) +  
  geom_point(size = 4) +  
  scale_color_manual(values = mypal[1:length(ordered.season)]) +  
  labs(x = "Number of Marker Genes", y = "Cell Type", title = "Marker Gene Count per Cell Type") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),  
    axis.text.x = element_text(size = 10),  
    legend.position = "none"  
  )
ggsave("marker_bangbangtu.pdf", width = 8, height = 6, units = "in", dpi = 300)

markergenes_df <- marker.genes[, c("cluster", "gene")]
marker_summary <- aggregate(gene ~ cluster, data = markergenes_df,
                            FUN = function(x) paste(x, collapse = ","))
colnames(marker_summary)[1] <- "celltype"
marker_summary$celltype <- factor(marker_summary$celltype, levels = ordered.season)
marker_summary <- marker_summary[order(marker_summary$celltype), ]
write.table(marker_summary, file = "markergenes.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


##### SupplementaryFigure6 c
# GO enrichment analysis
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(dplyr)

get_go_enrichment <- function(gene_list) {
  ego <- enrichGO(
    gene         = gene_list,
    OrgDb        = org.Mm.eg.db,  
    keyType      = "SYMBOL",      
    ont          = "BP",          
    pAdjustMethod = "BH",         
    pvalueCutoff  = 0.05, 
    qvalueCutoff  = 0.05
  )
  return(ego)
}
                            
go_results_list <- list()
for (cell_type in unique(marker.genes$cluster)) {
  genes <- marker.genes %>% filter(cluster == cell_type) %>% pull(gene)
  go_results <- get_go_enrichment(genes)
  go_results_list[[cell_type]] <- go_results
}

go_all_df <- do.call(rbind, lapply(names(go_results_list), function(cell_type) {
  go_result <- go_results_list[[cell_type]]
  if (!is.null(go_result) && nrow(go_result@result) > 0) {
    df <- go_result@result
    df$cluster <- cell_type
    return(df)
  } else {
    return(NULL)
  }
}))
write.table(go_all_df, file = "go_allcelltype_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

ordered_go_terms <- c(
  "AT1" = "cell-cell junction organization",
  "AT2" = "purine nucleotide metabolic process",
  "Ciliated" = "cilium organization",
  "Club" = "glutathione metabolic process",
  "gCap" = "regulation of angiogenesis",
  "aCap" = "blood vessel endothelial cell migration",
  "EC_Artery" = "regulation of blood circulation",
  "EC_Vein" = "endothelial cell migration",
  "EC_Lymph" = "regulation of vasculature development",
  "Pericyte" = "cell junction assembly",
  "Col13a1" = "extracellular matrix organization",
  "Col14a1" = "external encapsulating structure organization",
  "V_SMC" = "regulation of smooth muscle contraction",
  "A_SMC" = "actin-mediated cell contraction",
  "Mesothelial" = "extracellular matrix organization",
  "T" = "T cell differentiation",
  "B" = "B cell activation",
  "Neutrophil" = "myeloid leukocyte migration",
  "aM" = "phagocytosis",
  "NK" = "natural killer cell mediated cytotoxicity",
  "Innate_Lymphoid" = "cell activation involved in immune response",
  "Monocyte" = "mononuclear cell proliferation"
)

go_map <- data.frame(
  cluster = names(ordered_go_terms),
  target_description = unname(ordered_go_terms),
  stringsAsFactors = FALSE
)

go_selected_exact <- go_all_df %>%
  inner_join(go_map, by = c("cluster" = "cluster")) %>%
  filter(Description == target_description)

go_selected_exact$log10padj <- -log10(go_selected_exact$p.adjust)

mypal <- c("#AED6F1", "#6BB3E6", "#75A3D1", "#1E75B3", 
           "#B23E4D", "#E3A8A2", "#CC7177", "#CF5B5B", "#93344F", 
           "#BDE0B6", "#77BB88", "#448855", "#65B8AD", "#3A7E75", "#2A5B54", 
           "#f8b595", "#B65A20", "#643112", "#B5420D","#F86930", "#CE4A0D", 
           '#F59061')

go_selected_exact$cluster <- factor(go_selected_exact$cluster, levels = c(
  "AT1", "AT2", "Ciliated", "Club", "gCap", "aCap", "EC_Artery", "EC_Vein", 
  "EC_Lymph", "Pericyte", "Col13a1", "Col14a1", "V_SMC", "A_SMC", "Mesothelial", 
  "T", "B", "Neutrophil", "aM", "NK", "Innate_Lymphoid", "Monocyte"
))

ggplot(go_selected_exact, aes(x = log10padj, y = cluster, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) +  
  scale_fill_manual(values = mypal) +
  labs(x = "-log10(p.adjust)", y = "Cell Type", title = "GO Terms per Cell Type") + 
  theme_classic() + 
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) + 
  geom_text(
    aes(x = 0, label = Description),
    hjust = 0,
    nudge_x = 0.5,
    size = 3,
    color = "black"
  ) +   
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + 
  scale_y_discrete(limits = rev(c(
    "AT1", "AT2", "Ciliated", "Club", "gCap", "aCap", "EC_Artery", "EC_Vein", 
    "EC_Lymph", "Pericyte", "Col13a1", "Col14a1", "V_SMC", "A_SMC", "Mesothelial", 
    "T", "B", "Neutrophil", "aM", "NK", "Innate_Lymphoid", "Monocyte"
  )))
ggsave("marker_GO.pdf", width = 8, height = 6, units = "in", dpi = 300)

                            

##### SupplementaryFigure6 d
# marker gene ----
ordered.season <- c("AT1", "AT2", "Ciliated", "Club", 
                    "aCap", "gCap", "EC_Artery", "EC_Vein", "EC_Lymph", 
                    "Pericyte", "Col13a1", "Col14a1", "V_SMC", "A_SMC", 
                    "Mesothelial", "T", "B", "Neutrophil", "aM", "NK", 
                    "Innate_Lymphoid", "Monocyte")

Idents(lung_obj_fliter) <- factor(Idents(lung_obj_fliter), levels = ordered.season)
lung_obj_fliter@commands$NormalizeData.RNA 

DotPlot(lung_obj_fliter, 
        features = c(
          "Ager", "Hopx",      # AT1
          "Lamp3", "Sftpc",    # AT2
          "Foxj1", "Rsph1",    # Ciliated
          "Scgb1a1", "Scgb3a2",# Club
          "Ednrb", "Car4",     # aCap
          "Aplnr", "Gpihbp1",  # gCap
          "Gja5", "Olfm2",     # EC_Artery
          "Slc6a2", "Ephb4",   # EC_Vein
          "Prox1", "Flt4",     # EC_Lymph
          "Gucy1a1", "Notch3", # Pericyte
          "Cxcl14", "Col13a1", # Col13a1
          "Col14a1", "Pi16",   # Col14a1
          "Cnn1", "Ntrk3",     # V_SMC
          "Igf1", "Hhip",      # A_SMC
          "Msln", "Wt1",       # Mesothelial
          "Cd3e", "Cd4",       # T 
          "Cd79a", "Cd19",     # B 
          "Il1b", "Retnlg",    # Neutrophil
          "Ear2", "Abcg1",     # aM
          "Nkg7", "Klrd1",     # NK
          "Gata3", "Cd69",     # Innate Lymphoid Cells
          "Itgam", "Csf1r"     # Monocyte
        )) + 
  scale_color_gradient(low = "white", high = "#B23E4D") +
  coord_flip() + 
  theme_minimal() + 
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(), 
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5)
  )
ggsave("dotplot_2markerGenes.pdf", plot = last_plot(), width = 10, height = 12, units = "in", dpi = 300)



##### SupplementaryFigure6 e
# Protein ----
library(scales)
library(colorspace)
library(ggsignif)
                            
DefaultAssay(lung_obj_fliter) <- "Protein"
lung_obj_fliter <- NormalizeData(lung_obj_fliter, normalization.method = "CLR", margin = 2)                          
fill_colors <- c("#E3A8A2", "#75A3D1", "#E69D70", "#65B8AD")
border_colors <- darken(fill_colors, amount = 0.3)

lung_obj_fliter$cell_type <- ifelse(lung_obj_fliter$celltype %in% c("AT1", "AT2", "Ciliated", "Club"), "Epithelial", 
                                    ifelse(lung_obj_fliter$celltype %in% c("aCap", "gCap", "EC_Artery", "EC_Vein", "EC_Lymph"), "Endothelial", 
                                           ifelse(lung_obj_fliter$celltype %in% c("Col13a1", "Col14a1", "V_SMC", "A_SMC", "Mesothelial",'Pericyte'), "Mesenchymal", 
                                                  "Immune")))
lung_obj_fliter$celltype <- factor(lung_obj_fliter$celltype, levels = ordered.season)


# Epcam Protein
expression_data <- FetchData(lung_obj_fliter, vars = c("Epcam-940269-pAbO", "cell_type"))
comparisons <- list(
  c("Epithelial", "Endothelial"),
  c("Epithelial", "Immune"),
  c("Epithelial", "Mesenchymal")
)
ggplot(expression_data, aes(x = cell_type, y = `Epcam-940269-pAbO`, fill = cell_type)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  NoLegend()
ggsave("vln_Epcam_Pro_1_expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# CD31 Protein
expression_data <- FetchData(lung_obj_fliter, vars = c("CD31-940068-pAbO", "cell_type"))
comparisons <- list(
  c("Endothelial", "Epithelial"),
  c("Endothelial", "Immune"),
  c("Endothelial", "Mesenchymal")
)
ggplot(expression_data, aes(x = cell_type, y = `CD31-940068-pAbO`, fill = cell_type)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  NoLegend()
ggsave("vln_CD31_Pro_1-expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# CD11b Protein Neutrophil / otherImmune
IMM_clusters <- subset(lung_obj_fliter, subset = cell_type == "Immune")
expression_data <- FetchData(IMM_clusters, vars = c("CD11b-940380-pAbO", "celltype"))
comparisons <- list(
  c("Neutrophil", "T"),
  c("Neutrophil", "B"),
  c("Neutrophil", "NK"),
  c("Neutrophil", "Innate_Lymphoid"),
  c("Neutrophil", "aM"),
  c("Neutrophil", "Monocyte")
)
ggplot(expression_data, aes(x = celltype, y = `CD11b-940380-pAbO`, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  scale_x_discrete(limits = c("T", "B", "NK", "Innate_Lymphoid", "aM", "Neutrophil", "Monocyte")) +
  NoLegend()
ggsave("vln_CD11b_Pro_expression_NeutrophilvsOther.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)


# CD11b Protein Monocyte / otherImmune
IMM_clusters <- subset(lung_obj_fliter, subset = cell_type == "Immune")
expression_data <- FetchData(IMM_clusters, vars = c("CD11b-940380-pAbO", "celltype"))
comparisons <- list(
  c("Monocyte", "T"),
  c("Monocyte", "B"),
  c("Monocyte", "NK"),
  c("Monocyte", "Innate_Lymphoid"),
  c("Monocyte", "aM"),
  c("Monocyte", "Neutrophil")
)
ggplot(expression_data, aes(x = celltype, y = `CD11b-940380-pAbO`, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  scale_x_discrete(limits = c("T", "B", "NK", "Innate_Lymphoid", "aM", "Neutrophil", "Monocyte")) +
  NoLegend()
ggsave("vln_CD11b_Pro_expression_MonocytevsOther.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)


# NKX2.1 Protein
expression_data <- FetchData(lung_obj_fliter, vars = c("NKX2.1-940213-pAbO", "cell_type"))
comparisons <- list(
  c("Epithelial", "Endothelial"),
  c("Epithelial", "Immune"),
  c("Epithelial", "Mesenchymal")
)
ggplot(expression_data, aes(x = cell_type, y = `NKX2.1-940213-pAbO`, fill = cell_type)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  NoLegend()
ggsave("vln_NKX2.1_Pro_expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# Hopx Protein
Epi_clusters <- subset(lung_obj_fliter, subset = cell_type == "Epithelial")
expression_data <- FetchData(Epi_clusters, vars = c("Hopx-940299-pAbO", "celltype"))
comparisons <- list(
  c("AT1", "AT2"),
  c("AT1", "Ciliated"),
  c("AT1", "Club")
)
ggplot(expression_data, aes(x = celltype, y = `Hopx-940299-pAbO`, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  NoLegend()
ggsave("vln_Hopx_Pro_expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# Scgb1a1 Protein
Epi_clusters <- subset(lung_obj_fliter, subset = cell_type == "Epithelial")
expression_data <- FetchData(Epi_clusters, vars = c("Scgb1a1-940272-pAbO", "celltype"))
comparisons <- list(
  c("Club", "AT1"),
  c("Club", "AT2"),
  c("Club", "Ciliated")
)
ggplot(expression_data, aes(x = celltype, y = `Scgb1a1-940272-pAbO`, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  NoLegend()
ggsave("vln_Scgb1a1_Pro_expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)
