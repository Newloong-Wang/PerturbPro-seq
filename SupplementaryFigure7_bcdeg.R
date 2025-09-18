library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(paletteer)
library(viridis)
library(ggridges)

setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/figure/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")
K562<- Read10X("/public/home/wangxinlong/Project/crispr1/data/20250102ABFC20221434-76_K562/sc20241231K562_addsequence/Ab-1_L2_S2401B2401_DBEC_MolsPerCell_MEX")
K562_seurat<-CreateSeuratObject(counts =K562$`Gene Expression`)


##### SupplementaryFigure7 b
VlnPlot(crispr_obj, features = "nCount_sgRNA", group.by = "orig.ident", pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white") + 
  scale_fill_manual(values ='#65B8AD' ) + 
  NoLegend()
ggsave("allCell_nCount_sgRNA.pdf", plot = last_plot(), width = 4, height = 5, units = "in", dpi = 300)

##### SupplementaryFigure7 c
VlnPlot(K562_seurat,features=c("nCount_RNA"),ncol=1, pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")+ 
  scale_fill_manual(values ='#65B8AD' ) + 
  NoLegend()
median(K562_seurat$nCount_RNA)
ggsave("Nothresholdapplied_allCell_12468.nCount_RNA.pdf", plot = last_plot(), width = 4, height = 5, units = "in", dpi = 300)





##### SupplementaryFigure7 de
##### extract_nCountRNA_nFeatureRNA_percentmt
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

setwd("/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/seurat_meta/")

### process Chromium_standard_human
object <- Read10X("/public/home/wangxinlong/Project/crispr1/data/compare_seq_K562/Chromium_standard/cellranger_human_subsample_073/result/outs/filtered_feature_bc_matrix/")
obj <- CreateSeuratObject(counts = object)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
meta_info <- obj@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, percent.mt)
meta_info <- tibble::rownames_to_column(meta_info, var = "cell_id")
write.table(
  meta_info,
  file = "Chromiumstandardhuman_meta_info.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


### process PerturbPro
object <- Read10X("/public/home/wangxinlong/Project/crispr1/data/compare_seq_K562/PerturbPro/BD_subsample_074/074_results/WTA-2_subsampled_L2_Q0063W0072_RSEC_MolsPerCell_MEX/")
obj <- CreateSeuratObject(counts = object)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
meta_info <- obj@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, percent.mt)
meta_info <- tibble::rownames_to_column(meta_info, var = "cell_id")
write.table(
  meta_info,
  file = "PerturbPro_meta_info.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


### process Perturb
object <- Read10X("/public/home/wangxinlong/Project/crispr1/data/compare_seq_K562/Perturb/cellranger_count_A6_subsample_024/result/outs/filtered_feature_bc_matrix/")
obj <- CreateSeuratObject(counts = object)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
meta_info <- obj@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, percent.mt)
meta_info <- tibble::rownames_to_column(meta_info, var = "cell_id")
write.table(
  meta_info,
  file = "Perturb_meta_info.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


### process ECCITE
object <- Read10X("/public/home/wangxinlong/Project/crispr1/data/compare_seq_K562/ECCITE/cellranger_subsample_022/result/outs/filtered_feature_bc_matrix/")
obj <- CreateSeuratObject(counts = object)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
meta_info <- obj@meta.data %>%
  dplyr::select(nCount_RNA, nFeature_RNA, percent.mt)
meta_info <- tibble::rownames_to_column(meta_info, var = "cell_id")
write.table(
  meta_info,
  file = "ECCITE_meta_info.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



##### bam_summary                   
### The following section of code uses Python
import pandas as pd


### PerturbPro
file_path = "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/PerturbPro_bam_info.tsv"
PerturbPro_bam = pd.read_csv(file_path, sep="\t")
PerturbPro_bam = PerturbPro_bam[PerturbPro_bam["HI"].isin([0, 1])]
# nRead_RNA
cell_counts = PerturbPro_bam["cell_id"].value_counts().reset_index()
cell_counts.columns = ["cell_id", "nRead_RNA"]
# nRead_intron
introns_df = PerturbPro_bam[PerturbPro_bam["xf_tag"].str.startswith("__introns")]

introns_counts = introns_df["cell_id"].value_counts().reset_index()
introns_counts.columns = ["cell_id", "nRead_intron"]

bam_summary = pd.merge(cell_counts, introns_counts, on="cell_id", how="left")
bam_summary["nRead_intron"] = bam_summary["nRead_intron"].fillna(0).astype(int)
bam_summary["percent_intron"] = ((bam_summary["nRead_intron"] / bam_summary["nRead_RNA"]) * 100).round(2)
bam_summary = bam_summary.sort_values(by="cell_id")
bam_summary.to_csv(
    "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/PerturbPro_bam_summary.tsv",
    sep="\t",
    index=False
)


### Chromium_standard_human
file_path = "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/Chromiumstandardhuman_bam_info.tsv"
Chromiumstandardhuman_bam = pd.read_csv(file_path, sep="\t")
Chromiumstandardhuman_bam = Chromiumstandardhuman_bam.dropna(subset=["cell_id", "UMI"])
# nRead_RNA
cell_counts = Chromiumstandardhuman_bam["cell_id"].value_counts().reset_index()
cell_counts.columns = ["cell_id", "nRead_RNA"]
Chromiumstandardhuman_bam = Chromiumstandardhuman_bam.dropna(subset=["RE"])
# nRead_intron
introns_df = Chromiumstandardhuman_bam[Chromiumstandardhuman_bam["RE"] == "N"]

introns_counts = introns_df["cell_id"].value_counts().reset_index()
introns_counts.columns = ["cell_id", "nRead_intron"]

bam_summary = pd.merge(cell_counts, introns_counts, on="cell_id", how="left")
bam_summary["nRead_intron"] = bam_summary["nRead_intron"].fillna(0).astype(int)
bam_summary["percent_intron"] = ((bam_summary["nRead_intron"] / bam_summary["nRead_RNA"]) * 100).round(2)
bam_summary = bam_summary.sort_values(by="cell_id")
bam_summary.to_csv(
    "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/Chromiumstandardhuman_bam_summary.tsv",
    sep="\t",
    index=False
)


### Perturb
file_path = "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/Perturb_bam_info.tsv"
Perturb_bam = pd.read_csv(file_path, sep="\t")
Perturb_bam = Perturb_bam.dropna(subset=["cell_id", "UMI"])
# nRead_RNA
cell_counts = Perturb_bam["cell_id"].value_counts().reset_index()
cell_counts.columns = ["cell_id", "nRead_RNA"]
Perturb_bam = Perturb_bam.dropna(subset=["RE"])
# nRead_intron
introns_df = Perturb_bam[Perturb_bam["RE"] == "N"]

introns_counts = introns_df["cell_id"].value_counts().reset_index()
introns_counts.columns = ["cell_id", "nRead_intron"]

bam_summary = pd.merge(cell_counts, introns_counts, on="cell_id", how="left")
bam_summary["nRead_intron"] = bam_summary["nRead_intron"].fillna(0).astype(int)
bam_summary["percent_intron"] = ((bam_summary["nRead_intron"] / bam_summary["nRead_RNA"]) * 100).round(2)
bam_summary = bam_summary.sort_values(by="cell_id")
bam_summary.to_csv(
    "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/Perturb_bam_summary.tsv",
    sep="\t",
    index=False
)


### ECCITE
file_path = "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/ECCITE_bam_info.tsv"
ECCITE_bam = pd.read_csv(file_path, sep="\t")
ECCITE_bam = ECCITE_bam.dropna(subset=["cell_id", "UMI"])
# nRead_RNA
cell_counts = ECCITE_bam["cell_id"].value_counts().reset_index()
cell_counts.columns = ["cell_id", "nRead_RNA"]
ECCITE_bam = ECCITE_bam.dropna(subset=["RE"])
# nRead_intron
introns_df = ECCITE_bam[ECCITE_bam["RE"] == "N"]
introns_counts = introns_df["cell_id"].value_counts().reset_index()
introns_counts.columns = ["cell_id", "nRead_intron"]

bam_summary = pd.merge(cell_counts, introns_counts, on="cell_id", how="left")
bam_summary["nRead_intron"] = bam_summary["nRead_intron"].fillna(0).astype(int)
bam_summary["percent_intron"] = ((bam_summary["nRead_intron"] / bam_summary["nRead_RNA"]) * 100).round(2)
bam_summary = bam_summary.sort_values(by="cell_id")
bam_summary.to_csv(
    "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/ECCITE_bam_summary.tsv",
    sep="\t",
    index=False
)



### The following section of code uses R
##### merge nCountRNA_nFeatureRNA_percentmt and bam_summary
setwd("/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/merge/")


### Chromium_standard_human
bam_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/Chromiumstandardhuman_bam_summary.tsv"
meta_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/seurat_meta/Chromiumstandardhuman_meta_info.tsv"
Chromiumstandardhuman_bam <- read.delim(bam_file, stringsAsFactors = FALSE)
Chromiumstandardhuman_meta <- read.delim(meta_file, stringsAsFactors = FALSE)
Chromiumstandardhuman_merged <- inner_join(Chromiumstandardhuman_meta, Chromiumstandardhuman_bam, by = "cell_id")
write.table(
  Chromiumstandardhuman_merged,
  file = "Chromiumstandardhuman_merged.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


### PerturbPro
bam_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/PerturbPro_bam_summary.tsv"
meta_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/seurat_meta/PerturbPro_meta_info.tsv"
PerturbPro_bam <- read.delim(bam_file, stringsAsFactors = FALSE)
PerturbPro_meta <- read.delim(meta_file, stringsAsFactors = FALSE)
PerturbPro_merged <- inner_join(PerturbPro_meta, PerturbPro_bam, by = "cell_id")
write.table(
  PerturbPro_merged,
  file = "PerturbPro_merged.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


### Perturb
bam_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/Perturb_bam_summary.tsv"
meta_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/seurat_meta/Perturb_meta_info.tsv"
Perturb_bam <- read.delim(bam_file, stringsAsFactors = FALSE)
Perturb_meta <- read.delim(meta_file, stringsAsFactors = FALSE)
Perturb_merged <- inner_join(Perturb_meta, Perturb_bam, by = "cell_id")
write.table(
  Perturb_merged,
  file = "Perturb_merged.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


### ECCITE
bam_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/extract_bam/ECCITE_bam_summary.tsv"
meta_file <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/seurat_meta/ECCITE_meta_info.tsv"
ECCITE_bam <- read.delim(bam_file, stringsAsFactors = FALSE)
ECCITE_meta <- read.delim(meta_file, stringsAsFactors = FALSE)
ECCITE_merged <- inner_join(ECCITE_meta, ECCITE_bam, by = "cell_id")
write.table(
  ECCITE_merged,
  file = "ECCITE_merged.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


### merge 4 data
PerturbPro_merged$cell_id <- as.character(PerturbPro_merged$cell_id)

Chromiumstandardhuman_merged <- Chromiumstandardhuman_merged %>%
  mutate(group = "Chromium")

PerturbPro_merged <- PerturbPro_merged %>%
  mutate(group = "PerturbPro")

Perturb_merged <- Perturb_merged %>%
  mutate(group = "Perturb")

ECCITE_merged <- ECCITE_merged %>%
  mutate(group = "ECCITE")

combined_merged <- bind_rows(
  Chromiumstandardhuman_merged,
  GSM4367984_merged,
  PerturbPro_merged,
  Perturb_merged,
  ECCITE_merged
)

write.table(
  combined_merged,
  file = "combined_merged_all.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



##### plot
setwd("/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/figure_2/")

file_path <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/merge/combined_merged_all.tsv"

combined_merged <- read.delim(file_path, stringsAsFactors = FALSE)

combined_merged$group <- factor(combined_merged$group, levels = c("Chromium", "ECCITE", "PerturbPro", "Perturb"))
combined_merged <- combined_merged %>% arrange(group)

custom_colors <- c(
  "Chromium" = "#F4A88E",
  "PerturbPro" = "#5E89BA",
  "Perturb" = "#9594D1",
  "ECCITE" = "#85C1BB"
)


##### SupplementaryFigure7 d
# x = nRead_RNA, y = nFeature_RNA
ggplot(combined_merged, aes(x = nRead_RNA, y = nFeature_RNA, color = group)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    x = "Read count",
    y = "Number of genes",
    color = "Group"
  ) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
ggsave("nReadRNA_nFeatureRNA.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)

# x = nRead_RNA, y = nCount_RNA
ggplot(combined_merged, aes(x = nRead_RNA, y = nCount_RNA, color = group)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    x = "Read count",
    y = "Number of UMIs",
    color = "Group"
  ) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
ggsave("nReadRNA_nCountRNA.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)


##### SupplementaryFigure7 e
# boxplot UMIs
ggplot(combined_merged, aes(x = group, y = nCount_RNA, fill = group)) +
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
    y = "Number of UMIs"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave("boxplot_UMIs.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)


# boxplot percent.mt
ggplot(combined_merged, aes(x = group, y = percent.mt, fill = group)) +
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
    y = "Percent_MT"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave("boxplot_Percent_MT.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)


# boxplot percent_intron
ggplot(combined_merged, aes(x = group, y = percent_intron, fill = group)) +
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
    y = "Percent_Intron"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave("boxplot_Percent_Intron.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)





##### SupplementaryFigure7 g
# Load packages.
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(zoo) 
library(tibble)
library(tidyr)


setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNAcounts_mRNA_sgRNAcounts_Protein")
# Load object.
crispr_obj <- readRDS("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")

DefaultAssay(crispr_obj) <- "Protein"
crispr_obj <- subset(crispr_obj,subset = nCount_Protein<150)

### Custom function
fun_cellnumber <- function(crispr_obj = crispr_obj,sgRNA_select,protein_select,rna_select){
    subset_obj <- subset(crispr_obj, subset = sgRNA_identity %in% c(sgRNA_select))
    sgRNA_data <- GetAssayData(subset_obj, assay = "sgRNA", slot = "counts")
    sgRNA_sum <- colSums(sgRNA_data)
    tmp <- sgRNA_data  %>% 
        as.data.frame()  %>% t() %>% 
        as.data.frame() %>% 
        rownames_to_column('cells')
    subset_obj$sgRNA_sum <- tmp[,sgRNA_select]#sgRNA_sum
    subset_obj$sgRNA_rank <- rank(-subset_obj$sgRNA_sum, ties.method = "first")
    # Extract normalized expression data for Protein and RNA
    DefaultAssay(subset_obj) <- "Protein"
    protein_expression <- GetAssayData(subset_obj, assay = "Protein", slot = "data")[protein_select, ]
    DefaultAssay(subset_obj) <- "RNA"
    rna_expression <- GetAssayData(subset_obj, assay = "RNA", slot = "data")[rna_select, ]
    # Extract metadata
    meta <- FetchData(subset_obj, vars = c("sgRNA_sum", "sgRNA_rank"))
    data <- cbind(meta, protein = protein_expression, rna = rna_expression)
    data <- data[order(data$sgRNA_rank), ]
    data_res <- data.frame(
        mRNA = (data %>% filter(rna>0) %>% nrow())/nrow(data),
        Protein = (data %>% filter(protein>0) %>% nrow())/nrow(data)
    ) %>% 
        pivot_longer(cols = names(.)) %>% 
        mutate(
            name = factor(name,levels = c('mRNA','Protein')),
            sgRNA = sgRNA_select,protein_select = protein_select,
            rna_select = rna_select
        )
    return(data_res)
}

data_use <- list(
    c('NT1','ITGB1','CD29-pAbO'),
    c('NT2','ITGB1','CD29-pAbO'),
    c('AAVS','ITGB1','CD29-pAbO'),
    c('NT1','GATA1','GATA1-pAbO'),
    c('NT2','GATA1','GATA1-pAbO'),
    c('AAVS','GATA1','GATA1-pAbO')
)

data_res <- lapply(data_use,function(x){
    # print(x[2])
    fun_cellnumber(crispr_obj = crispr_obj,sgRNA_select = x[1],rna_select = x[2],protein_select = x[3])
})

data_plot <- data_res %>% do.call(rbind,.)
data_plot

data_plot %>% group_by(name) %>% summarise(value_median = median(value))
p <- ggboxplot(
  data_plot, x="name", y="value",
  color ="name",
  width = 0.6,
  palette = paletteer_d("ggsci::nrc_npg"),
  add = "jitter",
  xlab = F,  bxp.errorbar=T,
  bxp.errorbar.width=0.5,
  size=0.5, outlier.shape=NA,legend = "right") +
  guides(color = guide_legend(title = 'Group'))+
  stat_compare_means(
    label = "p.format",size = 8,
    comparisons = list(c(
      "mRNA","Protein"
    )),
    method = "wilcox.test",
  ) +
  scale_color_manual(values = c('#B25B69','#4D83B8')) +
  # facet_wrap(~metabolites,nrow = 2,strip.position = 'left',scales = 'free') +
  theme(
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(1.2, "cm"),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 20,face = 'italic',vjust = 0.7),
    strip.placement = 'outside'
  )
p
ggsave("NT_Proportionofcells.pdf", plot = p, width = 7, height = 8, units = "in", dpi = 300)