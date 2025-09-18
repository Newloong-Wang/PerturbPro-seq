##### SupplementaryFigure8 a
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
scale_to_range <- function(x, new_min = -2, new_max = 2) {
  old_min <- min(x, na.rm = TRUE)
  old_max <- max(x, na.rm = TRUE)
  scaled <- (x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
  return(scaled)
}
data_plot <- function(crispr_obj = crispr_obj,sgRNA_identity_select = NULL,sgRNA_select,protein_select,rna_select,window_size = NULL){
    sgRNA_identity_select <- if(is.null(sgRNA_identity_select)){
        sgRNA_identity
    }else{
        sgRNA_identity_select
    }
    cell_id <- crispr_obj@meta.data %>% 
        filter(grepl(pattern = paste(sgRNA_identity_select,collapse = '|'),sgRNA_identity)) %>% 
        rownames()
    subset_obj <- subset(crispr_obj, cells = cell_id)
    sgRNA_data <- GetAssayData(subset_obj, assay = "sgRNA", slot = "counts")
    sgRNA_sum <- colSums(sgRNA_data)
    tmp <- sgRNA_data  %>% 
        as.data.frame()  %>% t() %>% 
        as.data.frame() %>% 
        rownames_to_column('cells')
    subset_obj$sgRNA_sum <- tmp %>% select(all_of(sgRNA_select)) %>% rowSums()#sgRNA_sum
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
    # data %>% head()
    print('鉴定为sgRNA的细胞数：')
    print(nrow(data))
    print('sgRNA有表达的细胞数：')
    table(sgRNA_sum>0) %>% print()
    print('对应RNA有表达的细胞数：')
    table(data$rna>0) %>% print()
    print('对应Protein有表达的细胞数：')
    table(data$protein>0)%>% print()
    if(is.na(window_size)){
        window_size <- (151 / 1000 * nrow(data)) %>%
            round(0) %>% { if (. %% 2 == 0) . + 1 else . }
    }
    data_plot <- data %>% 
        mutate(
            rna_smoothed = rollmedian(rna, k = window_size, fill = NA, align = "center"),
            protein_smoothed = rollmedian(protein, k = window_size, fill = NA, align = "center")
        ) %>% 
        mutate(across(c(protein,rna,rna_smoothed, protein_smoothed), ~ scale_to_range(.x, -2, 2)))
    p1 <- ggplot(data_plot, aes(x = sgRNA_rank)) +
        geom_line(aes(y = sgRNA_sum), size = 1,color = '#663399') +  
        geom_point(aes(y = rescale(rna, to = range(sgRNA_sum))), size = 0.3, alpha = 0.3,color = '#49525E') +  
        geom_smooth(
            aes(y = rescale(rna_smoothed, to = range(sgRNA_sum))),
            method = "loess", span = 0.9,formula = 'y ~ x',
            color = '#A90C38',
            show.legend = FALSE,se = FALSE
        ) +
        # geom_line(aes(y = rescale(rna_smoothed, to = range(sgRNA_sum)), color = "mRNA"), size = 1) +  
        # geom_point(aes(y = sgRNA_sum), size = 0.5, alpha = 0.5) +  
        # geom_point(aes(y = rescale(rna_smoothed, to = range(sgRNA_sum))), size = 1, alpha = 0.5,color = '#49525E') +  
        scale_y_continuous(
          name = "sgRNA counts", 
          sec.axis = sec_axis(~ rescale(., from = range(sgRNA_sum), to = range(data_plot$rna_smoothed,na.rm = TRUE)), name = "mRNA")
        ) +
        # scale_color_manual(values = c("sgRNA counts" = "#DDA0DD", "mRNA" = "#87CEFA")) +
        labs(
          x = "Cells ranked by sgRNA counts",
          color = NULL,  
          title = paste(sgRNA_select," - mRNA (Window =", window_size, ")")
        ) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
          axis.title.y.left = element_text(color = "#663399", size = 12),
          axis.title.y.right = element_text(color = "#A90C38", size = 12),  
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "bottom",  
          legend.text = element_text(size = 10)
        )
    p2 <- ggplot(data_plot, aes(x = sgRNA_rank)) +
        geom_line(aes(y = sgRNA_sum), size = 1,color = '#663399') +  
        geom_point(aes(y = rescale(protein, to = range(sgRNA_sum))), size = 0.3, alpha = 0.3,color = '#49525E') +  
        geom_smooth(
            aes(y = rescale(protein_smoothed, to = range(sgRNA_sum))),
            method = "loess", span = 0.9,formula = 'y ~ x',
            color = '#2E5A87',
            show.legend = FALSE,se = FALSE
        ) +
        # geom_line(aes(y = rescale(rna_smoothed, to = range(sgRNA_sum)), color = "mRNA"), size = 1) +  
        # geom_point(aes(y = sgRNA_sum), size = 0.5, alpha = 0.5) +  
        # geom_point(aes(y = rescale(protein_smoothed, to = range(sgRNA_sum))), size = 1, alpha = 0.5,color = '#49525E') +  
        scale_y_continuous(
          name = "sgRNA counts", 
          sec.axis = sec_axis(~ rescale(., from = range(sgRNA_sum), to = range(data_plot$protein_smoothed,na.rm = TRUE)), name = "Protein")
        ) +
        # scale_color_manual(values = c("sgRNA counts" = "#DDA0DD", "Protein" = "#228B22")) +
        labs(
          x = "Cells ranked by sgRNA counts",
          color = NULL,  
          title = paste(sgRNA_select," - Protein (Window =", window_size, ")")
        ) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
          axis.title.y.left = element_text(color = "#663399", size = 12),
          axis.title.y.right = element_text(color = "#2E5A87", size = 12),  
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "bottom",  
          legend.text = element_text(size = 10)
        )
    data_plot_cellnum <- data.frame(
        mRNA = (data %>% filter(rna>0) %>% nrow())/nrow(data),
        Protein = (data %>% filter(protein>0) %>% nrow())/nrow(data)
    ) %>% 
        pivot_longer(cols = names(.)) %>% 
        mutate(
            name = factor(name,levels = c('mRNA','Protein'))
        )
    p3 <- ggplot(data = data_plot_cellnum,aes(x = name,y = value,fill = name)) +
        geom_bar(stat = 'identity',show.legend = FALSE) +
        scale_fill_manual(values = c('#A90C38','#2E5A87')) +
        scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
        theme_bw() +
        theme(
            plot.background = element_blank(),
            panel.grid = element_blank(),
            axis.line.y = element_line(),
            axis.line.x = element_line(),
            axis.text.y = element_text(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_blank()
        )
    data_res <- list(p1,p2,p3)
    return(data_res)
}

    
### ITGB1-sg1
sgRNA_select <- "ITGB1-sg1"
rna_select <- rownames(crispr_obj@assays$RNA)[grepl(rownames(crispr_obj@assays$RNA),pattern = 'ITGB1',ignore.case = TRUE)]
protein_select <- rownames(crispr_obj@assays$Protein)[grepl(rownames(crispr_obj@assays$Protein),pattern = 'CD29',ignore.case = TRUE)]
#
sgRNA_select
rna_select
protein_select
rna_select <- rna_select[1]
rna_select

data_res <- data_plot(
    crispr_obj = crispr_obj,
    sgRNA_identity_select = sgRNA_select,
    sgRNA_select = sgRNA_select,
    protein_select = protein_select,
    rna_select = rna_select,
    window_size = 21
)
data_res[[1]]/data_res[[2]]
ggsave("ITGB1_sg1_mRNA.pdf", plot = data_res[[1]], width = 7, height = 5, dpi = 300)
ggsave("ITGB1_sg1_Protein.pdf", plot = data_res[[2]], width = 7, height = 5, dpi = 300)


### GATA1-sg3
sgRNA_select <- "GATA1-sg3"
rna_select <- rownames(crispr_obj@assays$RNA)[grepl(rownames(crispr_obj@assays$RNA),pattern = 'GATA1',ignore.case = TRUE)]
protein_select <- rownames(crispr_obj@assays$Protein)[grepl(rownames(crispr_obj@assays$Protein),pattern = 'GATA1',ignore.case = TRUE)]
#
sgRNA_select
rna_select
protein_select

data_res <- data_plot(
    crispr_obj = crispr_obj,
    sgRNA_identity_select = sgRNA_select,
    sgRNA_select = sgRNA_select,
    protein_select = protein_select,
    rna_select = rna_select,
    window_size = 21
)
data_res[[1]]/data_res[[2]]
ggsave("GATA1_sg3_mRNA.pdf", plot = data_res[[1]], width = 7, height = 5, dpi = 300)
ggsave("GATA1_sg3_Protein.pdf", plot = data_res[[2]], width = 7, height = 5, dpi = 300)





##### SupplementaryFigure8 b
# Load packages.
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(purrr)
library(tidyr)
library(zoo) 
library(tibble)

setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNAcounts_mRNA_sgRNAcounts_Protein")
# Load object.
crispr_obj <- readRDS("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")

DefaultAssay(crispr_obj) <- "Protein"
crispr_obj <- subset(crispr_obj,subset = nCount_Protein<150)


### custom function
crispr_obj@meta.data <- crispr_obj@meta.data %>% 
    mutate(
        cell_type = case_when(
            sgRNA_identity %in% c('AAVS','NT1','NT2') ~ 'NT',
            TRUE ~ sgRNA_identity
        )
    )
crispr_obj$cell_type %>% table()

fun_p_value_singleSample <- function(crispr_obj = crispr_obj,cell_id_list = cell_id_list,sgRNA_select,protein_select,rna_select){
    data_res <- lapply(seq(10,100,10),FUN = function(cell_num){
        # downsample cells
        cell_id <- cell_id_list[[as.character(cell_num)]]
        subset_obj_cellnumber <- subset(subset_obj, cells = cell_id)
        # Extract the normalized expression data of Protein and RNA
        DefaultAssay(subset_obj_cellnumber) <- "Protein"
        protein_expression <- FetchData(object = subset_obj_cellnumber,vars = c(protein_select,'cell_type')) %>% 
            rename(Protein = !!sym(protein_select))
        DefaultAssay(subset_obj_cellnumber) <- "RNA"
        rna_expression <- FetchData(object = subset_obj_cellnumber,vars = c(rna_select,'cell_type')) %>% 
            rename(RNA = !!sym(rna_select))
        ## Integrate the expression matrix to calculate P-values
        data_exp <- protein_expression %>% 
            rownames_to_column('sample_id') %>% 
            left_join(rna_expression %>% rownames_to_column('sample_id'),by = c('sample_id','cell_type')) %>% 
            column_to_rownames('sample_id')
        result_pvalues <- map_dfr(sgRNA_select[1:3], function(target_sg) {
          # Extract the data for the target sg and NT
          df_sub <- data_exp %>%
            filter(cell_type %in% c(target_sg, 'NT')) %>%
            mutate(group = ifelse(cell_type == 'NT', 'control', 'target'))
          # Ensure that both groups have data
          if(length(unique(df_sub$group)) < 2) {
            return(tibble(
              sgRNA_identity = target_sg,
              p_raw_protein = NA,
              p_raw_RNA = NA
            ))
          }
          # Perform the Wilcoxon test for Protein and RNA
          p_protein <- tryCatch({
            wilcox.test(Protein ~ group, data = df_sub, exact = FALSE)$p.value %>% log10()
          }, error = function(e) NA)
          p_rna <- tryCatch({
            wilcox.test(RNA ~ group, data = df_sub, exact = FALSE)$p.value %>% log10()
          }, error = function(e) NA)
          tibble(
            cell_type = target_sg,
            p_raw_protein = p_protein,
            p_raw_RNA = p_rna
          )
        }) %>% 
            mutate(cell_num = cell_num)
    }) %>% do.call(rbind,.)
    return(data_res)
}

fun_p_value <- function(crispr_obj = crispr_obj,cell_id_all,sgRNA_select,protein_select,rna_select){
    sample_random = length(cell_id_all)
    data_res <- lapply(1:sample_random,function(random_i){
        set.seed(as.integer(Sys.time()) + random_i)
        data_tmp <- fun_p_value_singleSample(
            crispr_obj = crispr_obj,
            cell_id_list = cell_id_all[[(random_i) %>% as.character()]],
            sgRNA_select,protein_select,rna_select
        ) %>% 
            mutate(random = random_i)
    }) %>% do.call(rbind,.)
    return(data_res)
}


### GATA1
sample_random <- 10
sgRNA_select = c('GATA1-sg1','GATA1-sg2','GATA1-sg3','NT')
protein_select = 'GATA1-pAbO'
rna_select = 'GATA1'

# Select specific molecules
subset_obj <- subset(crispr_obj, subset = cell_type %in% c(sgRNA_select))
meta_data <- subset_obj@meta.data 
cell_id_all <- list()
for(i in 1:sample_random){
    cell_id_list <- lapply(seq(10,100,10),FUN = function(cell_num){
        # downsample cells
        cell_id <- meta_data %>%
          rownames_to_column('sample_id') %>%
          group_by(cell_type) %>%
          group_modify(~ {
            if(nrow(.x) >= cell_num){
              slice_sample(.x, n = cell_num)
            } else {
              .x
            }
          }) %>%
          ungroup() %>%
          pull(sample_id)
        return(cell_id)
    })
    names(cell_id_list) <- (seq(10,100,10) %>% as.character())
    cell_id_all[[i]] <- cell_id_list
}
names(cell_id_all) <- (1:sample_random  %>% as.character())

data_res <- fun_p_value(
    crispr_obj = crispr_obj,
    cell_id_all = cell_id_all,
    sgRNA_select = sgRNA_select,
    protein_select = protein_select,
    rna_select = rna_select
)

data_plot <- data_res %>% 
    pivot_longer(cols = c('p_raw_protein','p_raw_RNA'),names_to = 'Type',values_to = 'Expression') %>% 
    mutate(
        cell_num = factor(cell_num,levels = cell_num %>% unique()),
        Type = ifelse(grepl('protein',Type),yes = 'Protein',no = 'mRNA'),
        group = paste(Type,cell_type,sep = '-'),
        group = factor(group,levels = group %>% unique())
    )

options(repr.plot.width = 22,repr.plot.height = 8)
ggplot(data = data_plot,aes(x = cell_num,y = Expression,color = group)) +
    geom_boxplot(position = position_dodge(width = 0.7),outlier.shape = NA)+
    geom_point(aes(shape = Type, group = group),position = position_dodge(width = 0.7),size = 3) +
    # geom_dotplot(
    #     binaxis = "y", stackdir = "center",
    #     position = position_dodge(1)) +
    scale_color_manual(values = c('#209874','#6AC1A6','#D05E15','#EF955D','#6F6BAA','#9D9BC5')) +
    scale_shape_manual(values = c(1,16)) +
    guides(shape = guide_legend(title = 'Type',order = 1,override.aes = list(size = 5)),color = guide_legend(order = 2)) +
    labs(x = 'Number of Cells',y = 'Log10(P-Value)',color = 'Group') +
    theme(
        plot.title = element_text(hjust = 0.5,size = 32),
        panel.background = element_blank(),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.position = c(0.1,0.4),
        axis.line = element_line(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 24)
    )
ggsave("GATA1.pdf", plot = last_plot(), width = 18, height = 10, units = "in", dpi = 300)


### ITGB1
sample_random <- 10
sgRNA_select = c('ITGB1-sg1','ITGB1-sg2','ITGB1-sg3','NT')
protein_select = 'CD29-pAbO'
rna_select = 'ITGB1'

# Select specific molecules
subset_obj <- subset(crispr_obj, subset = cell_type %in% c(sgRNA_select))
meta_data <- subset_obj@meta.data 
cell_id_all <- list()
for(i in 1:sample_random){
    cell_id_list <- lapply(seq(10,100,10),FUN = function(cell_num){
        # downsample cells
        cell_id <- meta_data %>%
          rownames_to_column('sample_id') %>%
          group_by(cell_type) %>%
          group_modify(~ {
            if(nrow(.x) >= cell_num){
              slice_sample(.x, n = cell_num)
            } else {
              .x
            }
          }) %>%
          ungroup() %>%
          pull(sample_id)
        return(cell_id)
    })
    names(cell_id_list) <- (seq(10,100,10) %>% as.character())
    cell_id_all[[i]] <- cell_id_list
}
names(cell_id_all) <- (1:sample_random  %>% as.character())

data_res <- fun_p_value(
    crispr_obj = crispr_obj,
    cell_id_all = cell_id_all,
    sgRNA_select = sgRNA_select,
    protein_select = protein_select,
    rna_select = rna_select
)

data_plot <- data_res %>% 
    pivot_longer(cols = c('p_raw_protein','p_raw_RNA'),names_to = 'Type',values_to = 'Expression') %>% 
    mutate(
        cell_num = factor(cell_num,levels = cell_num %>% unique()),
        Type = ifelse(grepl('protein',Type),yes = 'Protein',no = 'mRNA'),
        group = paste(Type,cell_type,sep = '-'),
        group = factor(group,levels = group %>% unique())
    )

options(repr.plot.width = 22,repr.plot.height = 8)
ggplot(data = data_plot,aes(x = cell_num,y = Expression,color = group)) +
    geom_boxplot(position = position_dodge(width = 0.7),outlier.shape = NA)+
    geom_point(aes(shape = Type, group = group),position = position_dodge(width = 0.7),size = 3) +
    scale_color_manual(values = c('#209874','#6AC1A6','#D05E15','#EF955D','#6F6BAA','#9D9BC5')) +
    scale_shape_manual(values = c(1,16)) +
    guides(shape = guide_legend(title = 'Type',order = 1,override.aes = list(size = 5)),color = guide_legend(order = 2)) +
    labs(x = 'Number of Cells',y = 'Log10(P-Value)',color = 'Group') +
    theme(
        plot.title = element_text(hjust = 0.5,size = 32),
        panel.background = element_blank(),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.position = c(0.1,0.4),
        axis.line = element_line(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 24)
    )
ggsave("ITGB1.pdf", plot = last_plot(), width = 18, height = 10, units = "in", dpi = 300)





##### SupplementaryFigure8 c
### DEG_upset
library(dplyr)
library(tidyr)
library(UpSetR)
library(ggplot2)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/DEG_upset/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/DEG_upset/")

# all_positive
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/GATA1_NT_go_kegg/GATA1_NT_DEGs_positive.tsv"
deg <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_positive <- deg$gene_ID
# all_negative
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/GATA1_NT_go_kegg/GATA1_NT_DEGs_negative.tsv"
deg <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_negative <- deg$gene_ID
# all_both
all_both <- union(all_positive, all_negative)

# removetop20_positive
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/go_kegg_mixscale_GATA1_top20middlebottoom20/GATA1_NT_DEGs_positive.tsv"
deg <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
removetop20_positive <- deg$gene_ID
# removetop20_negative
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/go_kegg_mixscale_GATA1_top20middlebottoom20/GATA1_NT_DEGs_negative.tsv"
deg <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
removetop20_negative <- deg$gene_ID
# removetop20_both
removetop20_both <- union(removetop20_positive, removetop20_negative)

# removetop30_positive
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/go_kegg_mixscale_GATA1_top30middlebottoom30/GATA1_NT_DEGs_positive.tsv"
deg <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
removetop30_positive <- deg$gene_ID
# removetop30_negative
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/go_kegg_mixscale_GATA1_top30middlebottoom30/GATA1_NT_DEGs_negative.tsv"
deg <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
removetop30_negative <- deg$gene_ID
# removetop30_both
removetop30_both <- union(removetop30_positive, removetop30_negative)

gene_sets <- list(
  All_positive = all_positive,
  All_negative = all_negative,
  All_both = all_both,
  RemoveTop20_positive = removetop20_positive,
  RemoveTop20_negative = removetop20_negative,
  RemoveTop20_both = removetop20_both,
  RemoveTop30_positive = removetop30_positive,
  RemoveTop30_negative = removetop30_negative,
  RemoveTop30_both = removetop30_both
)

upset_input <- fromList(gene_sets)


pdf("DEG_UpSet.pdf", width = 8, height = 6)
upset(upset_input,
      sets = c(
        "RemoveTop30_negative","RemoveTop30_positive",
        "RemoveTop20_negative", "RemoveTop20_positive",
        "All_negative", "All_positive"
      ),
      keep.order = TRUE,
      order.by = "freq",
      sets.bar.color = "skyblue",
      main.bar.color = "steelblue",
      text.scale = c(1.5, 1.2, 1.5, 1.2, 1.2, 1.5))
dev.off()





##### SupplementaryFigure8 d
### go_plot
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GO/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GO/")



### positive

# All
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/GATA1_NT_go_kegg/GATA1_NT_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "positive regulation of mononuclear cell proliferation",
  "regulation of mononuclear cell proliferation",
  "positive regulation of leukocyte proliferation"
)
item <- go_data %>% filter(Cluster == "positive" & Description %in% target_descriptions)
item$group <- "All"
go_positive <- item

# remove_top20
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/go_kegg_mixscale_GATA1_top20middlebottoom20/mixscale_GATA1_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "positive regulation of mononuclear cell proliferation",
  "positive regulation of leukocyte proliferation",
  "regulation of mononuclear cell proliferation"
)
item <- go_data %>% filter(Cluster == "positive" & Description %in% target_descriptions)
item$group <- "remove_top20"
go_positive <- rbind(go_positive, item)

# remove_top30
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/go_kegg_mixscale_GATA1_top30middlebottoom30/mixscale_GATA1_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "positive regulation of mononuclear cell proliferation",
  "positive regulation of leukocyte proliferation",
  "regulation of mononuclear cell proliferation"
)
item <- go_data %>% filter(Cluster == "positive" & Description %in% target_descriptions)
item$group <- "remove_top30"
go_positive <- rbind(go_positive, item)

write.table(go_positive,
            file = "go_positive_3group.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# UpregulatedGOtermsinGATA1sg_3group
go_positive$group <- factor(go_positive$group, levels = c("All", "remove_top20", "remove_top30"))
go_positive <- go_positive[order(go_positive$Description, go_positive$group), ]

go_positive$y_label <- factor(seq_len(nrow(go_positive)), levels = rev(seq_len(nrow(go_positive))))

group_colors <- c(
  "All" = "#FAD4D4",         
  "remove_top20" = "#F08080", 
  "remove_top30" = "#DC143C"
)

go_positive$label_combined <- paste0(
  go_positive$Description,
  str_dup(" ", 5),
  "p.adjust=", signif(go_positive$p.adjust, 2)
)

ggplot(go_positive, aes(x = Count, y = y_label, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label_combined),
            x = 0.04, hjust = 0, vjust = 0.5, size = 4, color = "black") +
  geom_text(aes(label = geneID),
            x = 0.04, hjust = 0, vjust = 4, size = 3.5,
            color = "#EF3B2C") +
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
  labs(title = "Upregulated GO terms in GATA1sg", x = "Count", y = "") +
  coord_cartesian(clip = "off", expand = TRUE)
ggsave(filename = "UpregulatedGOtermsinGATA1sg_3group.pdf", plot = last_plot(), width = 10, height = 8.3, dpi = 300)
                                  