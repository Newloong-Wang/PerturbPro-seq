library(Seurat)
library(dplyr)
library(patchwork)
library(limma)
library(ggplot2)
library(paletteer)
library(viridis)

setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_before/")



crispr<-Read10X("/public/home/wangxinlong/Project/crispr1/data/20250102ABFC20221434-76_K562/sc20241231K562_addsequence/Ab-1_L2_S2401B2401_DBEC_MolsPerCell_MEX/")
crispr_obj<-CreateSeuratObject(counts =crispr$`Gene Expression`)
crispr_obj[['Protein']] = CreateAssayObject(counts = crispr$`Antibody Capture`)
cell_id <- colnames(crispr_obj)



sg_counts <- read.csv("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/sg-3_L2_Q0258W0072_umi_counts.csv")


sg_counts <- sg_counts[, !(names(sg_counts) %in% "barcode")]
# Keep only the cells present in cell_id
sg_counts <- sg_counts[sg_counts$id %in% cell_id, ]

all_sgRNA <- unique(sg_counts$sgRNA)
# Create a data frame containing all possible IDs and sgRNAs
complete_sgRNA_cell <- expand.grid(id = cell_id, sgRNA = all_sgRNA)
sg_counts <- merge(complete_sgRNA_cell, sg_counts, by = c("id", "sgRNA"), all.x = TRUE)
sg_counts$counts[is.na(sg_counts$counts)] <- 0
sg_counts <- reshape2::dcast(sg_counts, sgRNA ~ id, value.var = "counts", fill = 0)
rownames(sg_counts) <- sg_counts$sgRNA
sg_counts <- sg_counts[, -1]
sg_type_counts <- as.matrix(sg_counts)



###### Count the number of different sgRNAs present in each cell
sgRNA_count_per_cell <- apply(sg_type_counts, 2, function(x) sum(x > 0))
data <- data.frame(sgRNA_Count = sgRNA_count_per_cell)

summary_data <- data %>%
  group_by(sgRNA_Count) %>%
  summarise(Count = n())

ggplot(summary_data, aes(x = factor(sgRNA_Count), y = Count, fill = factor(sgRNA_Count))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  scale_fill_manual(values = c("#1f77b4",  
                               "#ff7f0e",  
                               "#2ca02c",  
                               "#d62728",  
                               "#9467bd",  
                               "#8c564b",  
                               "#e377c2",  
                               "#7f7f7f",  
                               "#bcbd22",
                               "#17becf",
                               "#fdae6b",
                               "#6baed6",
                               "#bcbddc")) + 
  labs(y = "Cell count", x = "Number of sgRNA types", fill = "sgRNA Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"))
ggsave("CellcontainssgRNAtypes.pdf", plot = last_plot(), width = 8, height = 8, units = "in", dpi = 300)


# sgRNA_Count 1 2 ≥3
summary_data <- summary_data %>%
  mutate(Group = factor(case_when(
    sgRNA_Count == 1 ~ "1",
    sgRNA_Count == 2 ~ "2",
    sgRNA_Count >= 3 ~ "≥3"
  ), levels = c("1", "2", "≥3")))

summary_data <- summary_data %>%
  filter(!is.na(Group))

grouped_data <- summary_data %>%
  group_by(Group) %>%
  summarise(Count = sum(Count))

ggplot(grouped_data, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  scale_fill_manual(values = c("#ff7f0e",
                               "#2ca02c",
                               "#1f77b4")) +
  labs(y = "Cell count", x = "Number of sgRNA types", fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"))
ggsave("CellcontainssgRNAtypes_123_filterbefore.jpg", plot = last_plot(), width = 4, height = 5, units = "in", dpi = 300)



##### sgRNA box plot
sg_type_counts_df <- as.data.frame(sg_type_counts)
sg_type_counts_df$sgRNA <- rownames(sg_type_counts_df)
long_data <- reshape2::melt(sg_type_counts_df, 
                  id.vars = "sgRNA",      
                  variable.name = "Cell_ID", 
                  value.name = "Counts") 
long_data <- long_data[long_data$Counts > 0, ]
long_data$sgRNA <- factor(long_data$sgRNA, 
                          levels = c("NT1", "NT2", "AAVS", 
                                     "ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3",
                                     "RPS6-sg1", "RPS6-sg2", "RPS6-sg3",
                                     "GATA1-sg1", "GATA1-sg2", "GATA1-sg3"))

ggplot(long_data, aes(x = sgRNA, y = Counts)) +
  geom_boxplot(fill = "#C8D7EB", color = "black", outlier.color = "black", outlier.shape = 1) +
  stat_summary(
    fun = median,  
    geom = "text",  
    aes(label = ..y..),  
    vjust = -0.5,  
    color = "blue"  
  ) +
  labs(x = "sgRNA Type", y = "sgRNA Counts") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)  
  ) 
ggsave("sgRNAboxplot_filterbefore.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 500)

# sgRNA box plot Y-axis 0-150
ggplot(long_data, aes(x = sgRNA, y = Counts)) +
  geom_boxplot(fill = "#C8D7EB", color = "black", outlier.color = "black", outlier.shape = 1) +
  stat_summary(
    fun = median,  
    geom = "text",  
    aes(label = ..y..),  
    vjust = -0.5,  
    color = "blue" 
  ) +
  labs(x = "sgRNA Type", y = "sgRNA Counts") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)  
  ) +
  coord_cartesian(ylim = c(0, 150))
ggsave("sgRNAboxplot_filterbefore0_150.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 500)



##### Cells with dominant sgRNA detected
sg_type_counts_df <- as.data.frame(sg_type_counts)

# Calculate the maximum value of each column and add it as a new row
top_sgRNA_number <- apply(sg_type_counts_df, 2, max)
sg_type_counts_df["top_sgRNA_number",] <- top_sgRNA_number
# Calculate the sum of each column and add it as a new row
sum_values <- colSums(sg_type_counts_df[1:12,])
sg_type_counts_df["sum_sgRNA_number", ] <- sum_values
# Divide the maximum value of each column by its sum and add the result as a new row
sg_type_counts_df['ratio_top_sum', ] <- sg_type_counts_df['top_sgRNA_number', ] / sg_type_counts_df['sum_sgRNA_number', ]

df_transposed <- as.data.frame(t(sg_type_counts_df))

ggplot(df_transposed, aes(x = sum_sgRNA_number, y = top_sgRNA_number, color = ratio_top_sum)) +
  geom_point(size = 1) +
  scale_color_viridis_c(option = "viridis") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Total sgRNA counts in cell",
    y = "Counts top sgRNA in cell",
    title = "Cells with dominant sgRNA detected",
    color = "ratio_top_sum"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  ) +
  coord_fixed(ratio = 1)
ggsave("Cells with dominant sgRNA detected.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 500)



##### Plot the cumulative density distribution
sgRNA_list <- rownames(sg_counts)

# Iterate over each sgRNA
for (sgRNA in sgRNA_list) {
  sgRNA_values <- sg_type_counts_df[sgRNA, ]
  sgRNA_df <- data.frame(value = as.numeric(sgRNA_values))
  p <- ggplot(sgRNA_df, aes(x = value)) +
    stat_ecdf(geom = "step", color = "orange") +
    labs(
      title = paste(sgRNA, "Cumulative Density Distribution"),
      x = "Value",
      y = "Cumulative Density"
    ) +
    theme_minimal()
  ggsave(
    filename = paste0(sgRNA, "_cumulativedensitydistribution.pdf"),
    plot = p, width = 5, height = 5, units = "in", dpi = 500
  )
  print(paste("Processed:", sgRNA))
}

# Iterate over each sgRNA
for (sgRNA in sgRNA_list) {
  sgRNA_values <- sg_type_counts_df[sgRNA, ]
  sgRNA_df <- data.frame(value = as.numeric(sgRNA_values))
  # X-axis 0-200
  p <- ggplot(sgRNA_df, aes(x = value)) +
    stat_ecdf(geom = "step", color = "orange") +
    scale_x_continuous(limits = c(0, 200)) +
    labs(
      title = paste(sgRNA, "Cumulative Density Distribution"),
      x = "Value",
      y = "Cumulative Density"
    ) +
    theme_minimal()
  ggsave(
    filename = paste0(sgRNA, "_cumulativedensitydistribution_0-200.pdf"),
    plot = p, width = 5, height = 5, units = "in", dpi = 500
  )
  print(paste("Processed:", sgRNA))
}



##### Cumulative density plot of all sgRNAs
all_data_df <- data.frame()

for (sgRNA in sgRNA_list) {
  sgRNA_values <- sg_type_counts_df[sgRNA, ]
  sgRNA_df <- data.frame(value = as.numeric(sgRNA_values), type = sgRNA)
  all_data_df <- rbind(all_data_df, sgRNA_df)
}

ggplot(all_data_df, aes(x = value, color = type)) +
  stat_ecdf(geom = "step") +
  labs(
    title = "Cumulative Density Distribution of All sgRNAs",
    x = "Value",
    y = "Cumulative Density",
    color = "sgRNA Type"
  ) +
  theme_minimal()
ggsave("All_sgRNA_Cumulative_Density_Distribution.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)

# Cumulative density plot of all sgRNAs X-axis 0-300
ggplot(all_data_df, aes(x = value, color = type)) +
  stat_ecdf(geom = "step") +
  scale_x_continuous(limits = c(0, 300)) +
  labs(
    title = "Cumulative Density Distribution of All sgRNAs (0-300 Zoomed)",
    x = "Value",
    y = "Cumulative Density",
    color = "sgRNA Type"
  ) +
  theme_minimal()
ggsave("All_sgRNA_Cumulative_Density_Distribution_X0-300.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)



##### Subtract the inflection point value of the cumulative density plot
dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_after/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_after/")

sg_type_counts <- as.matrix(sg_counts)
sg_type_counts_df <- as.data.frame(sg_type_counts)

subtract_values <- c(
  "AAVS" = 90,
  "GATA1-sg1" = 11,
  "GATA1-sg2" = 60,
  "GATA1-sg3" = 60,
  "ITGB1-sg1" = 90,
  "ITGB1-sg2" = 180,
  "ITGB1-sg3" = 80,
  "NT1" = 150,
  "NT2" = 150,
  "RPS6-sg1" = 11,
  "RPS6-sg2" = 11,
  "RPS6-sg3" = 11
)

sgTypeCounts_df_filter <- sg_type_counts_df
for (row_name in rownames(sg_type_counts_df)) {
  sgTypeCounts_df_filter[row_name, ] <- sgTypeCounts_df_filter[row_name, ] - subtract_values[row_name]
}
sgTypeCounts_df_filter[sgTypeCounts_df_filter < 0] <- 0



##### Transfer the sgRNA count data into the assay of the Seurat object 
sg_matrix <- as.matrix(sgTypeCounts_df_filter)
sgRNA_assay <- CreateAssayObject(counts = sg_matrix)
crispr_obj[["sgRNA"]] <- sgRNA_assay

saveRDS(crispr_obj, file = "crispr_obj.rds")



###### After subtracting the inflection point value from the cumulative density plot, count the number of different sgRNAs introduced into each cell
sgRNA_count_per_cell <- apply(sgTypeCounts_df_filter, 2, function(x) sum(x > 0))

sgRNA_group <- data.frame(
  Group = factor(case_when(
    sgRNA_count_per_cell == 1 ~ "1",
    sgRNA_count_per_cell == 2 ~ "2",
    sgRNA_count_per_cell >= 3 ~ "≥3"
  ), levels = c("1", "2", "≥3"))
)

summary_group <- sgRNA_group %>%
  group_by(Group) %>%
  summarise(Count = n())

summary_group <- summary_group %>%
  filter(!is.na(Group))

ggplot(summary_group, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5) +
  scale_fill_manual(values = c("#ff7f0e",  
                               "#2ca02c",  
                               "#1f77b4")) + 
  labs(y = "Cell count", x = "Number of sgRNA types", fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(face = "bold"))
ggsave("CellcontainssgRNAtypes_123_filterafter.jpg", plot = last_plot(), width = 4, height = 5, units = "in", dpi = 300)



#####  After subtracting the inflection point value from the cumulative density plot, sgRNA box plot
sgTypeCounts_df_filter$sgRNA <- rownames(sgTypeCounts_df_filter)
long_data <- reshape2::melt(sgTypeCounts_df_filter, 
                  id.vars = "sgRNA",      
                  variable.name = "Cell_ID", 
                  value.name = "Counts")
long_data <- long_data[long_data$Counts > 0, ]
long_data$sgRNA <- factor(long_data$sgRNA, 
                          levels = c("NT1", "NT2", "AAVS", 
                                     "ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3",
                                     "RPS6-sg1", "RPS6-sg2", "RPS6-sg3",
                                     "GATA1-sg1", "GATA1-sg2", "GATA1-sg3"))

ggplot(long_data, aes(x = sgRNA, y = Counts)) +
  geom_boxplot(fill = "#C8D7EB", color = "black", outlier.color = "black", outlier.shape = 1) +
  stat_summary(
    fun = median,  
    geom = "text",
    aes(label = ..y..),
    vjust = -0.5, 
    color = "blue" 
  ) +
  labs(x = "sgRNA Type", y = "sgRNA Counts") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)  
  ) 
ggsave("sgRNAboxplot_filterafter.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 500)


# sgRNA box plot Y-axis 0-300
ggplot(long_data, aes(x = sgRNA, y = Counts)) +
  geom_boxplot(fill = "#C8D7EB", color = "black", outlier.color = "black", outlier.shape = 1) +
  stat_summary(
    fun = median, 
    geom = "text",  
    aes(label = ..y..),  
    vjust = -0.5,  
    color = "blue"
  ) +
  labs(x = "sgRNA Type", y = "sgRNA Counts") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  ) +
  coord_cartesian(ylim = c(0, 300))
ggsave("sgRNAboxplot_filterafter0_300.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 500)

sgTypeCounts_df_filter$sgRNA <- NULL



##### Cells with dominant sgRNA detected
# Calculate the maximum value of each column and add it as a new row
top_sgRNA_number <- apply(sgTypeCounts_df_filter, 2, max)
sgTypeCounts_df_filter["top_sgRNA_number",] <- top_sgRNA_number
# Calculate the sum of each column and add it as a new row
sum_values <- colSums(sgTypeCounts_df_filter[1:12,])
sgTypeCounts_df_filter["sum_sgRNA_number", ] <- sum_values
# Divide the maximum value of each column by its sum and add the result as a new row
sgTypeCounts_df_filter['ratio_top_sum', ] <- sgTypeCounts_df_filter['top_sgRNA_number', ] / sgTypeCounts_df_filter['sum_sgRNA_number', ]

nan_count <- sum(is.na(sgTypeCounts_df_filter["ratio_top_sum", ]))
nan_count

df_transposed <- as.data.frame(t(sgTypeCounts_df_filter))

ggplot(df_transposed, aes(x = sum_sgRNA_number, y = top_sgRNA_number, color = ratio_top_sum)) +
  geom_point(size = 1) +
  scale_color_viridis_c(option = "viridis") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Total sgRNA counts in cell",
    y = "Counts top sgRNA in cell",
    title = "Cells with dominant sgRNA detected",
    color = "ratio_top_sum"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6)
  ) +
  coord_fixed(ratio = 1)
ggsave("Cells with dominant sgRNA detected.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 500)



##### Histogram of ratio_top_sum
row_values <- as.numeric(sgTypeCounts_df_filter["ratio_top_sum", ])
row_values[is.na(row_values)] <- 0
row_values_df <- data.frame(Value = row_values)


# Histogram
ggplot(row_values_df, aes(x=Value)) +
  geom_histogram(binwidth=0.1, fill="skyblue", color="black") +
  labs(title="ratio_top_sumHistogram1", x="Value", y="number")  
ggsave("ratio_top_sumHistogram1.pdf", plot = last_plot(), width = 3, height = 3, units = "in", dpi = 300)


# Histogram
row_values_df <- row_values_df %>%
  mutate(Group = cut(Value, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                     include.lowest = TRUE,
                     right = FALSE))

levels(row_values_df$Group) <- c("[0, 0.1]", "(0.1, 0.2]", "(0.2, 0.3]", "(0.3, 0.4]", "(0.4, 0.5]", "(0.5, 0.6]", "(0.6, 0.7]", "(0.7, 0.8]", "(0.8, 0.9]", "(0.9, 1]")

ggplot(row_values_df, aes(x=Group)) +
  geom_bar(fill="skyblue", color="black") +
  labs(title="ratio_top_sumHistogram2", x="Value", y="number") +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("ratio_top_sumHistogram2.pdf", plot = last_plot(), width = 3, height = 3, units = "in", dpi = 300)


# Histogram
row_values_df <- row_values_df %>%
  mutate(Group = cut(Value, breaks = seq(0, 1, by = 0.01),
                     include.lowest = TRUE,
                     right = FALSE))

labels <- paste0("(", seq(0, 0.99, by = 0.01), ", ", seq(0.01, 1, by = 0.01), "]")
labels[1] <- "[0, 0.01]"

ggplot(row_values_df, aes(x=Group)) +
  geom_bar(fill="skyblue", color="black") +
  labs(title="ratio_top_sumHistogram3", x="Value", y="Frequency") +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("ratio_top_sumHistogram3.pdf", plot = last_plot(), width = 10, height = 5, units = "in", dpi = 300)



##### Retain cells with ratio_top_sum > 0.99
sgTypeCounts_df_filter

sgTypeCounts_df_filter["ratio_top_sum", ] <- ifelse(is.na(sgTypeCounts_df_filter["ratio_top_sum", ]), 0, sgTypeCounts_df_filter["ratio_top_sum", ])

sgTypeCounts_df_filter2 <- sgTypeCounts_df_filter[, sgTypeCounts_df_filter["ratio_top_sum", ] > 0.99]

dim(sgTypeCounts_df_filter2)

sgTypeCounts_df_filter2["sgRNA_identity", ] <- apply(sgTypeCounts_df_filter2[1:12, ], 2, function(col) rownames(sgTypeCounts_df_filter2)[which.max(col)])

write.csv(sgTypeCounts_df_filter2, file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_after/cell_sgRNA_identity.csv", row.names = TRUE, quote = FALSE)
