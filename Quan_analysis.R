library(Seurat)
library(tidyverse)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
library(tibble)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
set.seed(42)
#Read data
setwd('/Users/lequanthinh/Desktop/GSE130973')
matrix_obj <- ReadMtx(mtx = 'GSE130973_matrix_filtered.mtx', cells= "GSE130973_barcodes_filtered.tsv", features="GSE130973_genes_filtered.tsv")
#Create Seurat object
skin_mtx<-CreateSeuratObject(count=matrix_obj, project = "skin aging", min.cells = 3,min.features = 200)
#Mitochondrial gene percentages
skin_mtx[["percent.mt"]]<-PercentageFeatureSet(skin_mtx,pattern = "^MT-")
#Visualize QC metrics
View(skin_mtx@meta.data)
VlnPlot(skin_mtx,features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
FeatureScatter(skin_mtx,feature1 = "nCount_RNA",feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
#Filter low quality cells
skin_mtx<-subset(skin_mtx,subset = nFeature_RNA>500 & nFeature_RNA <7500 & percent.mt <25)
#Normalize data
skin_mtx<-NormalizeData(skin_mtx)
#Identify highly variable genes
skin_mtx<-FindVariableFeatures(skin_mtx,selection.method = "vst", nfeatures = 2000)
top10<-head(VariableFeatures(skin_mtx),10)
plot1<-VariableFeaturePlot(skin_mtx)
LabelPoints(plot = plot1,points = top10, repel = TRUE)
#Scale data
skin_mtx<-ScaleData(skin_mtx)
#Principal component analysis (PCA) and Elbowplot
skin_mtx<-RunPCA(skin_mtx,features = VariableFeatures(object = skin_mtx))
DimHeatmap(skin_mtx,dims = 1, cells = 500,balanced = TRUE)
ElbowPlot(skin_mtx)
#Find cell cluster
skin_mtx<-FindNeighbors(skin_mtx, dims=1:20)
skin_mtx<-FindClusters(skin_mtx, resolution = 0.4)
DimPlot(skin_mtx,group.by="RNA_snn_res.0.4", label = TRUE)
#UMAP visualization
skin_mtx<-RunUMAP(skin_mtx,dims = 1:20)
DimPlot(skin_mtx, label = TRUE) + ggtitle("Unsupervised Cell Clusters")
#Identify cluster marker genes
markers<-FindAllMarkers(skin_mtx,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#Heatmap of cluster markers
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
DoHeatmap(skin_mtx,features = top10$gene)
#Identify cell types using known markers (either run manual or automatic anotation)
##Manual Anotation
marker_list <- list(
  Pericytes                = c("ACTA2", "RGS5", "PDGFRB"),
  Epidermal_stem_cells     = c("KRT5", "KRT14", "TP63", "ITGB1", "ITGA6"),
  Differentiated_keratinocytes = c("KRT1", "KRT10", "SBSN", "KRTDAP"),
  Fibroblasts               = c("PDGFRA", "LUM", "DCN", "VIM", "COL1A2"),
  Myeloid_cells             = c("AIF1", "LYZ", "HLA-DRA", "CD68", "ITGAX"),
  T_cells                    = c("CD3D", "CD3G", "CD3E", "LCK"),
  Vascular_ECs               = c("SELE", "CLDN5", "VWF", "CDH5"),
  Lymphatic_ECs               = c("PROX1", "CLDN5", "LYVE1"),
  Erythrocytes                = c("HBA1", "HBA2", "HBB"),
  Melanocytes                  = c("PMEL", "MLANA", "TYRP1", "DCT")
)
### check all genes actually exist in your object (avoids silent drops/errors)
all_genes <- unlist(marker_list)
missing_genes <- setdiff(all_genes, rownames(skin_mtx))
if (length(missing_genes) > 0) {
  message("Missing from object: ", paste(missing_genes, collapse = ", "))
}
for (ct in names(marker_list)) {
  genes_present <- intersect(marker_list[[ct]], rownames(skin_mtx))
  skin_mtx <- AddModuleScore(skin_mtx,
                               features = list(genes_present),
                               name = paste0(ct, "_score"))
}
score_cols <- grep("_score1$", colnames(skin_mtx@meta.data), value = TRUE)
# average score per cluster -> helps you assign objectively rather than by eye
score_summary <- skin_mtx@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean)) %>%
  as.data.frame()
print(score_summary)
# for each cluster, which signature scores highest
score_matrix <- score_summary[,-1]
rownames(score_matrix) <- score_summary$seurat_clusters
best_call <- apply(score_matrix, 1, function(x) names(x)[which.max(x)])
best_call <- sub("_score1$", "", best_call)
data.frame(cluster = rownames(score_matrix), top_signature = best_call)
cluster_annotations <- c("0" = "Fibroblasts","1" = "Fibroblasts","2" = "Myeloid_cells","3" = "Vascular_ECs","4" = "Fibroblasts","5" = "Fibroblasts",
  "6" = "Pericytes","7" = "T_cells","8" = "Differentiated_keratinocytes","9" = "Epidermal_stem_cells","10" = "Epidermal_stem_cells",
  "11" = "Myeloid_cells","12" = "Lymphatic_ECs","13" = "T_cells","14" = "Pericytes","15" = "Melanocytes","16" = "Pericytes")
skin_mtx$cluster_celltype <- unname(cluster_annotations[as.character(skin_mtx$seurat_clusters)])
Idents(skin_mtx) <- "cluster_celltype"
DimPlot(skin_mtx, reduction = "umap", group.by = "cluster_celltype",
        label = TRUE, repel = TRUE) + NoLegend()
#-----------------------------------------------------------------------------
#Add age metadata
skin_mtx$age <- ifelse(grepl("-1$", colnames(skin_mtx)) | grepl("-2$", colnames(skin_mtx)),"YOUNG","OLD")
#Compare cell type population
cell_counts <- table(skin_mtx$cluster_celltype, skin_mtx$age)
cell_df <- as.data.frame(cell_counts)
# Normalize by total cells per age group (calculate proportions)
cell_df <- cell_df %>%
  group_by(Var2) %>%
  mutate(Proportion = 100 * Freq / sum(Freq)) %>%
  ungroup()
# Plot normalized proportions
ggplot(cell_df, aes(x = Var2, y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(title = "Normalized Cell Type Composition (Old vs Young)", y = "Proportion")
#Global aging gene analysis
skin_mtx$sampleid <- dplyr::case_when(
  grepl("-1$", colnames(skin_mtx)) ~ "Y1",
  grepl("-2$", colnames(skin_mtx)) ~ "Y2",
  grepl("-3$", colnames(skin_mtx)) ~ "O1",
  grepl("-4$", colnames(skin_mtx)) ~ "O2",
  grepl("-5$", colnames(skin_mtx)) ~ "O3",
  TRUE ~ NA_character_
)
pseudo_mtx <- AggregateExpression(
  skin_mtx,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("age", "sampleid")
)
Cells(pseudo_mtx)
Idents(pseudo_mtx) <- "age"
aging_genes <- FindMarkers(
  pseudo_mtx,
  ident.1 = "OLD",
  ident.2 = "YOUNG",
  test.use = "DESeq2",
  min.cells.group = 2
)
head(aging_genes)
#Cell type specific aging genes
Idents(skin_mtx) <- "cluster_celltype"
celltypes <- levels(skin_mtx)
aging_by_celltype <- list()
for (ct in celltypes) {
  cells <- subset(skin_mtx, idents = ct)
  if(length(unique(cells$age)) < 2) next
  Idents(cells) <- "age"
  deg <- FindMarkers(
    cells,
    ident.1 = "OLD",
    ident.2 = "YOUNG")
  deg$celltype <- ct
  aging_by_celltype[[ct]] <- deg
}
aging_celltype_results <- do.call(rbind, aging_by_celltype)
#Volcanic plot of DEG in young vs old
## pick your cell type
my_celltype <- "Myeloid_cells"   # replace with whatever cell type we want
deg_subset <- aging_celltype_results %>%
  filter(celltype == my_celltype) %>%
  rownames_to_column("gene") %>%
  mutate(gene = sub("\\.\\.\\.\\d+$", "", gene))  # strip rowname suffixes like "MERTK...123" if rbind added them
## handle p_val_adj = 0
min_nonzero_p <- min(deg_subset$p_val_adj[deg_subset$p_val_adj > 0], na.rm = TRUE)
deg_subset$p_val_adj_plot <- ifelse(deg_subset$p_val_adj == 0, min_nonzero_p, deg_subset$p_val_adj)

deg_subset <- deg_subset %>%
  mutate(
    neglog10p = -log10(p_val_adj_plot),
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 1  ~ "Up in OLD",
      p_val_adj < 0.05 & avg_log2FC < -1 ~ "Up in YOUNG",
      TRUE ~ "NS"
    )
  )
top_genes <- deg_subset %>%
  filter(sig != "NS") %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 8) %>%
  pull(gene)
deg_subset$label <- ifelse(deg_subset$gene %in% top_genes, deg_subset$gene, NA)
ggplot(deg_subset, aes(x = avg_log2FC, y = neglog10p, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up in OLD" = "firebrick",
                                "Up in YOUNG" = "steelblue",
                                "NS" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_label_repel(aes(label = label), size = 3, na.rm = TRUE, show.legend = FALSE,
                   max.overlaps = 20) +
  theme_bw() +
  labs(title = paste("OLD vs YOUNG —", my_celltype),
       x = "log2FC (OLD/YOUNG)", y = "-log10(adj p-value)")
#Save results
saveRDS(skin_mtx, "/Users/lequanthinh/Desktop/GSE130973/Outputs/skin_analysis_new.rds")
write.csv(markers, "/Users/lequanthinh/Desktop/GSE130973/Outputs/cluster_markers.csv")
write.csv(aging_genes, "/Users/lequanthinh/Desktop/GSE130973/Outputs/aging_genes.csv")
write.csv(aging_celltype_results, "/Users/lequanthinh/Desktop/GSE130973/Outputs/aging_by_celltype.csv")
#subcluster based on signature score--------------------------------------------------
## 1. Subset
myeloids <- subset(skin_mtx, subset = cluster_celltype == "Myeloid_cells")
myeloids <- NormalizeData(myeloids)
## 2. Re-find variable features within just this subset
myeloids <- FindVariableFeatures(myeloids, nfeatures = 2000)
## 3. Re-scale (on the subset's own variable features)
myeloids <- ScaleData(myeloids)
## 4. Re-run PCA
myeloids <- RunPCA(myeloids, npcs = 30)
## Check how many PCs actually capture variance in this subset — 
# with fewer cells and less heterogeneity you often need fewer PCs than the full dataset
ElbowPlot(myeloids)
library(harmony)
myeloids<-RunHarmony(myeloids, group.by.vars = "sampleid", reduction.use = "pca",
                     reduction.save = "harmony")
## 5. Neighbors, clustering, UMAP — using harmony embedding if you ran it
myeloids <- FindNeighbors(myeloids, dims = 1:10, reduction = "harmony")
myeloids <- FindClusters(myeloids, resolution = 0.3)
myeloids <- RunUMAP(myeloids, dims = 1:10, reduction = "harmony")
## 6. Visualize
DimPlot(myeloids, label = TRUE)
Idents(myeloids)="seurat_clusters"
myeloids_markers<-FindAllMarkers(myeloids,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
myeloids_top10 <- myeloids_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
DoHeatmap(myeloids,features = myeloids_top10$gene)
#Subtype of myeloids
##SingleR automatic annotation using published dataset
library(SingleR)
library(celldex)
library(scran)
library(scuttle)
library(scRNAseq)
ref<-fetchDataset("he-organs-2020", "2023-12-21",path = "skin")
ref <- logNormCounts(ref)
expr<-as.SingleCellExperiment(myeloids)
colLabels(expr) <- myeloids$seurat_clusters
predictions<-SingleR(test = expr, ref = ref, labels = ref$reclustered.fine,clusters = colLabels(expr))
cl <- as.character(colLabels(expr))
myeloids$singler_type <- predictions$labels[match(cl, rownames(predictions))]
Idents(myeloids) <- "singler_type"
DimPlot(myeloids, reduction = "umap", group.by = "singler_type", label = TRUE, repel = TRUE)
singler_markers<-FindAllMarkers(myeloids,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#DEGs aging genes on myeloids subsets
subset <- levels(myeloids)
aging_by_subset <- list()
for (st in subset) {
  cells <- subset(myeloids, idents = st)
  if(length(unique(cells$age)) < 2) next
  Idents(cells) <- "age"
  deg <- FindMarkers(
    cells,
    ident.1 = "OLD",
    ident.2 = "YOUNG")
  deg$subset <- st
  aging_by_subset[[st]] <- deg
}
aging_subset_results <- do.call(rbind, aging_by_subset)
#Save
write.csv(myeloids_markers, "/Users/lequanthinh/Desktop/GSE130973/Outputs/myeloids_clusters.csv")
write.csv(aging_subset_results, "/Users/lequanthinh/Desktop/GSE130973/Outputs/myeloids_subset.csv")
##Heatmap of myeloid markers
myeloid_order <- c("cDC1", "cDC2", "Langerhans")
myeloid_type_markers <- c("CPVL","CLEC9A","WDFY4","XCR1","CLEC10A","CD1C","CD1A","FCER1A","CD207","EPCAM")
Idents(myeloids) <- factor(myeloids$singler_type, levels = myeloid_order)
myeloid_avg <- AverageExpression(myeloids, features = myeloid_type_markers,
                                 group.by = "ident", slot = "data")$RNA
myeloid_avg <- myeloid_avg[myeloid_type_markers, , drop = FALSE]
myeloid_mat <- t(scale(t(as.matrix(myeloid_avg))))
colnames(myeloid_mat) <- sub("^g", "", colnames(myeloid_mat))
myeloid_mat <- myeloid_mat[, colSums(!is.na(myeloid_mat)) > 0, drop = FALSE]
col_fun <- colorRamp2(c(min(myeloid_mat, na.rm = TRUE), 0, max(myeloid_mat, na.rm = TRUE)),
                      c("#2166AC", "white", "#B2182B"))
Heatmap(myeloid_mat,
        name = "Scaled\nexpression",
        col = col_fun,
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_title = "myeloid marker expression",
        row_names_gp = gpar(fontsize = 10),
        rect_gp = gpar(col = "white", lwd = 0.5))
#---------------------------------------------------------------------------------------------------
#Signatures check
##Plot cell-type specific signature score-----------------------------
Efferocytosis_list<-c("NR1H3", "SIRPB1", "MERTK", "CEBPB", "CPT1C", "CPT1A", "CPT1B",
                      "CREB1", "CRK", "CRKL", "MAPK14", "CX3CR1", "AGER", "DNMT3A", "DNMT3B", "DOCK1", "DUSP2",
                      "DUSP4", "DUSP5", "DUSP7", "DUSP8", "S1PR1", "ANO6", "ANO5", "EPO", "EPOR", "STAB1", "ATP11B", "ATP11A",
                      "VPS39", "VPS8", "SIRT1", "PLA2G15", "ALOX5", "ALOX15", "HAVCR1", "VPS41", "XKR6", "ATP11C", "GPR132",
                      "HIF1A", "RAB7B", "ANO9", "XKR7", "IL10", "ITGAV", "ITGB3", "ITGB5", "JAK2", "ARG1", "ARG2", "XKR9", "LIPA", "LRP1",
                      "ARNT", "NFATC1", "NFATC2", "NFATC3", "NFATC4", "ATP2A1", "ATP2A2", "ATP2A3", "ODC1", "P2RY2", "P2RY6", "ANO7",
                      "PBX1", "GULP1", "RAB14", "ATP8A2", "PPARD", "PPARG", "SLC66A1", "XKR8", "SIRPG", "STAB2", "TMEM30A", "AXL",
                      "VPS11", "MAPK1", "MAPK3", "MAPK11", "MAPK13", "MAP2K1", "MAP2K2", "PTGER2", "PTGER4", "PTGS2", "PTK2",
                      "ADGRB1", "VPS18", "RAB5A", "RAB5B", "RAB5C", "RAC1", "RXRA", "MAPK12", "ANO3", "RAB17",
                      "VPS16", "P2RY12", "VPS33A", "AD AM17", "TGFB1", "C1QB", "C1QC",
                      "TYRO3", "NR1H2", "UQCRFS1", "RAB7A", "DUSP16", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
                      "CASP7", "MEGF11", "MEGF10", "HAVCR2", "SCARF1", "ADAM9", "SPHK1", "CH25H", "TIMD4", "MAPKAPK2",
                      "TGFBRAP1", "CD36", "BCAR1", "ELMO1") 
                      #The gene list were from KEGG efferocytosis pathway, excluded don't eat me, glucose uptake elements and included elements from efferocytes only
skin_mtx<-AddModuleScore(skin_mtx,features = list(Efferocytosis_list),name = "Efferocytosis_Score")
Idents(skin_mtx)="cluster_celltype"
FeaturePlot(skin_mtx,features = "Efferocytosis_Score1", label = TRUE)
VlnPlot(skin_mtx, features = "Efferocytosis_Score1", group.by = "cluster_celltype",pt.size = 0) + 
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.9),outlier.size = 0.6, alpha = 0.7, linewidth = 0.3) +
  geom_hline(yintercept = c(0), linetype = "dashed", color = "red", linewidth = 0.4) +
  ggtitle("Efferocytosis score")
#Complex heatmap--------------------------------------------------------------------
celltype_order <- c("Myeloid-cells","Fibroblasts","Lymphatic-ECs","Vascular-ECs")
# Average expression per gene x celltype x age
avg <- AverageExpression(skin_mtx, features = Efferocytosis_list,
                         group.by = c("cluster_celltype", "age"),
                         assays = "RNA", slot = "data")$RNA
missing_genes <- setdiff(Efferocytosis_list, rownames(avg))
genes_present <- intersect(Efferocytosis_list, rownames(avg))
# Row z-score
mat <- t(scale(t(as.matrix(avg))))[genes_present, , drop = FALSE]
mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
# Reorder columns: celltype block, then Y then O within each
info <- data.frame(col = colnames(mat),
                   ct = sub("_(YOUNG|OLD)$","",colnames(mat)),
                   age = sub("^.*_","",colnames(mat))) %>%
  filter(ct %in% celltype_order) %>%                      # exclude any ct not in your defined order
  mutate(ct = factor(ct, levels = celltype_order[celltype_order %in% ct]),
         age = factor(age, levels = c("YOUNG","OLD"))) %>%
  arrange(ct, age)
mat <- mat[, info$col, drop = FALSE]
Heatmap(mat,
        name = "Gene expression",
        col = colorRamp2(c(-2,0,2), c("red","white","blue")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        column_split = info$ct, column_gap = unit(2,"mm"),
        column_labels = ifelse(info$age=="YOUNG","Y","O"),
        column_names_rot = 0, row_names_side = "right",
        row_names_gp = gpar(fontface = "italic"),
        rect_gp = gpar(col = "white", lwd = 0.5),
        heatmap_legend_param = list(direction = "horizontal"))
