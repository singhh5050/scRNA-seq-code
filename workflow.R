# ----
# Load packages
library(Seurat)
library(dplyr)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(data.table) #fread
library(Matrix) #sparse matrix
library(rhdf5) #to read in .h5 files
library(fgsea)
library(tibble)
library(tidyr)
library(patchwork)
library(pheatmap)
library(decoupleR)

# ----
# QC Preprocessing - Read in data

# Control samples
# Get list of all files
vector_of_files <- list.files(
  path = "/Users/hs/Desktop/breast-cancer-analysis/data/GSE164898_RAW",
  recursive = TRUE, 
  full.names = TRUE
)

# Read all of them in
load_data <- function(x) {
  x <- Read10X_h5(x)
  CreateSeuratObject(counts = x, min.cells = 3, min.features = 200)
}
ctrl_seurat_list <- lapply(vector_of_files, load_data)

# Merge all Seurats in the list together
ctrl_seurat <- Reduce(merge, ctrl_seurat_list)

# Assign unique IDs for each control sample
for (i in seq_along(ctrl_seurat_list)) {
  ctrl_seurat_list[[i]]@meta.data[["orig.ident"]] <- paste0("normal", i)
}

ctrl_seurat@meta.data$orig.ident <- "Normal"

# Experimental Samples

series <- fread("/Users/hs/Desktop/breast-cancer-analysis/data/GSE176078_series_matrix.txt", fill = TRUE)
series <- series[-c(31, 68, 70),]
series_subset <- t(series[c(31, 38, 40),])[-1, ]
series_TNBC <- as.vector(series_subset[series_subset[, 3] == "clinical_subtype: TNBC", 1])

# Read them in
exp_counts <- readMM("/Users/hs/Desktop/breast-cancer-analysis/data/Wu_etal_2021_BRCA_scRNASeq/matrix.mtx")
exp_features <- fread("/Users/hs/Desktop/breast-cancer-analysis/data/Wu_etal_2021_BRCA_scRNASeq/features.tsv.gz")
exp_barcodes <- fread("/Users/hs/Desktop/breast-cancer-analysis/data/Wu_etal_2021_BRCA_scRNASeq/barcodes.tsv.gz")

# For some reason the first element of the features.tsv and barcodes.tsv files becomes a column name. I want to use the column name as the first row, 
colname_to_row <- function(x) {
  x <- rbind(colnames(x), x, fill = TRUE)
  x[1, 2] <- x[1, 1]
  pull(x, var = 2)
}

exp_features <- colname_to_row(exp_features)
exp_barcodes <- colname_to_row(exp_barcodes)

exp_counts@Dimnames[[1]] <- exp_features
exp_counts@Dimnames[[2]] <- exp_barcodes

exp_seurat <- CreateSeuratObject(
  counts = exp_counts,
  min.cells = 3,
  min.features = 200,
  assay = "RNA"
)

exp_seurat <- subset(exp_seurat, idents = series_TNBC)
exp_seurat@meta.data$orig.ident <- "Cancer"

ctrlmetadata <- ctrl_seurat@meta.data
expmetadata <- exp_seurat@meta.data
metadata <- rbind(ctrlmetadata, expmetadata)

# Create a merged Seurat object
merged_seurat <- merge(
  x = ctrl_seurat, 
  y = exp_seurat, 
  add.cell.id = c("ctrl", "exp")
)

head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# ----
# QC - Generating Quality Metrics

# Calculate the # of genes detected per UMI (novelty score)  (higher =  more complex)
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Calculate the mitochondrial ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# scRNA-seq analysis: Assessing QC metrics
# Note: we aren't filtering out doublets because often it removes viable cells

# Visualize the # of UMI (transcripts) per cell
# Expected: at least 500 at the very least
countdep <- merged_seurat@meta.data %>% 
  ggplot(aes(color = orig.ident, x = nCount_RNA, fill = orig.ident, ..scaled..)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  xlab("UMIs (Transcripts) Detected Per Cell") +
  ylab("Cell Density") +
  labs(color = "Condition", fill = "Condition") +
  ggtitle("Count Depth\nVisualization") +
  theme(
    plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 500, linetype = "dashed") +
  geom_vline(xintercept = 50000, linetype = "dashed")

# Visualize the distribution of genes detected per cell
# High quality: one single peak
# Small shoulder means either - cells failed or different cell types
genedist <- merged_seurat@meta.data %>% 
  ggplot(aes(color = orig.ident, x = nFeature_RNA, fill = orig.ident, ..scaled..)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  xlab("Genes Detected Per Cell") +
  ylab("Cell Density") +
  labs(color = "Condition", fill = "Condition") +
  ggtitle("Gene Distribution\nVisualization") +
  theme(
    plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 300, linetype = "dashed") +
  geom_vline(xintercept = 6500, linetype = "dashed")

# Visualize the overall complexity of gene expression by visualizing novelty scores (# of genes / # of UMIs)
# Good quality score > 0.8
trancomp <- merged_seurat@meta.data %>% 
  ggplot(aes(color = orig.ident, x = log10GenesPerUMI, fill = orig.ident, ..scaled..)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  xlab("Novelty Score\n(Genes Detected / UMIs Detected)") +
  ylab("Cell Density") +
  labs(color = "Condition", fill = "Condition") +
  ggtitle("Transcriptome Complexity\nVisualization") +
  theme(
    plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0.8, linetype = "dashed")

# Visualize the distribution of mitochondrial gene expression detected per cell 
# Dead or dying cells > 0.2
mitoexp <- merged_seurat@meta.data %>% 
  ggplot(aes(color = orig.ident, x = mitoRatio, fill = orig.ident, ..scaled..)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  xlab("Mitochondrial Ratio") +
  ylab("Cell Density") +
  labs(color = "Condition", fill = "Condition") +
  ggtitle("Mitochondrial Gene\nExpression Visualization") +
  theme(
    plot.title = element_text(size = 15, face = "bold", vjust = 3, hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0.2, linetype = "dashed")

qcplots <- plot_grid(countdep, genedist, trancomp, mitoexp)

qctitle <- ggdraw() + 
  draw_label("   Quality Control Metrics", fontface = 'bold', size = 20, x = 0, hjust = 0.5) +
  theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(qctitle, qcplots, ncol = 1, rel_heights = c(0.1, 1))

# ----
# Execute QC
filtered_seurat <- subset(
  x = merged_seurat,
  subset = (nCount_RNA >= 500) &
    (nCount_RNA <= 50000) &
    (nFeature_RNA >= 300) &
    (nFeature_RNA <= 6500) &
    (log10GenesPerUMI >= 0.8) &
    (mitoRatio <= 0.2)
)

# Visualize the correlation between number of genes detected and number of UMIs, while overlaying mitochondrial ratio
# Jointly considering factors to determine filtering
# Look at before vs after
plot_grid(
  merged_seurat@meta.data %>% 
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black", limits = c(0, 1), name = "Mitochondrial \nRatio") +
    stat_smooth(method = lm) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("UMIs Detected Per Cell") +
    ylab("Genes Detected Per Cell") +
    theme_classic() +
    theme(plot.title = element_text(size = 15, face = "bold")) +
    geom_vline(xintercept = 500, linetype = "dashed") +
    geom_hline(yintercept = 300, linetype = "dashed") +
    geom_vline(xintercept = 50000, linetype = "dashed") +
    geom_hline(yintercept = 6500, linetype = "dashed") +
    geom_vline(xintercept = 150, linetype = "blank") +
    geom_hline(yintercept = 130, linetype = "blank") + 
    geom_vline(xintercept = 305000, linetype = "blank") +
    geom_hline(yintercept = 17500, linetype = "blank") +
    facet_wrap(~orig.ident, scales = 'free') +
    ggtitle("Unfiltered Data"),
  
  filtered_seurat@meta.data %>% 
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black", limits = c(0, 1), name = "Mitochondrial \nRatio") +
    stat_smooth(method = lm) +
    scale_x_log10() +
    scale_y_log10() +
    xlab("UMIs Detected Per Cell") +
    ylab("Genes Detected Per Cell") +
    theme_classic() +
    theme(plot.title = element_text(size = 15, face = "bold")) +
    geom_vline(xintercept = 500, linetype = "dashed") +
    geom_hline(yintercept = 300, linetype = "dashed") +
    geom_vline(xintercept = 50000, linetype = "dashed") +
    geom_hline(yintercept = 6500, linetype = "dashed") +
    geom_vline(xintercept = 150, linetype = "blank") +
    geom_hline(yintercept = 130, linetype = "blank") + 
    geom_vline(xintercept = 305000, linetype = "blank") +
    geom_hline(yintercept = 17500, linetype = "blank") +
    facet_wrap(~orig.ident, scales = 'free') +
    ggtitle("Filtered Data"),
  ncol = 1
)

# Gene-Level Filtering

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts (transcripts) per cell
nonzero <- counts > 0

# Let's keep only genes which are expressed in 10 or more cells

# Sums all TRUE values and returns TRUE only if there are 10 or more TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keep genes in 10 or more cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered Seurat
saveRDS(filtered_seurat, file = "/Users/hs/Desktop/breast-cancer-analysis/data/filtered_seurat.RDS")

# ----
# Normalization

# Split seurat object by condition so we can treat normalization/SCT separately
split_seurat <- SplitObject(filtered_seurat, split.by = "orig.ident")
split_seurat <- split_seurat[c("ctrl", "exp")]

# Adjust the limit for allowable object sizes within R
options(future.globals.maxSize = 4000 * 1024^2)

# Perform the sctransform on each sample
split_seurat[[1]] <- SCTransform(split_seurat[[1]], vars.to.regress = "mitoRatio")
split_seurat[[2]] <- SCTransform(split_seurat[[2]], vars.to.regress = "mitoRatio")

# Check which assays are stored in objects
split_seurat$ctrl@assays

# Save the split seurat object
save(split_seurat, file = "/Users/hs/Desktop/breast-cancer-analysis/data/split_seurat.RData")
saveRDS(split_seurat, "/Users/hs/Desktop/breast-cancer-analysis/data/split_seurat.rds")

# ----
# Dataset integration across different conditions

# Select the most variable features (first 3000) to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)

# Perform CCA, find best buddies! May take a while to run. Progress bar may stay at 0%, but it is actually running
integ_anchors <- FindIntegrationAnchors(
  object.list = split_seurat,
  normalization.method = "SCT",
  anchor.features = integ_features
)
saveRDS(integ_anchors, "/Users/hs/Desktop/breast-cancer-analysis/data/integ_anchors.rds")

# Integrate across conditions (Harmony - alternative)
seurat_integrated <- IntegrateData(
  anchorset = integ_anchors,
  normalization.method = "SCT"
)
saveRDS(seurat_integrated, "/Users/hs/Desktop/breast-cancer-analysis/data/seurat_integrated.rds")

# scRNA-seq analysis: PCA/UMAP dimensionality reduction
seurat_integrated@meta.data[["orig.ident"]] <- replace(
  seurat_integrated@meta.data[["orig.ident"]],
  seurat_integrated@meta.data[["orig.ident"]] == "ctrl",
  "Normal"
)

seurat_integrated@meta.data[["orig.ident"]] <- replace(
  seurat_integrated@meta.data[["orig.ident"]],
  seurat_integrated@meta.data[["orig.ident"]] == "exp",
  "Cancer"
)

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
pca_plot <- PCAPlot(seurat_integrated, group.by = "orig.ident") +
  labs(
    color = "Condition", 
    title = "Principal Component Analysis (PCA) Plot"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 18, margin = margin(b = 12), hjust = 0.5)
  )

# Explore heatmap of PCs
DimHeatmap(seurat_integrated, dims = 1:9, cells = 500, balanced = TRUE)

# Print out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], dims = 1:10, nfeatures = 5)

# Construct elbow plot
elbowplot <- ElbowPlot(object = seurat_integrated, ndims = 50)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40, reduction = "pca")

# Plot UMAP
DimPlot(seurat_integrated)

# Plot UMAP split by sample
umap_plot <- DimPlot(seurat_integrated, group.by = "orig.ident") +
  labs(
    color = "Condition", 
    title = "Uniform Manifold Approximation
Projection (UMAP) Plot"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 18, margin = margin(b = 12), hjust = 0.5)
  )

# Save integrated object
saveRDS(seurat_integrated, "/Users/hs/Desktop/breast-cancer-analysis/data/seurat_integrated.rds")

# ----
# Clustering
seurat_integrated@meta.data[["orig.ident"]] <- replace(
  seurat_integrated@meta.data[["orig.ident"]],
  seurat_integrated@meta.data[["orig.ident"]] == "ctrl",
  "Normal"
)

seurat_integrated@meta.data[["orig.ident"]] <- replace(
  seurat_integrated@meta.data[["orig.ident"]],
  seurat_integrated@meta.data[["orig.ident"]] == "exp",
  "Cancer"
)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:40)

# Determine the clusters for various resolutions
seurat_integrated <- FindClusters(object = seurat_integrated, resolution = c(1))

# Explore resolutions
seurat_integrated@meta.data %>% View()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1"

# Plot the UMAP
k_means_plot <- DimPlot(
  seurat_integrated,
  reduction = "umap",
  label = TRUE,
  label.size = 6
) +
  labs(title = "K-Nearest Neighbors
(KNN) Graph") + 
  theme(legend.position = "right") +
  theme(
    plot.title = element_text(face = "bold", size = 18, margin = margin(b = 12), hjust = 0.5)
  )

# ----
# Marker identification (database - https://panglaodb.se/)

# Endothelial cell markers
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("PECAM1", "EGFL7", "ID3", "MMRN1"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Clusters 7, 14, 17, 35, 37, & 42
endo <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("PECAM1"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# Epithelial cell markers
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("EPCAM", "KRT14", "CLDN1", "MUC1"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Clusters 2, 5, 12, 13, 20, 22, 25, 27, 28, 30, 31, & 38
epi <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("EPCAM"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# T cell markers
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("CD3D", "TRBC2"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Clusters 0, 3, 8, 19, 21, 34, & 40
tcells <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("CD3D"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# NK cell markers
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("KLRD1", "GNLY", "NKG7"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Clusters 4, 9, & 29
nkcells <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("NKG7"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# B cell markers
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("MS4A1", "CD79A"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Clusters 23, 24, 36, & 44 (24 and 44 = plasmablasts)
bcells <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("CD79A"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# Mesenchymal cell population markers (includes MSCs, pericytes and fibroblasts)
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("PDGFRB", "ACTA2", "COL6A2", "LUM"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Clusters 6, 10, 16, 18, 26, 32, 41 & 43
mesenchymal <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("COL6A2"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# Monocyte cell markers
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("APOBEC3A", "CFP"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Cluster 11
monocyte <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("CFP"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# Macrophage cell markers
FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("CD68", "TREM2"),
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)
# Clusters 1, 15, 33, & 39
macrophage <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = c("CD68"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE
)

# Let's clarify what these cycling cells actually are.
library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
pred.seurat <- SingleR(
  test = seurat_integrated@assays[["RNA"]]@counts, 
  ref = hpca.se, 
  assay.type.test = 1,
  labels = hpca.se$label.main
)

# Summarizing the distribution:
table(pred.seurat$labels)
metadata <- seurat_integrated@meta.data
pred.seurat <- pred.seurat[rownames(pred.seurat) %in% rownames(metadata),]
predicted_celltype <- pred.seurat$labels
metadata <- cbind(metadata, predicted_celltype)
metadata$predicted_celltype[metadata$predicted_celltype %in% names(which(table(metadata$predicted_celltype) <= 1000L))] <- "Rare"
seurat_integrated@meta.data <- metadata
DimPlot(seurat_integrated, reduction = "umap", group.by = "predicted_celltype")

markidentplots <- plot_grid(endo, epi, tcells, nkcells, bcells, mesenchymal, monocyte, macrophage, nrow = 2, ncol = 4)
marktitle <- ggdraw() + 
  draw_label("Marker Identification Process", fontface = 'bold', size = 17, x = 0, hjust = -1) +
  theme(plot.margin = margin(0, 0, 0, 0))
plot_grid(marktitle, markidentplots, ncol = 1, rel_heights = c(0.2, 1))

# Rename the clusters after their cell type identity
new_ids <- c(
  "T cells",
  "Macrophages",
  "Epithelial cells",
  "T cells", 
  "NK cells", 
  "Epithelial cells", 
  "Mesenchymal cells", 
  "Endothelial cells", 
  "T cells", 
  "NK cells", 
  "Mesenchymal cells", 
  "Monocytes", 
  "Epithelial cells", 
  "Epithelial cells", 
  "Endothelial cells", 
  "Macrophages", 
  "Mesenchymal cells", 
  "Endothelial cells", 
  "Mesenchymal cells", 
  "T cells", 
  "Epithelial cells", 
  "T cells",
  "Epithelial cells",
  "B cells",
  "B cells",
  "Epithelial cells",
  "Mesenchymal cells",
  "Epithelial cells",
  "Epithelial cells",
  "NK cells",
  "Epithelial cells",
  "Epithelial cells",
  "Mesenchymal cells",
  "Macrophages",
  "T cells",
  "Endothelial cells",
  "B cells",
  "Endothelial cells",
  "Epithelial cells",
  "Macrophages",
  "T cells",
  "Mesenchymal cells",
  "Endothelial cells",
  "Mesenchymal cells",
  "B cells"
)
names(new_ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new_ids)

# Plot the new UMAP
celltype_annotation <- DimPlot(
  seurat_integrated,
  reduction = "umap"
) + ggtitle("Cell Type Annotation") + theme(
  plot.title = element_text(face = "bold", size = 17, hjust = 0.5)
)

# ----
# Downstream Analysis (fgsea + Gene Ontology)

# Load in fgsea
library(fgsea)

# Generate list of differentially expressed genes between cancer and control conditions
DefaultAssay(seurat_integrated) <- "integrated"
dge_ecs <- FindMarkers(
  seurat_integrated,
  ident.1 = "Cancer",
  ident.2 = "Normal",
  group.by = "orig.ident",
  subset.ident = "Macrophages",
  only.pos = FALSE,
  logfc.threshold = 0
)

# Import gene sets used in fgsea
hallmarks <- fgsea::gmtPathways("/Users/hs/Desktop/breast-cancer-analysis/data/hallmark.genesets.v6.1.symbols.gmt") 
kegg <- fgsea::gmtPathways("/Users/hs/Desktop/breast-cancer-analysis/data/kegg.genesets.v6.1.symbols.gmt")
go <- fgsea::gmtPathways("/Users/hs/Desktop/breast-cancer-analysis/data/GOTerms.BP.v6.1.symbols.gmt")
reactome <- fgsea::gmtPathways("/Users/hs/Desktop/breast-cancer-analysis/data/reactome.genesets.v6.1.symbols.gmt")
gene_sets <- go

# Organize genes and rank them by FC
dge_ecs$gene <- rownames(dge_ecs)
dge_ecs <- dge_ecs %>% arrange(desc(avg_log2FC))
fold_changes <- dge_ecs$avg_log2FC
names(fold_changes) <- dge_ecs$gene

# Run preranked gene set enrichment analysis
gsea_ecs <- fgsea(
  pathways = gene_sets,
  stats = fold_changes,
  minSize = 5,
  maxSize = 5000,
  nproc = 2
)

# Create separate matrix where you have only statistically significant gene sets and they are arranged by normalized enrichment score
gsea_sig <- filter(gsea_ecs, padj <= 0.05) %>%
  arrange(NES)

# Score each cell for gene set activity and see how the gene set is expressed across a population of cells (optional)
seurat_integrated <- AddModuleScore(
  seurat_integrated,
  features = gene_sets["GO_ANGIOGENESIS"],
  name = "Angiogenesis_GS"
)

# Visualize gene set
hist(seurat_integrated$Angiogenesis_GS1, breaks = 50)
FeaturePlot(seurat_integrated, features = "Angiogenesis_GS1", cols = c('lightgrey', 'red'), order = TRUE, split.by = "orig.ident")
VlnPlot(seurat_integrated, features = "Angiogenesis_GS1", pt.size = 0.05)
VlnPlot(
  seurat_integrated,
  features = "Angiogenesis_GS1", 
  idents = "Endothelial cells",
  split.by = "orig.ident"
)

gsea_sig <- gsea_sig[!6:24,]
gsea_sig$pathway <- c(
  "Protein-DNA Complex Subunit Organization",
  "DNA Conformation Change",
  "Chromatin Assembly or Disassembly",
  "DNA Packaging",
  "Chromatin Organization",
  "Tube Morphogenesis",
  "Angiogenesis",
  "Circulatory System Development",
  "Vasculature Development",
  "Blood Vessel Morphogenesis"
)

options(warn = 1)
bar_gsea <- ggplot(gsea_sig, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score") +
  theme_minimal() +
  ggtitle("Top 5 Most Positively and Negatively 
Enriched Gene Sets In Cancer 
Phenotype (p-value < 0.05)") +
  theme(
    axis.text.y.left = element_text(face = c('plain', 'plain', 'plain', 'plain', 'plain', 'bold', 'bold', 'bold', 'bold', 'bold'), size = 15),
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    axis.title = element_text(size = 13),
    axis.text.x.bottom = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(0.5, 1, 0.25, 0.25, "cm")
  )

tubemorph <- plotEnrichment(
  pathway = gene_sets[["GO_TUBE_MORPHOGENESIS"]],
  stats = fold_changes,
  gseaParam = 1
) +
  xlab("Ranked List of Genes") +
  ylab("Enrichment Score") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, face = "bold"), 
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  ggtitle("Tube Morphogenesis") +
  scale_x_continuous(
    breaks = c(0, 500, 1000, 1500),
    labels = c("0\nC", "500", "1000", "1500\nN")
  )

circsys <- plotEnrichment(
  pathway = gene_sets[["GO_CIRCULATORY_SYSTEM_DEVELOPMENT"]],
  stats = fold_changes,
  gseaParam = 1
) +
  xlab("Ranked List of Genes") +
  ylab("Enrichment Score") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, face = "bold"), 
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  ggtitle("Circulatory System Development") +
  scale_x_continuous(
    breaks = c(0, 500, 1000, 1500),
    labels = c("0\nC", "500", "1000", "1500\nN")
  )

bloodmorph <- plotEnrichment(
  pathway = gene_sets[["GO_BLOOD_VESSEL_MORPHOGENESIS"]],
  stats = fold_changes,
  gseaParam = 1
) +
  xlab("Ranked List of Genes") +
  ylab("Enrichment Score") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, face = "bold"), 
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  ggtitle("Blood Vessel Morphogenesis") +
  scale_x_continuous(
    breaks = c(0, 500, 1000, 1500),
    labels = c("0\nC", "500", "1000", "1500\nN")
  )

angio <- plotEnrichment(
  pathway = gene_sets[["GO_ANGIOGENESIS"]],
  stats = fold_changes,
  gseaParam = 1
) +
  xlab("Ranked List of Genes") +
  ylab("Enrichment Score") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, face = "bold"), 
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  ggtitle("Angiogenesis") +
  scale_x_continuous(
    breaks = c(0, 500, 1000, 1500),
    labels = c("0\nC", "500", "1000", "1500\nN")
  )

vascsys <- plotEnrichment(
  pathway = gene_sets[["GO_VASCULATURE_DEVELOPMENT"]],
  stats = fold_changes,
  gseaParam = 1
) +
  xlab("Ranked List of Genes") +
  ylab("Enrichment Score") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, face = "bold"), 
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  ggtitle("Vasculature Development") +
  scale_x_continuous(
    breaks = c(0, 500, 1000, 1500),
    labels = c("0\nC", "500", "1000", "1500\nN")
  )

gsea_plots <- plot_grid(bloodmorph, vascsys, circsys, angio, tubemorph)

gseatitle <- ggdraw() + draw_label("   Gene Set Enrichment Plots", fontface = 'bold', size = 17, x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 3))
plot_grid(gseatitle, gsea_plots, ncol = 1, rel_heights = c(0.1, 1))

gsea_rel <- gsea_sig[6:10, ]
gsea_rel <- gsea_rel %>% arrange(desc(NES))
as.data.frame(gsea_rel)

# Leading Edge Analysis
leadblood <- as.vector(unlist((gsea_rel[1, 8][[1]])))
leadvasc <- as.vector(unlist((gsea_rel[2, 8][[1]])))
leadcirc <- as.vector(unlist((gsea_rel[3, 8][[1]])))
leadang <- as.vector(unlist((gsea_rel[4, 8][[1]])))
leadtube <- as.vector(unlist((gsea_rel[5, 8][[1]])))

# For each vector, get a vector of values without duplicates
deduplicated_vectors <- lapply(list(leadtube, leadcirc, leadblood, leadang, leadvasc), unique)

# Flatten the lists, then sort and use rle to determine how many
# lists each value appears in
rl <- rle(sort(unlist(deduplicated_vectors)))

# Get the values that appear in one or more lists
lead1 <- rl$values[rl$lengths == 1]
lead2 <- rl$values[rl$lengths == 2]
lead3 <- rl$values[rl$lengths == 3]
lead4 <- rl$values[rl$lengths == 4]
lead5 <- rl$values[rl$lengths == 5]
leadnumbers <- list(lead1, lead2, lead3, lead4, lead5)

sharedleadedge <- Reduce(intersect, list(leadtube, leadcirc, leadblood, leadang, leadvasc))
sharedleadedge

# ----
# Downstream Analysis (NicheNet)

# We’ll perform a NicheNet analysis to infer ligand-receptor interactions
# where the sender cell population is Macrophages and the receiver is Endothelial cells.

# Load NicheNet
library(nichenetr)

# 1) Define the "receiver" cell type (Endothelial cells) and the "sender" cell type (Macrophages).
receiver_celltype <- "Endothelial cells"
sender_celltypes <- c("Macrophages")

# 2) Identify DE genes in the receiver cell population of interest
#    under the condition of interest ("Cancer" vs "Normal").
Idents(seurat_integrated) <- "orig.ident"
receiver_dge <- FindMarkers(
  seurat_integrated,
  ident.1 = "Cancer",
  ident.2 = "Normal",
  subset.ident = receiver_celltype,
  min.pct = 0.1
)
receiver_dge$gene <- rownames(receiver_dge)

# Define the set of genes that are significantly upregulated in our receiver population.
geneset_oi <- receiver_dge %>%
  dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
  dplyr::pull(gene) %>%
  unique()

# 3) Define a background set of genes expressed by the receiver population.
#    We'll pick all genes detected above a minimal threshold in Endothelial cells.
Idents(seurat_integrated) <- receiver_celltype
avg_expr_receiver <- AverageExpression(seurat_integrated)$RNA
background_expressed_genes <- names(which(avg_expr_receiver > 0.1))

# 4) Run the NicheNet “aggregate” pipeline, focusing on Endothelial cells as receiver,
#    Macrophages as sender, “Cancer” as condition of interest, and “Normal” as reference.
#    This aggregates expression at the cell-type level automatically.
Idents(seurat_integrated) <- "seurat_clusters"  # or whichever metadata column best separates your clusters
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = seurat_integrated,
  receiver = receiver_celltype,
  sender = sender_celltypes,
  assay_oi = "RNA",
  organism = "human",
  condition_colname = "orig.ident",
  condition_oi = "Cancer",
  condition_reference = "Normal",
  geneset = geneset_oi,
  background = background_expressed_genes
)

# 5) Inspect the top ligands
head(nichenet_output$ligand_activities)

# 6) Visualize the top ligand-target interactions
#    (Here we pick the top 10 ligands and plot their target regulatory potential scores in a heatmap.)
top_ligands <- nichenet_output$top_ligands[1:10]
ligand_target_heatmap <- nichenet_output$ligand_target_matrix %>%
  .[top_ligands, ] %>%
  as.matrix()

pheatmap::pheatmap(
  ligand_target_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top 10 Ligand-Target Interactions\n(Macrophages -> Endothelial Cells)"
)

# 7) Optionally, create a dot plot showing how these ligands are expressed in Macrophages,
#    confirming whether they can produce the relevant ligands.
DotPlot(
  object = seurat_integrated,
  features = unique(top_ligands),
  group.by = "ident"  # or your chosen metadata column
) +
  coord_flip() +
  ggtitle("Ligand Expression in Macrophages")
