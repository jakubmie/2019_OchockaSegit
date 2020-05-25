####################################
# sc.mm.gam.gender                 #
# data preprocessing and analysis  #
####################################

# Load packages
library(parallel)
library(future)
library(Seurat)
library(Matrix)
options(future.globals.maxSize=256*1024^3)
plan(multiprocess)
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)
library(purrr)

# Analysis parameters
mc <- 32                     # number of given cores
anchor_dims <- 30            # number of anchor dimensions used for biological replicates integration
pca_dims <- 30               # number of PCA dimensions to compute and use in tSNE, UMAP                             
umap_n_neighbors <- 30       # UMAP parameter,                                                                    
clustering_resolution <- 0.3 # Resolution parameter for Seurat clustering
n_features <- 2000

######## Load external data ########
# Read gene annotations
annot <- read.csv("mm_annot_20190131.txt", sep = "\t", stringsAsFactors = F)

# Read a list of cell cycle markers, from Tirosh et al, 2015
cc_genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
cc_genes <- paste0(substring(cc_genes, 1, 1), tolower(substring(cc_genes, 2, nchar(cc_genes))))
cc_genes_annot <- dplyr::left_join(x = as.data.frame(cc_genes, stringsAsFactors=F), y = annot, 
                                   by = c("cc_genes"="Gene.name"))

# Divide cell cycle genes list into markers of G2/M phase and markers of S phase
s_genes <- cc_genes_annot[1:43, "Gene.stable.ID"]
s_genes <- s_genes[!is.na(s_genes)]
g2m_genes <- cc_genes_annot[44:97, "Gene.stable.ID"]
g2m_genes <- g2m_genes[!is.na(g2m_genes)]

######## Load and preprocess expression data ########
# Define selected (previously reported) microglia and macrophages markers
microglia_markers <- annot[match(c("Tmem119", "P2ry12", "Sall1", "Pros1", "Crybb1"), annot$Gene.name), 
                           "Gene.stable.ID"]
macrophages_markers <- annot[match(c("Itga4", "Tgfbi", "Cxcl2", "Ccr2", "Il10", "Fgr"), annot$Gene.name), 
                             "Gene.stable.ID"]

# Read raw data (gene/cell count matrix from cellranger, filtered: use only detected cellular barcodes)
samples <- dir(path=".", pattern="filtered_feature_bc_matrix$", full.names=T, recursive=T, include.dirs=T)
samples_raw_data <- mclapply(samples, function(s) {
    matrix_dir = s
    
    # Set paths for files with: barcodes, features (gene names), expression (count) matrix
    barcode_path <- paste0(matrix_dir, "/", "barcodes.tsv.gz")
    features_path <- paste0(matrix_dir, "/", "features.tsv.gz")
    matrix_path <- paste0(matrix_dir, "/", "matrix.mtx.gz")
    
    # Read expression matrix, feature (gene) names and barcodes  
    mat <- readMM(file = matrix_path)
    feature_names = read.delim(features_path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode_names = read.delim(barcode_path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    
    # Set cell barcodes as column names and gene ids (Ensembl stable ID) as row names       
    colnames(mat) = barcode_names$V1
    rownames(mat) = feature_names$V1 # V1 - ENSMUSG, V2 - gene_names
    
    # Output gene/cell count matrix - genes in rows, cells in columns
    mat
}, mc.cores = mc)
names(samples_raw_data) <- substr(samples, 3, nchar(samples) - 27)

# Set up Seurat objects
samples_objects <- mclapply(seq_along(samples_raw_data), function(i) {
    seurat_object <- CreateSeuratObject(counts = samples_raw_data[[i]], 
                                        project = paste0("mm.gam.gender.", names(samples_raw_data)[i]), 
                                        min.cells = 5)
    
    # Add 'condition' value based on sample name
    if(grepl("ctrl", names(samples_raw_data)[i])) {
        seurat_object$condition <- "ctrl"
    } else {
        seurat_object$condition <- "tumor"
    }
    seurat_object
}, mc.cores = mc)
names(samples_objects) <- names(samples_raw_data)

# Analyze percentage of mitochondrial genes in cells
mito_features <- annot[grep(pattern = "^mt-", annot$Gene.name), ][, 1]
percent_mito <- mclapply(seq_along(samples_objects), function(i) {
    Matrix::colSums(x = GetAssayData(object = samples_objects[[i]], slot = 'counts')
                    [rownames(samples_objects[[i]]) %in% mito_features, ]) / 
        Matrix::colSums(x = GetAssayData(object = samples_objects[[i]], slot = 'counts'))
})

samples_objects <- mclapply(seq_along(samples_objects), function(i) {
    samples_objects[[i]][['percent_mito']] <- percent_mito[[i]]	
    samples_objects[[i]]
}, mc.cores = mc)
names(samples_objects) <- names(samples_raw_data)

# Filter cells (based on percent_mito and nFeature_RNA), normalize (log-normalize in columns) 
# and find n_features most variable features/genes (using Variance-stabilizing transformation)
# to the list of top variable genes add markers of microglia, macrophages and cell cycle
samples_objects <- mclapply(seq_along(samples_objects), function(i) {
    samples_objects[[i]] <- subset(x = samples_objects[[i]],
                                   subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent_mito < 0.05)
    samples_objects[[i]] <- NormalizeData(object = samples_objects[[i]], verbose = FALSE)
    samples_objects[[i]] <- FindVariableFeatures(object = samples_objects[[i]],
                                                 selection.method = "vst", 
                                                 nfeatures = n_features, verbose = FALSE)
    # Add sample shortID
    samples_objects[[i]]$shortID <- substr(as.character(samples_objects[[i]]$orig.ident), 15,
                                                        nchar(as.character(samples_objects[[i]]$orig.ident)))
    
    # Add meta data: sex and sex_condition
    samples_objects[[i]]$sex <- substr(as.character(samples_objects[[i]]$shortID), 1, 1)
    samples_objects[[i]]$sex_condition <- paste0(samples_objects[[i]]$sex, "_", samples_objects[[i]]$condition)
    
    # Add selected (previously reported) genes to var.features
    samples_objects[[i]]@assays$RNA@var.features <- unique(c(samples_objects[[i]]@assays$RNA@var.features, 
                                                           s_genes, g2m_genes, microglia_markers, macrophages_markers))
    samples_objects[[i]]
}, mc.cores = mc)
names(samples_objects) <- names(samples_raw_data)

# Integrate replicates within conditions, scale, regress out unwanted sources of variation, 
# calculate PCA, t-SNE and find cell clusters
sex_condition_objects <- mclapply(c(1,3,5,7), function(i) { # f_ctrl_1, f_tumor_1, m_ctrl_1, m_tumor_1
    sex_condition <- unique(samples_objects[[i]]$sex_condition)
    print(paste0("Find Integration Anchors: ", sex_condition))
    samples_anchors <- FindIntegrationAnchors(object.list = list(samples_objects[[i]], 
                                                                 samples_objects[[i + 1]]), 
                                                                 dims = 1:anchor_dims,
                                                                 anchor.features = n_features)
    print(paste0("Integrate Data: ", sex_condition))
    samples_integrated <- IntegrateData(anchorset = samples_anchors, dims = 1:anchor_dims)
    DefaultAssay(object=samples_integrated) <- "integrated"
    
    print(paste0("Cell Cycle Scoring:", sex_condition))
    samples_integrated <- CellCycleScoring(object = samples_integrated, s.features = s_genes, 
                                           g2m.features = g2m_genes, set.ident = TRUE)  
    
    print(paste0("Scale Data: ", sex_condition, 
                 " regress out: nCount_RNA, percent_mito, CC_difference"))

    # Calculate difference between G2M and S phase scores to separate non-cycling and cycling cells
    # Approch described in Seurat's Cell-Cycle Scoring and Regression vignette 
    samples_integrated$CC_Difference <- samples_integrated$S.Score - samples_integrated$G2M.Score
    samples_integrated <- ScaleData(object = samples_integrated, verbose = FALSE, 
                                    vars.to.regress = c("nCount_RNA", 
                                                        "percent_mito", 
                                                        "CC_Difference"))  
    
    print(paste0("Run PCA: ", sex_condition))
    samples_integrated <- RunPCA(object = samples_integrated, npcs = pca_dims, verbose = FALSE)
    
    print(paste0("Run t-SNE:", sex_condition))
    samples_integrated <- RunTSNE(object = samples_integrated, reduction = "pca", dims=1:pca_dims)
    
    print(paste0("Find Neighbors: ", sex_condition))
    samples_integrated <- FindNeighbors(object = samples_integrated, dims = 1:pca_dims)
    
    print(paste0("Find Clusters: ", sex_condition))
    samples_integrated <- FindClusters(object = samples_integrated, resolution = clustering_resolution) 
    samples_integrated
}, mc.cores = mc)
names(sex_condition_objects) <- substring(names(samples_objects), 1, nchar(names(samples_objects))-2)[c(1,3,5,7)]

# Find markers for obtained clusters
markers_sex_condition_objects <- mclapply(sex_condition_objects, function(x) {
    FindAllMarkers(object = x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
}, mc.cores = mc)

# Select cells (Microglia, Macrophages and BAMs) for downstream analysis
cell_types<-read.csv("cell_type_index.csv", row.names = 1, sep=";")

sex_condition_objects <- lapply(sex_condition_objects, function(x) {
    x$full_cluster_id <- paste(substring(x$shortID,1,1), x$condition, Idents(x), sep="_")
    x$cell_type <- cell_types[match(x$full_cluster_id, cell_types$cluster), "cell_type"]
    x$cell_type <- factor(x$cell_type, levels= c("micro", "pre-micro", "macro", "BAM", "NKT", "NK",
                                                 "B-cells", "T-cells","Ncam1+", "DC", "other"))
    x$cell_type_selection <- ""
    x$cell_type_selection[x$cell_type %in% c("micro", "pre-micro")] <- "Microglia"
    x$cell_type_selection[x$cell_type == "macro"] <- "Macrophages"
    x$cell_type_selection[x$cell_type == "BAM"] <- "BAM"
    x
})

# Create one Seurat object with selected cells from all samples
# Define clusters to be merged
# Microglia + macrophages + BAM
clusters_to_take <- list(
    f_ctrl = c(0,1,2,3,4,5,6,8),
    f_tumor = c(0,1,2,3,5,6,7),
    m_ctrl = c(0,1,2,3,6,7),
    m_tumor = c(0,1,2,3,4,5,6,7)
)

# Fetch cells from clusters defined in clusters_to_take list
object_list <- sex_condition_objects
clusters_taken_list <- unlist(mclapply(object_list[names(object_list) %in% names(clusters_to_take)], function(condition_object) {
    unique_shortIDs <- unique(condition_object$shortID) # Fetch sample ids present in the condition
    samples_expr_in_clusters <- unlist(mclapply(unique_shortIDs, function(sample_name) {
        # Condition name
        condition_name <- substring(sample_name, 1, nchar(sample_name) - 2)  
        
        # Which clusters to take for a given condition
        clusters <- clusters_to_take[[condition_name]] # same clusters for both replicates
        
        # Subset object to take only selected clusters
        sample_cluster_expr <- SubsetData(object = condition_object, ident.use = clusters, 
                                          subset.name = "shortID", accept.value = sample_name)
        
        # Set default assay as RNA and remove previously calculated elements
        DefaultAssay(object = sample_cluster_expr) <- "RNA"
        sample_cluster_expr@assays$integrated  <- NULL
        sample_cluster_expr@tools  <- list()
        sample_cluster_expr@commands  <- list()
        sample_cluster_expr@reductions  <- list()
        sample_cluster_expr@neighbors  <- list()
        sample_cluster_expr@graphs  <- list()
        
        # Seurat object subsetted to defined clusters, only RNA assay and metadata
        sample_cluster_expr
    }, mc.cores = mc))
    names(samples_expr_in_clusters) <- unique_shortIDs
    samples_expr_in_clusters
}, mc.cores = mc))
names(clusters_taken_list) <- substring(names(clusters_taken_list), regexpr("[.]", names(clusters_taken_list))+1)

# Normalize (log-normalize in columns) 
# Find n_features most variable features/genes (using Variance-stabilizing transformation)
clusters_objects <- mclapply(clusters_taken_list, function(cluster_sample_object) {
    cluster_sample_object <- NormalizeData(object = cluster_sample_object, verbose = FALSE)
    cluster_sample_object <- FindVariableFeatures(object = cluster_sample_object, 
                                                  selection.method = "vst", 
                                                  nfeatures = n_features, verbose = FALSE)

    cluster_sample_object@assays$RNA@var.features <- unique(c(cluster_sample_object@assays$RNA@var.features, 
                                                              s_genes, g2m_genes))

    cluster_sample_object
}, mc.cores = mc)

# Merge cells from selected clusters from all samples 
seu_object <- merge(x = clusters_objects[[1]], y = clusters_objects[2:8], add.cell.ids = c(paste0("_", 1:8)))
seu_object <- FindVariableFeatures(object = seu_object, 
                                                 selection.method = "vst", 
                                                 nfeatures = n_features, verbose = FALSE)
seu_object@assays$RNA@var.features <- unique(c(seu_object@assays$RNA@var.features, 
                                                             s_genes, g2m_genes))

print(paste0("Cell Cycle Scoring:"))
seu_object <- CellCycleScoring(object = seu_object, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)  

print(paste0("Scale Data:", " regress out: nCount_RNA, percent_mito, CC_Difference"))
seu_object$CC_Difference <- seu_object$S.Score - seu_object$G2M.Score
seu_object <- ScaleData(object = seu_object, verbose = FALSE, 
                                      vars.to.regress = c("nCount_RNA", 
                                                          "percent_mito",
                                                          "CC_Difference"))  

print(paste0("Run PCA:"))
seu_object <- RunPCA(object = seu_object, npcs = pca_dims, verbose = FALSE)

print(paste0("Run UMAP:"))
seu_object <- RunUMAP(object = seu_object, reduction = "pca", 
                                    dims=1:pca_dims, n.neighbors = umap_n_neighbors, min.dist = 0.5)

print(paste0("Find Neighbors:"))
seu_object <- FindNeighbors(object = seu_object, dims = 1:pca_dims)
print(paste0("Find Clusters:")) # double resolution for merged dataset
seu_object <- FindClusters(object = seu_object, resolution = 2*clustering_resolution) 

markers_seu_object <- FindAllMarkers(object = seu_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

seu_object$cell_type_4_groups <- plyr::mapvalues(Idents(seu_object), 
                                           from=c(0:14), 
                                           to = c("h.Microglia",  "h.Microglia",  "h.Microglia",  "Macrophages",  "a.Microglia", "h.Microglia",  "a.Microglia", "BAM", "h.Microglia", "a.Microglia", "Macrophages", "a.Microglia", "Macrophages", "Macrophages",  "h.Microglia"))
seu_object$cell_type_4_groups <- factor(seu_object$cell_type_4_groups, levels = levels(seu_object$cell_type_4_groups)[c(4,2,3,1)])

seu_object$cell_type_3_groups <- plyr::mapvalues(Idents(seu_object), 
                                                  from=c(0:14), 
                                                  to = c("Microglia",  "Microglia",  "Microglia",  "Macrophages",  "Microglia", "Microglia",  "Microglia", "BAM", "Microglia", "Microglia", "Macrophages", "Microglia", "Macrophages", "Macrophages",  "Microglia"))
seu_object$cell_type_3_groups <- factor(seu_object$cell_type_3_groups, levels = levels(seu_object$cell_type_3_groups)[c(3,2,1)])

Idents(seu_object) <- seu_object$cell_type_4_groups

markers_cell_types_3_groups <- FindAllMarkers(object = seu_object, only.pos = TRUE, min.pct = 0.25,
                                              logfc.threshold = 0.25)

get_top10_markers_gnames <- function(markers)  {
    markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% arrange(cluster, desc(avg_logFC))
    markers <- dplyr::left_join(x = markers, y = annot, by = c("gene"="Gene.stable.ID"))
    markers
}

markers_cell_types_3_groups_top10 <- get_top10_markers_gnames(markers_cell_types_3_groups)

# Gene Ontology analysis
library(clusterProfiler)
library(org.Mm.eg.db)

# Load annotation file with NCBI ids
annotEnsemblNCBI <- read.csv("mm_annot_entrezID_20190130.txt", sep = "\t", stringsAsFactors = F)

addGeneNames <- function(markers)  {
    markers <- dplyr::left_join(x = markers, y = annot, by = c("gene"="Gene.stable.ID"))
    markers
}

markers_ActMG_vs_HomMG <- FindMarkers(object = seu_object, ident.1 = "a.Microglia", 
                                      ident.2 = "h.Microglia", only.pos = TRUE, min.pct = 0.25, 
                                      logfc.threshold = 0.25)
markers_ActMG_vs_HomMG$gene <- rownames(markers_ActMG_vs_HomMG)
markers_ActMG_vs_HomMG <- addGeneNames(markers_ActMG_vs_HomMG)

markers_MoM_vs_ActMG <- FindMarkers(object = seu_object, ident.1 = "Macrophages", 
                                    ident.2 = "a.Microglia", only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
markers_MoM_vs_ActMG$gene <- rownames(markers_MoM_vs_ActMG)
markers_MoM_vs_ActMG <- addGeneNames(markers_MoM_vs_ActMG)

pvalue_threshold <- 0.05
object_tmp <- list(markers_ActMG_vs_HomMG, markers_MoM_vs_ActMG)
names(object_tmp) <- c("ActMG_vs_HomMG", "MoM_vs_ActMG")
GO_markers_ActMG_MoM <- lapply(names(object_tmp), function(group) {
    markers <- object_tmp[[group]]
    markers_to_GO <- markers
    markers_to_GO <- dplyr::left_join(x = markers_to_GO, y = annotEnsemblNCBI[, c("Gene.stable.ID", "NCBI.gene.ID")], by = c("gene"="Gene.stable.ID"))
    
    markers_GO_NCBI <- markers_to_GO[markers_to_GO$p_val_adj < pvalue_threshold, ]$NCBI.gene.ID
    
    universe_GO <- data.frame(gene = unique(seu_object@assays$RNA@data %>% rownames), stringsAsFactors = F) # all genes as universe
    universe_GO <- dplyr::left_join(x = universe_GO, y = annotEnsemblNCBI[, c("Gene.stable.ID", "NCBI.gene.ID")], by = c("gene"="Gene.stable.ID"))
    
    g <- as.character(unique(unlist(markers_GO_NCBI)))
    g <- g[!is.na(g)]
    u <- as.character(unique(unlist(universe_GO$NCBI.gene.ID)))
    u <- u[!is.na(u)]
    
    group %>% print
    
    if((length(g) > 0) && (length(u) > 0)) {
        GO_upregulated <-  enrichGO(gene = g[!is.na(g)], pvalueCutoff = pvalue_threshold, OrgDb = org.Mm.eg.db, universe = u[!is.na(u)], ont = "BP") %>% setReadable(OrgDb = org.Mm.eg.db)
    } else {
        GO_upregulated <- NULL
    }
    GO_upregulated
})
names(GO_markers_ActMG_MoM) <- names(object_tmp)

GO_markers_ActMG_MoM_simplify <- lapply(GO_markers_ActMG_MoM, function(GO_results) {
    GO_results@result <- GO_results@result[GO_results@result$p.adjust <= 0.05, ]
    simplify(GO_results, cutoff=0.7, by="p.adjust", select_fun=min)
})

gene_lists <- lapply(object_tmp, function(x) {
    result <- exp(1)^x$avg_logFC
    names(result) <- x$gene
    result <- sort(result, decreasing = TRUE)
})
names(gene_lists) <- names(object_tmp)









