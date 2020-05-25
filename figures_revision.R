####################################
# sc.mm.gam.gender                 #
# Figures - revision               #
####################################

#### colors
micro<-"#53AFE6"
pre_micro<-"#2DA7C8"
BAM<- "#0DD1AD"
UN<-"grey"
Mo<-"#FCE80C"
Mo_Mg<-"#FABF00"
Mg<-"#E98934"
NK<-"#8c42a3"
ncam<-"#C2B4FC"
NKT<-"#DFA5F2"
DC<-"#bf7a58"
Tcells<-"#94112f" 
Bcells<-"#EC5CA5"

col_Micro<-"#53AFE6"
col_Macro<-"#FABF00"
col_BAM<-"#0DD1AD"

col_ActMG<-"#2F8EA1"
col_HomMG<- "#AFCADB"

col_CTRL<- "#636775"
col_TUMOR<-"#F44686"

col_male<- "#4FCAFF"
col_female<- "#FFFF00"

col_m_CTRL<-"#3B47B5"
col_m_TUMOR<- "#057D9F"
col_f_CTRL<-"#FFA739"
col_f_TUMOR<-"#FFDC6B"

# Figure 1b
# f_ctrl
DimPlot(sex_condition_objects$f_ctrl)
# f_tumor
DimPlot(sex_condition_objects$f_tumor)
# m_ctrl
DimPlot(sex_condition_objects$m_ctrl)
# m_tumor
DimPlot(sex_condition_objects$m_tumor)

# Figure 1c
# ctrl objects
ctrl_gene_panel <- c("Cd14", "Tmem119", "P2ry12", "Csf1", "Crybb1", "Mcm5", "Ifit3", "Tgfbi", 
                     "Ifitm2", "Ifitm3", "S100a6", "Ly6c2", "Ccr2", "Mrc1", "Cd163", "Cd24a", 
                     "Ncam1")
DotPlot(sex_condition_objects$f_ctrl, features = 
            rev(annot[match(ctrl_gene_panel, annot$Gene.name), "Gene.stable.ID"])) + 
        RotatedAxis() +
        theme(legend.position="top")

DotPlot(sex_condition_objects$m_ctrl, features = 
            rev(annot[match(ctrl_gene_panel, annot$Gene.name), "Gene.stable.ID"])) + 
        RotatedAxis() +
        theme(legend.position="top")

# tumor objects
tumor_gene_panel <- c("Cd14", "Tmem119", "P2ry12", "Tgfbi", "Ifitm2", "Ifitm3", "S100a6", "Ly6c2",
                      "Ccr2", "Mrc1", "Cd163", "Cd24a", "Ncam1", "Klrk1", "Ncr1", "Cd2", "Cd3d",
                      "Cd4", "Cd8b1", "Ms4a1")
DotPlot(sex_condition_objects$f_tumor, features = 
            rev(annot[match(tumor_gene_panel, annot$Gene.name), "Gene.stable.ID"])) + 
            RotatedAxis() +
            theme(legend.position="top")

DotPlot(sex_condition_objects$m_tumor, features = 
            rev(annot[match(tumor_gene_panel, annot$Gene.name), "Gene.stable.ID"])) + 
            RotatedAxis() +
            theme(legend.position="top")

# Figure 1d
# Pie charts
freq_list <- lapply(sex_condition_objects, function(x) {
    freq <- data.frame(cell_type = x$cell_type)
    freq <- freq %>%
        group_by(cell_type) %>%
        count() %>%
        ungroup %>%
        mutate(per = `n`/sum(`n`))
    freq$cell_type <- factor(freq$cell_type, levels= c("micro", "pre-micro", "macro", "BAM", "NKT",
                                                       "NK", "B-cells", "T-cells","Ncam1+", "DC", 
                                                       "other"))
    freq$label <- scales::percent(freq$per)
    freq
})

cf<-ggplot(freq_list$f_ctrl, 
           aes(x="", y=per, fill=cell_type))+
    geom_bar(stat="identity", width=1, color="white")+
    coord_polar("y", start=0)+
    scale_fill_manual(values=c(micro, pre_micro, BAM, UN))+
    theme_light()+
    geom_label_repel(aes(label = label), size=3, show.legend = F, nudge_x = 1)

cm<-ggplot(freq_list$m_ctrl, 
           aes(x=" ", y=per, fill=cell_type))+
    geom_bar(stat="identity", width=1, color="white")+
    coord_polar("y", start=0)+
    scale_fill_manual(values=c(micro, pre_micro, Mo_Mg, BAM, NK, DC))+
    theme_light()+
    geom_label_repel(aes(label = label), size=3, show.legend = F, nudge_x = 1)

tf<-ggplot(freq_list$f_tumor, 
           aes(x="", y=per, fill=cell_type))+
    geom_bar(stat="identity", width=1, color="white")+
    coord_polar("y", start=0)+
    scale_fill_manual(values=c(micro, Mo_Mg, BAM,  NKT, NK, Bcells, Tcells, ncam, DC))+
    theme_light()+
    geom_label_repel(aes(label = label), size=3, show.legend = F, nudge_x = 1)

tm<-ggplot(freq_list$m_tumor, 
           aes(x="", y=per, fill=cell_type))+
    geom_bar(stat="identity", width=1, color="white")+
    coord_polar("y", start=0)+
    scale_fill_manual(values=c(micro, Mo_Mg, BAM,  NKT, ncam, DC, Tcells, Bcells))+
    theme_light()+
    geom_label_repel(aes(label = label), size=3, show.legend = F, nudge_x = 1)

ggarrange(cf, cm, tf, tm, ncol = 4)

# Figure 2a
# left panel plots
# cells selected from each sex-condition for further analysis, identified as MG, MoM, BAM
# f_ctrl
DimPlot(sex_condition_objects$f_ctrl, group.by = "cell_type_selection")
# f_tumor
DimPlot(sex_condition_objects$f_tumor, group.by = "cell_type_selection")
# m_ctrl
DimPlot(sex_condition_objects$m_ctrl, group.by = "cell_type_selection")
# m_tumor
DimPlot(sex_condition_objects$m_tumor, group.by = "cell_type_selection")

# right panel plot
DimPlot(seu_object, group.by = "cell_type_selection")

# Figure 2b
# Heatmap with expression of the top 10 markers (ordered by avg_logFC) for each cell type
library(ComplexHeatmap)
library(circlize)

markers_expression_data <- GetAssayData(object = seu_object, slot = "data")
markers_expression_data <- markers_expression_data[match(markers_cell_types_3_groups_top10$gene,
                                                         rownames(markers_expression_data)),
                                                   order(Idents(seu_object))]
markers_expression_data <- t(apply(markers_expression_data, 1, function(x) {
                                                               x/(quantile(x, probs=0.999))}))
rownames(markers_expression_data) <- annot[match(rownames(markers_expression_data), 
                                                 annot$Gene.stable.ID), "Gene.name"]

# Cap values
markers_expression_data[markers_expression_data > 1] <- 1  
markers_expression_data[markers_expression_data < 0] <- 0

col_annotations <- data.frame(cell_id = colnames(seu_object),
                              cluster_id = Idents(seu_object))
rownames(col_annotations) <- col_annotations$cell_id

col_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"))


ht_list <- Heatmap(markers_expression_data, name = "a", 
                   top_annotation = 
                      HeatmapAnnotation(
                         Cluster = col_annotations[colnames(markers_expression_data), "cluster_id"], 
                         col = list(Cluster = c("Microglia" = col_Micro, 
                                                "Macrophages" = col_Macro, 
                                                "BAM" = col_BAM)),
                         simple_anno_size = unit(3, "cm"),
                         annotation_legend_param = list(Cluster = list(direction = "horizontal"), 
                                                        title_gp = gpar(fontsize = 15), 
                                                        labels_gp = gpar(fontsize = 15), 
                                                        legend_height = unit(5, "cm"), ncol = 3)),
                   cluster_rows = F, cluster_columns = F, show_column_names = F, 
                   row_names_gp = gpar(fontsize = 25), col = col_fun, row_split = 
                                                                 rep(c("A", "B", "C"), each = 10),
                   heatmap_legend_param = list(title = "Expression", direction = "horizontal", 
                                               legend_height = unit(3, "cm"), 
                                               title_gp = gpar(fontsize = 15), 
                                               labels_gp = gpar(fontsize = 15), 
                                               at = c(0, 1), labels = c("Low", "High")), 
                                               use_raster = T, raster_device = "png"
)

draw(ht_list, merge_legend = T, heatmap_legend_side = "top", annotation_legend_side = "top")

# Figure 2c, 2d, 2f
# Feature Plots
# 2c: Tmem119, P2ry12
# 2d: Lgals3, Ptprc
# 2f: Ccr2, Ly6c2, Ly6i, Tgfbi, Ccl5, Cd274, Ccl22, Itga4
gene_to_plot <- "Tmem119"
FeaturePlot(seu_object, features = annot[annot$Gene.name == gene_to_plot, "Gene.stable.ID"]) + 
    scale_color_gradientn(colors=viridis::viridis(n=50)) +
    ggtitle(gene_to_plot)

# Figure 2g
# Density Plots
genes<-c("Ccr2", "Ly6c2","Ly6i","Tgfbi","Ccl5", "Cd274", "Ccl22", "Itga4")
gene_expression_data <- GetAssayData(object = seu_object, slot = "data")
gene_expression_data <- as.data.frame(t(gene_expression_data[annot[match(genes, annot$Gene.name),
                                                                   "Gene.stable.ID"], ]))
colnames(gene_expression_data) <- genes
gene_expression_data$cell_types_3_groups <- Idents(seu_object)
gene_expression_data$clusters_0.6 <- seu_object$RNA_snn_res.0.6
gene_expression_data <- gene_expression_data[gene_expression_data$cell_types_3_groups == "Macrophages", ]

gene_expression_data$MoM_subtype <- NULL
gene_expression_data$MoM_subtype[gene_expression_data$clusters_0.6 %in% c(10,12)] <- "Mo"
gene_expression_data$MoM_subtype[gene_expression_data$clusters_0.6 == 2] <- "intMoM"
gene_expression_data$MoM_subtype[gene_expression_data$clusters_0.6 == 13] <- "M"
gene_expression_data$MoM_subtype <- factor(gene_expression_data$MoM_subtype)
gene_expression_data$MoM_subtype <- factor(gene_expression_data$MoM_subtype, 
                                      levels = levels(gene_expression_data$MoM_subtype)[c(2,3,1)])

graphs<-list()
for (i in seq_along(genes)){
    names<-c("Ccr2", "Ly6c2", "Ly6i", "Tgfbi",
             "Ccl5", "Cd274(PD-L1)", "Ccl22", "Itga4(CD49d)")
    plot<- ggplot(gene_expression_data[gene_expression_data$cell_types_3_groups == "Macrophages", ],
                  aes(x="", fill = MoM_subtype)) +
           geom_density(aes_string(x=as.name(genes[i])), trim=T, alpha=0.8, size=0.5)+
           scale_fill_manual(values= c( "#FFDB3F", "orange", "#F0627C"))+
           labs(title=names[i], x="", y="")+
           ylim(0,1)+
           theme_classic()+
           theme(legend.title = element_blank())
    graphs[[i]] <- plot
}
ggarrange(plotlist=graphs, ncol=2, nrow=4, common.legend = T, legend = "top")

# Figure 2j
# DimPlot (section) Mo, intMoM, M; 
DimPlot(seu_object, group.by="RNA_snn_res.0.6")

# Figure 3a
# DimPlot
DimPlot(seu_object, group.by = "condition")

# Figure 3b
gene_expression_data <- GetAssayData(object = seu_object, slot = "data")
genes_micro01<-c("Tmem119", "P2ry12", "Cx3cr1", "Olfml3", "Sparc","Gpr34")
gene_expression_data_micro <- gene_expression_data[annot[match(genes_micro01, annot$Gene.name),
                                                         "Gene.stable.ID"], ]
seu_object$micro_score <- colMeans(gene_expression_data_micro)

genes_macro01<-c("Ifitm2", "S100a6", "S100a11", "Lgals3", "Isg15", "Ms4a4c", "Crip1")
gene_expression_data_macro <- gene_expression_data[annot[match(genes_macro01, annot$Gene.name),
                                                         "Gene.stable.ID"], ]
seu_object$macro_score <- colMeans(gene_expression_data_macro)

# Feature Plots
FeaturePlot(seu_object, features = "micro_score")
FeaturePlot(seu_object, features = "macro_score")

# Figure 3c
# Violin Plots
cell_types_labels <- c("MG", "Mo/M")
names(cell_types_labels) <- c("Microglia", "Macrophages")

scores_data <- data.frame(micro01 = seu_object$micro_score, macro01 = seu_object$macro_score, 
                          condition = seu_object$condition,
                          cell_types_3_groups = Idents(seu_object))

MG<-ggplot(scores_data[scores_data$cell_types_3_groups %in% c("Microglia", "Macrophages"),], 
           aes(x=condition, y= micro01))+
    geom_violin(aes(fill=condition), scale = "area", trim=F, size=0.5)+
    facet_grid(.~cell_types_3_groups, labeller = labeller(cell_types_3_groups=cell_types_labels))+
    xlab("")+
    ylab("MG score")+
    scale_fill_manual(values=c(col_CTRL, col_TUMOR))+
    theme_classic(base_size=14)
MoM<-ggplot(scores_data[scores_data$cell_types_3_groups %in% c("Microglia", "Macrophages"),], 
            aes(x=condition, y= macro01))+
    geom_violin(aes(fill=condition), scale="area", trim=F, size=0.5)+
    facet_grid(.~cell_types_3_groups, labeller = labeller(cell_types_3_groups = cell_types_labels))+
    xlab("")+
    ylab("MoM score")+
    scale_fill_manual(values=c(col_CTRL, col_TUMOR))+
    theme_classic(base_size=14)

ggarrange(MG, MoM, nrow = 2, common.legend=T)

# Figure 3d
# Heatmap with selected, reported markers of MG and MoM
# Hierarchical clustering of MG and MoM cells
selected_genes <- c("Cx3cr1", "P2ry12", "P2ry13", "Tmem119", "Selplg", 
                    "Tgfbi", "Ifitm2", "S100a6", "S100a11")

exclude.clusters <- c(7, 14) # BAM
library(plyr)

object_tmp <- seu_object
Idents(object_tmp) <- object_tmp$RNA_snn_res.0.6
object_tmp <- SubsetData(object = object_tmp, ident.use = 
                             as.integer(setdiff(unique(Idents(object_tmp)), exclude.clusters)))
selected_genes_expression <- GetAssayData(object = object_tmp, slot = "data")
selected_genes_expression <- t(selected_genes_expression[match(
    annot[match(selected_genes, annot$Gene.name), "Gene.stable.ID"], rownames(selected_genes_expression)), ])
selected_genes_expression <- t(apply(selected_genes_expression, 2, function(x) x/(quantile(x, probs=0.999))))
rownames(selected_genes_expression) <- selected_genes

cluster_id <- object_tmp$RNA_snn_res.0.6

col_annotations <- data.frame(cell_id = colnames(object_tmp), 
                              cluster_id = cluster_id,
                              condition = object_tmp$condition)
rownames(col_annotations) <- col_annotations$cell_id
col_annotations$group_id <- plyr::mapvalues (col_annotations$cluster_id, from = c(0,1,2,3,4,5,6,8,9,10,11,12,13,15), 
                                             to = c("micro","micro","macro","micro","micro","micro","micro",
                                                    "micro","micro","macro","micro","macro","macro","micro"))

c_annotations_condition <- col_annotations$condition
names(c_annotations_condition) <- col_annotations$cell_id

c_annotations_group <- col_annotations$group_id
c_annotations_group <- as.character(c_annotations_group)
names(c_annotations_group) <- col_annotations$cell_id 

col_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#edf8fb", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"))

dend <- dendsort(hclust(dist(t(selected_genes_expression)), method = "ward.D2"), isReverse=T)

# Fisher test on two dendrogram branches vs (micro, macro+bam)
dend_2clusters <- cutree(dend, k=2)
fisher.test(dend_2clusters, c_annotations_group, alternative="g", workspace=2e8)

ht_list <- Heatmap(selected_genes_expression, name = "spectrum", top_annotation = 
                       HeatmapAnnotation(condition = c_annotations_condition, group = c_annotations_group,
                                         col = list(condition = c("ctrl" = "lawngreen", "tumor" = "orangered"), 
                                                    group = c("micro" = "#98ECFE", "macro" = "#FABF00")),
                                         simple_anno_size = unit(3, "cm"), 
                                         annotation_legend_param = list(title_gp = gpar(fontsize = 15), 
                                                                        labels_gp = gpar(fontsize = 15), 
                                                                        legend_height = unit(2, "cm"))
                       ),
                   cluster_rows = T, cluster_columns = dend, show_column_names = F, 
                   row_names_gp = gpar(fontsize = 30), col = col_fun, column_dend_height = unit(3.5, "cm"), 
                   show_row_dend = F,
                   heatmap_legend_param = list(col_fun = col_fun, title = "Expression", 
                                               legend_height = unit(2, "cm"),
                                               title_gp = gpar(fontsize = 15), 
                                               labels_gp = gpar(fontsize = 15),  
                                               at = c(0, 1), labels = c("Low", "High"))
)

pdf(paste0("plots/", "spectrum_heatmap.pdf"), width = 70, height = 12)
draw(ht_list)
decorate_heatmap_body("spectrum", {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
})
dev.off()

# Figure 4a
# left panel
# DimPlot with Hom-MG, Act-MG and MoM
DimPlot(seu_object, group.by = "cell_type_4_groups")

# Figure 4b
# scatter plot
genes <- unique(c(markers_ActMG_vs_HomMG$gene, markers_MoM_vs_ActMG$gene))
gene_expression_data <- GetAssayData(object = seu_object, slot = "data")
gene_expression_data <- as.data.frame(t(gene_expression_data[genes, ]))

common_genes <- intersect(markers_ActMG_vs_HomMG$gene, markers_MoM_vs_ActMG$gene)
markers_ActMG_only <- markers_ActMG_vs_HomMG$gene[!(markers_ActMG_vs_HomMG %in% common_genes)]
markers_MoM_only <- markers_MoM_vs_ActMG$gene[!(markers_MoM_vs_ActMG$gene %in% common_genes)]

genes <- unique(c(markers_ActMG_vs_HomMG$gene, markers_MoM_vs_ActMG$gene))
cell_types_selected <- c("h.Microglia", "a.Microglia", "Macrophages")
names(cell_types_selected) <- c("h.Microglia", "a.Microglia", "Macrophages")
genes_mean_expr <- sapply(cell_types_selected, function(cell_type) {
    Matrix::rowMeans(seu_object@assays$RNA@data[genes, 
                                                colnames(seu_object)[seu_object$cell_type_4_groups == cell_type]],
                     na.rm = T)
})

colnames(genes_mean_expr)[1:3]<-c("Hom_MG", "Act_MG", "MoMphi")
genes_mean_expr <- as.data.frame(genes_mean_expr)
genes_mean_expr$gene <- rownames(genes_mean_expr)

#load genes to labeled 
labels<-c("Ly6c2", "Ccl5", "Ly6i", "Lyz2",   "Lgals3","Ifitm2", 
          "Tgfbi","Tmsb10", "Il1rn","Ass1","Ifitm3","Il1b", 
          "Irf7","Ccr2", "H2-Aa","H2-Ab1", "H2-Eb1","Cd74","Ifit3",        
          "Il18bp", "Mif", "Apoe", "Stat1", "Ccl12",  "H2-D1",  "Ccl4", "Ccl3", "Ly86")

genes_mean_expr$color <- 0
genes_mean_expr[genes_mean_expr$gene %in% markers_ActMG_only, "color"]<-"Act-MG"
genes_mean_expr[genes_mean_expr$gene %in% markers_MoM_only, "color"]<-"MoM"
genes_mean_expr[genes_mean_expr$gene %in% common_genes, "color"]<-"common"
genes_mean_expr$color<-factor(genes_mean_expr$color, levels=c("Act-MG", "MoM", "common"))

rownames(genes_mean_expr) <- annot[match(rownames(genes_mean_expr), annot$Gene.stable.ID), "Gene.name"]
genes_mean_expr$gene <- rownames(genes_mean_expr)
genes_mean_expr_labeled <- genes_mean_expr[labels, ]

ggplot(genes_mean_expr, aes(x=Act_MG, y=MoMphi))+
    geom_jitter(aes(fill=color),shape=21, color="white", alpha=0.7, size=4)+
    geom_text_repel(data=genes_mean_expr_labeled, aes(label=gene), nudge_y=0.2, size=5,
                    direction="both")+
    geom_abline(intercept=0, slope=1)+
    scale_fill_manual(values=c(col_ActMG, col_Macro, "black"))+ 
    xlim(0,5)+
    ylim(0,5)+
    xlab("Act-MG")+
    ylab("MoMphi")+
    coord_fixed()+
    theme_bw(base_size = 18)+
    theme(panel.grid = element_blank())

# Figure 4c
# Heatmap with expression
genes_for_heatmap <- unique(c(markers_MoM_vs_ActMG$Gene.name[1:25], markers_ActMG_vs_HomMG$Gene.name[1:25]))
genes_mean_expr_heatmap <- genes_mean_expr[genes_for_heatmap, ]

mat_breaks <- seq(0, abs(max(genes_mean_expr_heatmap[, 1:3])), length=51)

library(pheatmap)
pheatmap(
    mat               = genes_mean_expr_heatmap[, 1:3],
    color             = colorRampPalette(rev(c("#810f7c", "#8856a7", "#8c96c6", "#b3cde3", "#edf8fb")))(length(mat_breaks) - 1),
    breaks            = mat_breaks,
    border_color      = NA,
    cluster_cols      = F,
    cluster_rows      = F,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    treeheight_col    = 0,
    treeheight_row    = 0,
    gaps_row          = 25,
    drop_levels       = TRUE,
    fontsize          = 10,
    angle_col         = 0,
    main              = "Top upregulated genes in Act-MG and MoM"
)

# Figure 4d
# GO plot
cnetplot(GO_markers_ActMG_MoM_simplify$ActMG_vs_HomMG, foldChange=gene_lists$ActMG_vs_HomMG)

# Figure 4e
# GO plot
cnetplot(GO_markers_ActMG_MoM_simplify$MoM_vs_ActMG, foldChange=gene_lists$MoM_vs_ActMG)

# Figure 4f, g
# Violin Plots with expression of selected genes in Hom-MG, Act-MG and MoM cells
# 4f: Cd52, Usp18, Isg15, Stat1, Ifitm3
# 4g: Il1b, Il1rn, Il18, Il18bp, Cd274
gene_to_plot <- "Cd274"
VlnPlot(seu_object, features = annot[annot$Gene.name == gene_to_plot, "Gene.stable.ID"], 
        idents = c("h.Microglia", "a.Microglia", "Macrophages"), pt.size = 0) + 
        ggtitle(gene_to_plot)

# Figure 4h, i
# plots generated using scSVA

# Figure 5a
# DimPlot coloured by sex-condition value (f_ctrl, f_tumor, m_ctrl, m_tumor)
DimPlot(seu_object, group.by = "sex_condition")

# Figure 5b
# Volvano Plots for Act-MG and MoM markers
library(ggrepel)
col_male<- "#4FCAFF"
col_female<- "#FFFF00"

object_tmp <- seu_object
object_tmp$sex_cell_type <- paste0(seu_object$sex, "_", as.character(Idents(seu_object)))
Idents(object_tmp) <- object_tmp$sex_cell_type

markers_male_ActMG <- FindMarkers(object = object_tmp, ident.1 = "m_a.Microglia", 
                                  ident.2 = "f_a.Microglia", only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.0001)

markers_female_ActMG <- FindMarkers(object = object_tmp, ident.1 = "f_a.Microglia", 
                                    ident.2 = "m_a.Microglia", only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.0001)
markers_female_ActMG$avg_logFC <- markers_female_ActMG$avg_logFC*-1

markers_male_MoM <- FindMarkers(object = object_tmp, ident.1 = "m_Macrophages", 
                                ident.2 = "f_Macrophages", only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.0001)

markers_female_MoM <- FindMarkers(object = object_tmp, ident.1 = "f_Macrophages", 
                                  ident.2 = "m_Macrophages", only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.0001)
markers_female_MoM$avg_logFC <- markers_female_MoM$avg_logFC*-1

markers_list <- list(markers_male_ActMG, markers_female_ActMG, markers_male_MoM, markers_female_MoM)
names(markers_list) <- c("markers_male_ActMG", "markers_female_ActMG", "markers_male_MoM", "markers_female_MoM")
markers_list <- lapply(markers_list, function(x) {
    x$gene <- rownames(x)
    x <- addGeneNames(x)
    x$padj_log10<- -log10(x$p_val_adj)
    x[x$padj_log10 == Inf, "padj_log10"]<-284
    x$log2FC <- log2(exp(x$avg_logFC))
    x <- x %>% mutate(threshold = c(p_val_adj < 0.05 & abs(log2FC) > 0.5 & (pct.1>=0.25 | pct.2 >=0.25)))
    x$color <- x$threshold
    x[x$color == F, "color"] <- "below threshold"
    x[x$color == T, "color"]<-"male"
    x[x$color=="male" & x$log2FC<0, "color"]<-"female"
    x$color<-factor(x$color, levels=c("female", "male", "below threshold"))
    x<-na.omit(x)
    x <-x[rev(order(x$color)),]
    rownames(x)<-seq_along(x$p_val)
    x
})

markers_ActMG_sex <- rbind(cbind(markers_list$markers_male_ActMG, sex="m"),
                           cbind(markers_list$markers_female_ActMG, sex = "f"))

markers_MoM_sex <- rbind(cbind(markers_list$markers_male_MoM, sex="m"),
                         cbind(markers_list$markers_female_MoM, sex = "f"))
markers_MoM_sex[markers_MoM_sex$padj_log10 > 149, "padj_log10"] <- 149

markers_ActMG_sex$color <- factor(markers_ActMG_sex$color, levels=c("female", "male", "below threshold"))
markers_MoM_sex$color <- factor(markers_MoM_sex$color, levels=c("female", "male", "below threshold"))

MG<-ggplot(markers_ActMG_sex, aes(x=log2FC, y=padj_log10, fill=color))+
    geom_jitter(size=3, shape=21, color="black", stroke=0.01)+
    geom_text_repel(data=markers_ActMG_sex[markers_ActMG_sex$threshold==T & markers_ActMG_sex$log2FC>0,],
                    aes(label=Gene.name), nudge_x = 2, nudge_y = 2,
                    color="black", size=6, force=5)+
    geom_text_repel(data=markers_ActMG_sex[markers_ActMG_sex$threshold==T & markers_ActMG_sex$log2FC<0,],
                    aes(label=Gene.name), nudge_x = -2, nudge_y = 2,
                    color="black", size=6, force=5)+
    scale_fill_manual(values=c(col_female,col_male, "gray90"))+
    scale_x_continuous(limits=c(-4,4))+
    labs(title= "Act-MG",x="avg log2FC", y="log10 padj")+
    theme_light(base_size=18)+
    theme(panel.grid = element_blank(), legend.title = element_blank(),
          legend.text = element_text(size=18))+
    guides(fill = guide_legend(override.aes = list(size=6)))

MoM<-ggplot(markers_MoM_sex, aes(x=log2FC, y=padj_log10, fill=color))+
    geom_jitter(size=3, shape=21, color="black")+
    geom_text_repel(data=markers_MoM_sex[markers_MoM_sex$threshold==T & markers_MoM_sex$log2FC>0 & markers_MoM_sex$Gene.name!="AB124611",],
                    aes(label=Gene.name), nudge_x = 2, nudge_y = 5,
                    color="black", size=6)+
    geom_text_repel(data=markers_MoM_sex[markers_MoM_sex$threshold==T & markers_MoM_sex$log2FC<0,],
                    aes(label=Gene.name), nudge_x = -2.5, nudge_y = 2,
                    color="black", size=6)+
    scale_fill_manual(values=c(col_female,col_male, "gray90"))+
    scale_x_continuous(limits=c(-4,4))+
    labs(title= "Mo/M ", x="avg log2FC", y="log10 padj")+
    theme_light(base_size=18)+
    theme(panel.grid = element_blank(), legend.title = element_blank(),
          legend.text = element_text(size=18))+
    guides(fill = guide_legend(override.aes = list(size=6)))

ggarrange(MG, MoM, ncol=2, common.legend = T)


# Figure 5c
# FeaturePlots for selected MHC II genes
# genes: Cd74, H2-Ab1, H2-Eb1, H2-Aa
gene_to_plot <- "Cd74"
FeaturePlot(seu_object, features = annot[annot$Gene.name == gene_to_plot, "Gene.stable.ID"]) + 
    scale_color_gradientn(colors=viridis::viridis(n=50)) +
    ggtitle(gene_to_plot)

# Figure 5d
# Density plots for the same selected genes as in 5c
# Divided by assigned cell type and sex
library(ggridges)
genes<-c("Cd74", "H2-Ab1", "H2-Eb1", "H2-Aa" )

gene_expression_data <- GetAssayData(object = seu_object, slot = "data")
gene_expression_data <- as.data.frame(t(gene_expression_data[annot[match(genes, annot$Gene.name),
                                                                   "Gene.stable.ID"], ]))
colnames(gene_expression_data) <- genes
gene_expression_data$cell_types_4_groups <- as.character(Idents(seu_object))
gene_expression_data$clusters_0.6 <- seu_object$RNA_snn_res.0.6

gene_expression_data$cell_types_6_groups <- gene_expression_data$cell_types_4_groups
gene_expression_data$cell_types_6_groups[gene_expression_data$clusters_0.6 %in% c(10,12)] <- "Mo"
gene_expression_data$cell_types_6_groups[gene_expression_data$clusters_0.6 == 2] <- "intMoM"
gene_expression_data$cell_types_6_groups[gene_expression_data$clusters_0.6 == 13] <- "M"
gene_expression_data$cell_types_6_groups <- factor(gene_expression_data$cell_types_6_groups)
gene_expression_data$cell_types_6_groups <- factor(gene_expression_data$cell_types_6_groups, 
                                                   levels = levels(gene_expression_data$cell_types_6_groups)[c(1,2,6,3,4,5)])

gene_expression_data$sex <- seu_object$gender

graphs<-list() 
for (i in seq_along(genes)){ 
    plot<-ggplot(gene_expression_data, aes(x="", y=cell_types_6_groups, fill=sex))+
        geom_density_ridges(aes_string(x=as.name(genes[i])), alpha=0.65, scale=1 )+
        scale_fill_manual(values=c( "#ffff00",  "#4FCAFF"))+
        xlab("expression level")+
        ylab("")+
        labs(title=genes[i])+
        theme_classic(base_size=18)+
        theme(axis.text = element_text(size=16), axis.text.y = element_text(face="bold", size=18),
              axis.title=element_text(size=16),
              title = element_text(size = 20), legend.text = element_text(size=18),
              legend.title = element_blank())
    graphs[[i]]<-plot
}
ggarrange(plotlist=graphs, ncol=4, nrow=1, common.legend = T, legend="top")

# Figure 5e
# upper panel
# Density plots for MHC II score (H2-Ab1, H2-Eb1 and H2-Aa genes) in intMoM and Act-MG
gene_expression_data$mhc2_score <- rowMeans(gene_expression_data[, c("H2-Ab1", "H2-Eb1", "H2-Aa")])
gene_expression_data <- gene_expression_data  %>% mutate(threshold = mhc2_score > 2)
gene_expression_data[gene_expression_data$threshold==T, "threshold"]<-"MHCII High"
gene_expression_data[gene_expression_data$threshold==F, "threshold"]<-"MHCII Low"
gene_expression_data$threshold <- factor(gene_expression_data$threshold, 
                                         levels = c("MHCII Low", "MHCII High"))

# add Ciita and Mif expression
gene_expression_data_all_genes <- GetAssayData(object = seu_object, slot = "data")
gene_expression_data_ciita_mif <- as.data.frame(t(gene_expression_data_all_genes[annot[match(c("Ciita", "Mif"),
                                                                                             annot$Gene.name), "Gene.stable.ID"], ]))
colnames(gene_expression_data_ciita_mif) <- c("Ciita", "Mif")
gene_expression_data <- cbind(gene_expression_data, gene_expression_data_ciita_mif)

gene_expression_data <- gene_expression_data[gene_expression_data$cell_types_6_groups %in% 
                                                 c("a.Microglia", "intMoM"), ]
gene_expression_data$cell_types_selected <- factor(gene_expression_data$cell_types_6_groups, 
                                                   levels = c("a.Microglia", "intMoM"))

ggplot(gene_expression_data, aes(x=mhc2_score, y=cell_types_selected, fill=sex))+
    geom_density_ridges(alpha=0.65, scale=1 )+
    scale_fill_manual(values=c( "#ffff00",  "#4FCAFF"))+
    geom_vline(xintercept=2, color="red", linetype="dashed", size=1)+
    xlab("MHCII score (H2-Ab1, H2-Eb1, H2-Aa)")+
    ylab("")+
    theme_classic(base_size=18)+
    theme(axis.text.y = element_text(face="bold"),
          title = element_text(size = 20), axis.title.y = element_text(size=16),
          legend.text = element_text(size=18), legend.title = element_blank(),
          axis.text = element_text(size=18), legend.position = "top")

# lower panel
# Violin plots for Mif expression in Act-MG and MoM cells, divided by MHC II score value (low/high)
ggplot(gene_expression_data, aes(x=cell_types_selected, y=Mif, fill=sex))+
    geom_violin(scale="count")+
    scale_fill_manual(values=c( "#ffff00",  "#4FCAFF"))+
    facet_grid(~threshold)+
    ylab("Mif expression level")+
    xlab("")+
    theme_classic(base_size=18)+
    theme(axis.text = element_text(size=18), legend.position = "none")

# Figure 5f
# JM
