#####scRNA-seq GBM LncRNA PAPER

###First verify the number of pairs available
library(openxlsx)
library(readxl)
library(tidyverse)
pairs_meta<-read.xlsx("C:\\Users\\matti\\OneDrive\\Desktop\\GBM_lncRNA_PAPER\\scRNA-seq\\43018_2022_475_MOESM2_ESM.xlsx", sheet = 2)
#pairs_meta<-read.csv("C:\\Users\\matti\\OneDrive\\Desktop\\GBM_lncRNA_PAPER\\scRNA-seq\\EGAD00001008811_metadata\\samples.csv", header = TRUE)
samples_available<-read.delim("C:\\Users\\matti\\OneDrive\\Desktop\\GBM_lncRNA_PAPER\\scRNA-seq\\samples_names.txt", header = FALSE)

#cleaning of samples name (v2 at the end)
#NB: SF10433v2 in the metadata seems to be in reality SF12408
samples_available$sample_clean<-sapply(samples_available$V1, function(x){
  if(isTRUE(grepl("SF10433v2", x))){
    x<-"SF12408"
  } else if (isTRUE(grepl("v2", x)))  {
    x<-gsub("v2", "", x)
  }
  x
})
table(pairs_meta$ID %in% samples_available$sample_clean)

pairs_meta<-pairs_meta %>% filter(ID %in% samples_available$sample_clean)
table(pairs_meta$`Pair#`)

#only paired samples (we are not going to analyze only these, but for DE between pseudobulk groups yes)
pairs_index<-pairs_meta %>% group_by(`Pair#`) %>% count() %>% filter(n==2) %>% pull(`Pair#`)
pairs_only<-pairs_meta[pairs_meta$`Pair#` %in% pairs_index,]
write.csv(pairs_only, "C:\\Users\\matti\\OneDrive\\Desktop\\GBM_lncRNA_PAPER\\scRNA-seq\\pairs_only_metadata.csv")

###Perfect, now define the functions that allows you to order the data properly (two distinct and parallel analysis, on one side you have primary on the other recurrent)

##Load the matrices from STARsolo
export R_LIBS_USER=/data/mpala/R_libs   
export TMPDIR=/data/mpala/R_tmp
export TEMP=/data/mpala/R_tmp
export TMP=/data/mpala/R_tmp
export R_LIBS_USER=/data/mpala/R_libs
conda activate /data/mpala/conda_envs/rnaseq

##NB: later on when you have to call the CNVs remember that the tool is stored into another conda env, so remember to save the object first.
library(SoupX)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(sctransform)
library(cowplot)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(harmony)
library(scales)

###LOADING DATA (STARsolo input)

##NB:for adjusted counts, this is going to be done eventually later, once the clustering has been performed, since most of the correcting ambient RNA tools
##perform better if the clustering information is supplied
dirs<-list.dirs(full.names = TRUE, recursive = TRUE)
folder<-grep("GeneFull$", dirs, value=TRUE)
names(folder)<-sapply(basename(dirname(dirname(dirname(folder)))), function(x){
  if (isTRUE(grepl("SF10433v2", x))){
    x<-"SF12408"
  } else if (isTRUE(grepl("v2", x))){
    x<-gsub("v2", "",x)
  }
  x
})

folder<-folder[names(folder) %in% meta$ID]

read_input_starsolo<-function(dirs, ismultimapping=FALSE){
  sc_matrix<-setNames(lapply(dirs, function(x){
    sample_name<-sapply(basename(dirname(dirname(dirname(x)))), function(y){
      if (isTRUE(grepl("SF10433v2", y))){
        y<-"SF12408"
      } else if (isTRUE(grepl("v2", y))){
        y<-gsub("v2", "",y)
      }
      y
    })
    cat("Loading data for sample:",sample_name, "\n")
    features_path<-paste0(x,"/filtered/features.tsv")
    mtx_path<-paste0(x,"/filtered/matrix.mtx")
    barcodes_path<-paste0(x,"/filtered/barcodes.tsv")
    data<-ReadMtx(mtx=mtx_path, features=features_path, cells=barcodes_path, feature.column=2)
    if (isTRUE(ismultimapping)){
      mtx_multi<-list.files(x, pattern="UniqueAndMult", full.names=TRUE, recursive=TRUE)
      barcodes_raw<-paste0(x, "/raw/barcodes.tsv")
      #here because the counts including unique and multimapping are within the raw folder
      data_multi<-ReadMtx(mtx=mtx_multi, features=features_path, cells=barcodes_raw, feature.column=2)
      data<-data_multi[, colnames(data_multi) %in% colnames(data)] #we essentially filter the raw matrix having the multimapping counts, with the barcodes of the filtered matrix
    } 
    obj <- CreateSeuratObject(data, project = sample_name)
    obj
  }), names(dirs))
}

obj<-read_input_starsolo(folder, ismultimapping = TRUE)
##assign to each sample,a metadata column representing it being either primary or recurrent
primary_id<-meta %>% dplyr::filter(Stage=="Primary") %>% pull("ID")
obj_up<-setNames(lapply(names(obj), function(sample){
  seu<-obj[[sample]]
  seu[[]]$timepoint<-ifelse(sample %in% primary_id, "Primary", "Recurrent")
  seu
}), names(obj))


###MERGE SEURAT OBJECT (Keeping track of cells of origin)

merging_seurat<-function(seu_objs){
  seu_obj <- seu_objs[sapply(seu_objs, function(x) ncol(x) >= 100)] #initial filtering for samples with less than 100 cells
  obj<-merge(seu_objs[[1]],y=seu_objs[-1], add.cell.ids=names(seu_objs))
  return(obj)
}

obj_merged<-merging_seurat(obj_up)

###REMOVE THE DOUBLETS 
doublets_removal<-function(obj_merged, coreparam){
  BPP <- MulticoreParam(coreparam)
  seu_split<-SplitObject(obj_merged, split.by = "orig.ident")
  seu_split_filt<-lapply(seu_split, function(seu) {
    sce <- as.SingleCellExperiment(seu)
    sce <- scDblFinder(sce, BPPARAM = BPP)
    print(table(sce$scDblFinder.class))
    # porta indietro i risultati come metadata (stessa cell order)
    seu$scDblFinder.class <- sce$scDblFinder.class
    seu$scDblFinder.score <- sce$scDblFinder.score
    seu
  })
  merge(seu_split_filt[[1]], y=seu_split_filt[-1]) 
}


obj_merged_ndb<-doublets_removal(obj_merged, 50)

###PERFORM QC
qc_seurat_filtering <- function(obj_merged_ndb, nFeature_RNA_filter, percent_mito_filter,
                                nCount_RNA_filter_up, nCount_RNA_filter_down,
                                outdir = "QC_plots") {
  
  obj_merged_ndb[["percent_mito"]] <- PercentageFeatureSet(obj_merged_ndb, pattern = "^MT-")
  obj_merged_ndb[["percent_ribo"]] <- PercentageFeatureSet(obj_merged_ndb, pattern = "^RP[LS]")
  
  cat("Standard filtering: nFeature, percent_mito, nCount bounds. Removing doublets to stabilize cutoffs.", "\n",
      "Removing samples with fewer than 100 cells.", "\n")
  
  qc_vars <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")
  
  ## ---------- PRE-QC ----------
  data_plot_preQC <- FetchData(obj_merged_ndb, vars = c(qc_vars, "orig.ident", "scDblFinder.class"))
  data_plot_preQC$stage <- "pre_QC"
  
  # P1: scatter nCount vs nFeature colored by mito
  p1 <- ggplot(data_plot_preQC) +
    geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mito), position = "jitter", size = 0.3) +
    scale_color_gradientn(colors = c("lightgrey", "blue")) +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Pre-QC: counts vs features (log)") +
    cowplot::theme_cowplot()
  
  # P2: mito vs nFeature (dying-cell population, top-left) with threshold lines
  p2 <- ggplot(data_plot_preQC) +
    geom_point(aes(x = nFeature_RNA, y = percent_mito, color = scDblFinder.class), size = 0.3) +
    geom_hline(yintercept = percent_mito_filter, linetype = "dashed", color = "red") +
    geom_vline(xintercept = nFeature_RNA_filter, linetype = "dashed", color = "red") +
    scale_x_log10() +
    labs(title = "Pre-QC: mito vs features (threshold lines)") +
    cowplot::theme_cowplot()
  
  # P3: doublets on the count/feature scatter (should sit at high nCount/nFeature)
  p3 <- ggplot(data_plot_preQC) +
    geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = scDblFinder.class), size = 0.3) +
    scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = c(singlet = "grey70", doublet = "red")) +
    labs(title = "Pre-QC: doublet localization") +
    cowplot::theme_cowplot()
  
  # P4: per-sample violins, pre-QC
  p4 <- VlnPlot(obj_merged_ndb, features = qc_vars, group.by = "orig.ident",
                pt.size = 0, ncol = 2) &
    theme(axis.text.x = element_text(angle = 90, size = 6))
  
  # P5: density of each metric with threshold lines
  thresholds <- data.frame(
    variable = c("nFeature_RNA", "percent_mito", "nCount_RNA", "nCount_RNA"),
    xint = c(nFeature_RNA_filter, percent_mito_filter, nCount_RNA_filter_down, nCount_RNA_filter_up)
  )
  long_pre <- tidyr::pivot_longer(data_plot_preQC, cols = all_of(qc_vars),
                                  names_to = "variable", values_to = "value")
  p5 <- ggplot(long_pre, aes(x = value)) +
    geom_density(fill = "grey80") +
    geom_vline(data = thresholds, aes(xintercept = xint), linetype = "dashed", color = "red") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Pre-QC: metric distributions with thresholds") +
    cowplot::theme_cowplot()
  
  ## ---------- FILTERING ----------
  obj_merged_ndb <- subset(obj_merged_ndb,
                           subset = nFeature_RNA > nFeature_RNA_filter &
                             percent_mito < percent_mito_filter &
                             nCount_RNA < nCount_RNA_filter_up &
                             nCount_RNA > nCount_RNA_filter_down)
  obj_merged_ndb <- subset(obj_merged_ndb, subset = scDblFinder.class == "singlet")
  
  cell_per_sample <- table(obj_merged_ndb$orig.ident)
  keep_samples <- names(cell_per_sample[cell_per_sample >= 100])
  obj_merged_ndb <- subset(obj_merged_ndb, subset = orig.ident %in% keep_samples)
  
  ## ---------- POST-QC ----------
  data_plot_postQC <- FetchData(obj_merged_ndb, vars = c(qc_vars, "orig.ident"))
  data_plot_postQC$stage <- "post_QC"
  
  # P6: scatter post-QC
  p6 <- ggplot(data_plot_postQC) +
    geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent_mito), position = "jitter", size = 0.3) +
    scale_color_gradientn(colors = c("lightgrey", "blue")) +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Post-QC: counts vs features (log)") +
    cowplot::theme_cowplot()
  
  # P7: per-sample violins, post-QC
  p7 <- VlnPlot(obj_merged_ndb, features = qc_vars, group.by = "orig.ident",
                pt.size = 0, ncol = 2) &
    theme(axis.text.x = element_text(angle = 90, size = 6))
  
  # P8: cells per sample pre vs post (paired-design safety check)
  pre_counts  <- as.data.frame(table(data_plot_preQC$orig.ident));  pre_counts$stage  <- "pre_QC"
  post_counts <- as.data.frame(table(data_plot_postQC$orig.ident)); post_counts$stage <- "post_QC"
  cell_counts <- rbind(pre_counts, post_counts)
  colnames(cell_counts) <- c("sample", "n_cells", "stage")
  cell_counts$stage <- factor(cell_counts$stage, levels = c("pre_QC", "post_QC"))
  p8 <- ggplot(cell_counts, aes(x = sample, y = n_cells, fill = stage)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
    labs(title = "Cells per sample: pre vs post QC (100-cell cutoff)") +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, size = 6))
  
  ## ---------- SAVE PLOTS (PNG) ----------
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(file.path(outdir, "01_scatter_pre.png"),      p1, width = 7,  height = 6, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "02_mito_vs_feature.png"),  p2, width = 7,  height = 6, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "03_doublets.png"),         p3, width = 7,  height = 6, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "04_violin_pre.png"),       p4, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "05_density_pre.png"),      p5, width = 9,  height = 7, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "06_scatter_post.png"),     p6, width = 7,  height = 6, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "07_violin_post.png"),      p7, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "08_cells_per_sample.png"), p8, width = 12, height = 6, dpi = 300, bg = "white")
  
  cat("QC plots saved to:", normalizePath(outdir), "\n")
  
  return(obj_merged_ndb)
}

#using Inf and 0 for nCounts is like we are not filtering by this parameter(which we do since even the original paper did not)
obj_postqc<-qc_seurat_filtering(obj_merged_ndb=obj_merged_ndb, nFeature_RNA_filter = 200, percent_mito_filter = 2.5, nCount_RNA_filter_up = Inf, nCount_RNA_filter_down = 0, outdir = "QC_plots")

###After the filtering you must verify that all the paired samples are mantained. If not, thus only one sample of the paired is mantained
###we are going to remove it.
samples<-sapply(Layers(obj_postqc), function(x){strsplit(x, "\\.")[[1]][[2]]})
meta_samples<-meta[meta$ID %in% samples,]
table(meta_samples$`Pair.`) #check which pair has been removed and remove it from the object
remove<-meta_samples %>% dplyr::filter(`Pair.` %in% c("18", "19", "30")) %>% pull("ID")
obj_postqc <- subset(obj_postqc, subset = orig.ident %in% remove, invert = TRUE)

####NORMALIZATION , CELL CYCLE METRICS and PCA

SPE <- function(seu_list, outdir = "SPE_plots") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  seu_tmp <- seu_list
  DefaultAssay(seu_tmp) <- "RNA"
  seu_tmp<-JoinLayers(seu_tmp) #joinlayers so that you can compute cell cycle score using the default normalization method
  seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE)
  seu_tmp <- CellCycleScoring(seu_tmp, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  # carry cell-cycle metadata back to the original object
  md_cc <- seu_tmp@meta.data[, c("S.Score", "G2M.Score", "Phase"), drop = FALSE]
  md_cc <- md_cc[colnames(seu_list), , drop = FALSE]
  seu_list <- AddMetaData(seu_list, metadata = md_cc)
  # log-transformed nCount for correlation and perform SCT for normalization
  seu_tmp <- SCTransform(seu_tmp, verbose = FALSE)
  DefaultAssay(seu_tmp) <- "SCT"
  seu_tmp$log_nCount <- log1p(seu_tmp$nCount_RNA)
  seu_tmp <- RunPCA(seu_tmp, verbose = FALSE)
  # P1: PCA by sample
  ggsave(file.path(outdir, "01_pca_by_sample.png"),
         DimPlot(seu_tmp, reduction = "pca", group.by = "orig.ident"),
         width = 8, height = 6, dpi = 300, bg = "white")
  # P2: elbow
  ggsave(file.path(outdir, "02_elbowplot.png"),
         ElbowPlot(seu_tmp, ndims = 50),
         width = 8, height = 6, dpi = 300, bg = "white")
  # P3: PCA regression check (phase, nCount, mito)
  p_pca <- wrap_plots(
    DimPlot(seu_tmp, reduction = "pca", group.by = "Phase"), FeaturePlot(seu_tmp, "nCount_RNA", reduction = "pca"),
    FeaturePlot(seu_tmp, "percent_mito", reduction = "pca"))
  ggsave(file.path(outdir, "03_pca_regression_check.png"), p_pca, width = 18, height = 6, dpi = 300, bg = "white")
  seu_tmp <- RunUMAP(seu_tmp, dims = 1:15, verbose = FALSE)
  # P4: UMAP regression check (phase, nCount, mito, ribo, nFeature)
  p_umap <- wrap_plots(
    DimPlot(seu_tmp, reduction = "umap", group.by = "Phase", pt.size = 0.7),
    FeaturePlot(seu_tmp, "nCount_RNA", reduction = "umap", pt.size = 0.7),
    FeaturePlot(seu_tmp, "percent_mito", reduction = "umap", pt.size = 0.7),
    FeaturePlot(seu_tmp, "percent_ribo", reduction = "umap", pt.size = 0.7),
    FeaturePlot(seu_tmp, "nFeature_RNA", reduction = "umap", pt.size = 0.7))
  ggsave(file.path(outdir, "04_umap_regression_check.png"), p_umap, width = 22, height = 12, dpi = 300, bg = "white")
  # PC vs covariate correlations (diagnostic for what to regress)
  emb <- Embeddings(seu_tmp, "pca")[, 1:5]
  pca_cor_mito   <- apply(emb, 2, function(x) cor(x, seu_tmp$percent_mito))
  pca_cor_nCount <- apply(emb, 2, function(x) cor(x, seu_tmp$log_nCount))
  pca_cor_s      <- apply(emb, 2, function(x) cor(x, seu_tmp$S.Score))
  pca_cor_g2m    <- apply(emb, 2, function(x) cor(x, seu_tmp$G2M.Score))
  pca_cor_ribo   <- apply(emb, 2, function(x) cor(x, seu_tmp$percent_ribo))
  cat("Correlation of first 5 PCs with candidate covariates (to decide what to regress out):\n")
  cat("percent_mito:", round(pca_cor_mito, 3), "\n")
  cat("log_nCount  :", round(pca_cor_nCount, 3), "\n")
  cat("S.Score     :", round(pca_cor_s, 3), "\n")
  cat("G2M.Score   :", round(pca_cor_g2m, 3), "\n")
  cat("percent_ribo:", round(pca_cor_ribo, 3), "\n")
  
  # optional: save the correlation table too
  cor_tab <- rbind(percent_mito = pca_cor_mito, log_nCount = pca_cor_nCount,
                   S.Score = pca_cor_s, G2M.Score = pca_cor_g2m, percent_ribo = pca_cor_ribo)
  write.csv(cor_tab, file.path(outdir, "05_pc_covariate_correlations.csv"))
  
  cat("SPE diagnostic plots saved to:", normalizePath(outdir), "\n")
  
  return(seu_list)
}


objs_sct_initial<-SPE(obj_postqc) 
saveRDS(objs_sct_initial, "toy.rds")
####InferCNV to call and distinguish between tumor and normal cells. 
###First you need to set reference cells (non tumor) thus you have to run a first round of clustering and cell type detection based only on the non-tumor cell markers.
###Here you don't need integration , all you need is to identify cells that show high expression of these cells. (is true that the integration would help to localize better these cells in the UMAP though)

get_normal_cells<-function(seu_list, dimensions, regress_out=NULL, integration_method, resolution_parameter){
  integration_name<-deparse(substitute(integration_method))
  seu_list<-SCTransform(seu_list, vars.to.regress=regress_out, verbose=FALSE)
  DefaultAssay(seu_list) <- "SCT"
  seu_list<-RunPCA(seu_list)
  cat("Start Integration","\n")
  seu_list_int<-IntegrateLayers(seu_list, method=integration_method, normalization.method="SCT", verbose=F, new.reduction="integrated.dr")
  for(d in dimensions){
    seu_list_int<-RunUMAP(seu_list_int, dims=1:d, reduction="pca", reduction.name=paste0("umap.pca.d_",d)) 
    seu_list_int<-RunUMAP(seu_list_int, dims=1:d, reduction="integrated.dr", reduction.name=paste0("umap.integrated.d_",d))
    seu_list_int<-FindNeighbors(seu_list_int, reduction="integrated.dr", dims=1:d)
    seu_list_int<-FindClusters(seu_list_int, resolution=resolution_parameter, cluster.name=paste0("integrated.cluster.d_", d))
    #compute the module Scores
    normal_cells_genes<-list("immune_score"=c("PTPRC", "MRC1", "TMEM119", "P2RY12", "CSF1R", "AIF1"),"oligo_score"=c("MOG", "MBP", "MAG", "PLP1"), "endothelial_score"=c("VWF", "CLDN5", "PDGFRB", "RGS5"))
    for (set in names(normal_cells_genes)) {
      seu_list_int <- AddModuleScore(object= seu_list_int, features  = list(normal_cells_genes[[set]]),ctrl= 100, nbin= 12, name= set)
    }
    cl_col <- paste0("integrated.cluster.d_", d)
    qc_vars <- c("immune_score1", "oligo_score1", "endothelial_score1")
    md <- seu_list_int@meta.data
    qc_vars <- qc_vars[qc_vars %in% colnames(md)]
    qc_medians <- aggregate( md[, qc_vars, drop = FALSE], by = list(cluster = md[[cl_col]]), FUN = median, na.rm = TRUE)
    # Add cluster sizes
    cl_sizes <- as.data.frame(table(md[[cl_col]]))
    colnames(cl_sizes) <- c("cluster", "n_cells")
    qc_medians$cluster <- as.character(qc_medians$cluster)
    cl_sizes$cluster <- as.character(cl_sizes$cluster)
    qc_summary <- merge(cl_sizes, qc_medians, by = "cluster", all.x = TRUE)
    qc_summary <- qc_summary[order(as.integer(qc_summary$cluster)), ]
    ## ---------- PLOTS: module scores, key markers, clusters ----------
    plot_dir <- paste0("normal_cells_plots/d_", d)
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    umap_int <- paste0("umap.integrated.d_", d)
    score_cols <- qc_vars
    # 1) UMAP colored by cluster
    p_clusters <- DimPlot(seu_list_int, reduction = umap_int, group.by = cl_col, label = TRUE, repel = TRUE) + ggtitle(paste0("Clusters (integrated, d=", d, ")"))
    ggsave(file.path(plot_dir, "01_umap_clusters.png"), p_clusters, width = 9, height = 7, dpi = 300, bg = "white")
    # 2) UMAP FeaturePlot of module scores (one panel per score)
    p_scores <- FeaturePlot(seu_list_int, features = score_cols, reduction = umap_int,order = TRUE, ncol = length(score_cols)) &scale_color_gradientn(colors = c("lightgrey", "red"))
    ggsave(file.path(plot_dir, "02_umap_module_scores.png"), p_scores, width = 7 * length(score_cols), height = 6, dpi = 300, bg = "white")
    # 3) Violin of module scores per cluster (where each normal type sits)
    p_score_vln <- VlnPlot(seu_list_int, features = score_cols, group.by = cl_col, pt.size = 0, ncol = 1)
    ggsave(file.path(plot_dir, "03_violin_module_scores_by_cluster.png"), p_score_vln, width = 12, height = 4 * length(score_cols), dpi = 300, bg = "white")
    # 4) DotPlot of individual key markers per cluster (validation of scores)
    all_markers <- unique(unlist(normal_cells_genes))
    all_markers <- all_markers[all_markers %in% rownames(seu_list_int)]
    p_dot <- DotPlot(seu_list_int, features = all_markers, group.by = cl_col) + RotatedAxis() + ggtitle(paste0("Normal-cell markers per cluster (d=", d, ")"))
    ggsave(file.path(plot_dir, "04_dotplot_markers_by_cluster.png"), p_dot, width = max(8, 0.5 * length(all_markers) + 4), height = 8, dpi = 300, bg = "white")
    # 5) FeaturePlot of individual key markers on UMAP
    p_markers_umap <- FeaturePlot(seu_list_int, features = all_markers, reduction = umap_int, order = TRUE, ncol = 3) & scale_color_gradientn(colors = c("lightgrey", "blue"))
    ggsave(file.path(plot_dir, "05_umap_key_markers.png"), p_markers_umap, width = 15, height = 5 * ceiling(length(all_markers) / 3), dpi = 300, bg = "white")
    cat("Normal-cell plots for d =", d, "saved to:", normalizePath(plot_dir), "\n")
  }
  return(seu_list_int)
}

##After you visualized the cells, you must decide which cells keep as reference. Do not forget to build the data required for the
obj_normal_including<-get_normal_cells(objs_sct_initial,20,NULL, "HarmonyIntegration", 0.8)
md<-obj_normal_including[[]]
md<-md %>% mutate(cell_type=case_when(integrated.cluster.d_20 %in% c("0", "35","31", "25", "20") ~ "Myeloid", integrated.cluster.d_20 %in% c("16", "29", "37", "1") ~ "Oligodendrocytes", TRUE ~ "rest"))
refcells<-md %>% filter(cell_type %in% c("Myeloid", "Oligodendrocytes")) 
refcells_df<-data.frame(barcodes = rownames(refcells), cell_type=refcells$cell_type)
refCells<-lapply(split(refcells_df, refcells_df$cell_type), function(x){x<-x %>% select(-"cell_type") %>% pull("barcodes")})

####RUN INFERCNA

#NB: you need to set the genome. Check if the default hg38 is fine, if not you need to built it manually.
##IMPORTANT: the size of the matrix is also huge, which means that you may want to either run individually for each sample (with the risk that the reference cells are too low)
##or you may want filtering by 0 expression genes, so you reduce the size
library(infercna)
library(stringr)
library(Seurat)
library(scalop)
library(ggplot2)
useGenome("hg38")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)                              # unisce i layer per-campione in uno
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")   # estrae la MATRICE (sparsa)
keep<-rowSums(counts>0)>10
counts<-counts[keep,]
library_size<-colSums(counts)
cpm <- Matrix::t( Matrix::t(counts) / library_size ) * 1e6
#scale by factor 10 and log transform
cpm_log<-log2((cpm/10)+1)
cpm_log<-as.matrix(cpm_log)
cna = infercna(m = cpm_log, refCells = refCells, n = 5000, noise = 0.1, isLog = TRUE, verbose = FALSE)
cnaM = cna[, !colnames(cna) %in% unlist(refCells)]
p<-ggcna(cnaM, genome="hg38")
ggsave("heatmap_cnv_toy.png", p)

##compute signal and correlation for that
cna_sig<-cnaSignal(cnaM,)
png("QC_plots/cnaScatter_combined.png", width = 1600, height = 700)
par(mfrow = c(1, 2))
cnaScatterPlot(cna = cna, gene.quantile = 0.9, samples = NULL,
               refCells = unlist(refcell))
cnaScatterPlot(cna = cna, gene.quantile = 0.9, samples = NULL,
               refCells = unlist(refcell), groups = unlist(refcell))
par(mfrow = c(1, 1))
dev.off()
