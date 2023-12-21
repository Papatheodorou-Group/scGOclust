# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song ysong@ebi.ac.uk


#' get requested ensembl ID to GO mapping table
#' @name ensemblToGo
#' @param species species name matching ensembl biomaRt naming, such as hsapiens, mmusculus
#' @param GO_type GO term type, choose among 'biological_process', 'molecular_function', 'cellular_component', default 'biological_process'
#' @param GO_linkage_type GO annotation evidence codes to include. Default is 'standard', which means only including manually checked records (excluding IEA) and excluding those inferred from gene expression experiments to maximally suffice the species expression independence assumption. 'Stringent' means only including those with experimental evidence, also not from gene expression experiment, or from manual curation with evidence (excluding those from mass-annotation pipelines). Choose among 'experimental', 'phylogenetic', 'computational', 'author', 'curator', 'electronic', 'standard', stringent'
#' @param ... additional params for useEnsembl function called in this function
#' @return a table with ensembl to GO terms mapping including requested linkage type
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' ensemblToGo("mmusculus", GO_type = "biological_process", GO_linkage_type = 'stringent')
#' }
#' @importFrom biomaRt useMart useEnsembl
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate case_when
#' @export
#'
#'
ensemblToGo <- function(species, GO_type = "biological_process", GO_linkage_type = c("standard"), ...) {

  ## GO source type code
  go_source <- list(
    experimental = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"),
    phylogenetic = c("IBA", "IBD", "IKR", "IRD"),
    computational = c("ISS", "ISO", "ISA", "ISM", "IGC", "RCA"),
    author = c("TAS", "NAS"),
    curator = c("IC", "ND"),
    electronic = c("IEA"),
    standard = c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "IBA", "IBD", "IKR", "IRD", "ISS", "ISO", "ISA", "ISM", "IGC", "RCA", "TAS", "IC"),
    stringent = c("EXP", "IDA", "IPI", "IMP", "IGI",  "HTP", "HDA", "HMP", "HGI", "IKR", "RCA", "TAS", "IC")
  )

  # Get ensembl gene ids and GO terms

  message("start query ENSEMBL via biomaRt")

  bm <- tryCatch(
    {
      useEnsembl(biomart = "ensembl", dataset = paste0(species, "_gene_ensembl"), ...)
    },
    warning = function(w) {
      message("ensembl biomaRt warning (useEnsembl):")
      message(w)
    },
    error = function(e) {
      message("ensembl biomaRt error (useEnsembl):")
      message(e)
    }

  )


  EG2GO <- tryCatch(
    {
      getBM(mart = bm, attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006", "go_linkage_type", "namespace_1003"))
    },
    warning = function(w) {
      message("ensembl biomaRt warning (getBM):")
      message(w)
    },
    error = function(e) {
      message("ensembl biomaRt error (getBM):")
      message(e)
    }
  )



  # Remove blank entries
  EG2GO <- EG2GO[EG2GO$go_id != "", ]

  GO_terms <- EG2GO %>%
    dplyr::filter(namespace_1003 == GO_type) %>%
    dplyr::filter(!(name_1006 %in% c("biological_process", "molecular_function", "cellular_component")))

  types_in_source <- list()

  for (source_type in names(go_source)) {
    type_in <- levels(factor(GO_terms$go_linkage_type))[levels(factor(GO_terms$go_linkage_type)) %in% go_source[[source_type]]]
    types_in_source[[source_type]] <- type_in
  }

  included_terms <- sapply(GO_linkage_type, function(x) paste(go_source[[x]], sep = ', ', collapse = ', '))
  message("including GO link types: ")
  for (name in names(included_terms)){
    message(paste0(name, ": "))
    message(included_terms[name])
  }

  use <- as.character(unname(unlist(sapply(GO_linkage_type, function(x) go_source[[x]]))))

  GO_terms <- GO_terms %>% dplyr::filter(go_linkage_type %in% use)

  ## add source type to GO terms
  GO_terms <- GO_terms %>%
    mutate(evidence_type = case_when(
      go_linkage_type %in% go_source$experimental ~ "experimental",
      go_linkage_type %in% go_source$phylogenetic ~ "phylogenetic",
      go_linkage_type %in% go_source$computational ~ "computational",
      go_linkage_type %in% go_source$author ~ "author",
      go_linkage_type %in% go_source$curator ~ "curator",
      go_linkage_type %in% go_source$electronic ~ "electronic",
      TRUE ~ "unknown"
    ))

  return(GO_terms)
}




#' create a seurat object with GO terms
#' @name makeGOSeurat
#' @param ensembl_to_GO ensembl_to_go mapping table from function ensemblToGo
#' @param seurat_obj count matrix with genes to cells
#' @param feature_type feature type of count matrix, choose from ensembl_gene_id, external_gene_name, default ensembl_gene_id
#' @return a seurat object with GO terms as features
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#' }
#' @importFrom biomaRt useMart useDataset getBM
#' @importFrom Matrix Matrix
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom Seurat CreateSeuratObject
#' @importFrom dplyr filter mutate arrange select
#' @export
#'
#'

makeGOSeurat <- function(ensembl_to_GO, seurat_obj, feature_type = "ensembl_gene_id") {
  message("collect data")
  counts <- as.matrix(seurat_obj@assays$RNA@counts)

  if(!(feature_type %in% colnames(ensembl_to_GO))){
    stop(paste0(feature_type, " is not in colnames(ensembl_to_GO), please check the var type"))
  }

  ## pivot GO to feature type table to matrix
  go_matrix <- ensembl_to_GO %>%
    dplyr::mutate(placehold = 1) %>%
    tidyr::pivot_wider(id_cols = eval(feature_type), names_from = name_1006, values_from = placehold, values_fill = 0, values_fn = dplyr::first)

  shared <- intersect(go_matrix[[feature_type]], rownames(counts))

  go_matrix <- go_matrix %>%
    dplyr::filter(get(feature_type) %in% shared) %>%
    dplyr::arrange(get(feature_type))
  go_matrix <- go_matrix %>% tibble::column_to_rownames(eval(feature_type))

  counts <- t(counts)
  counts <- counts %>% as.data.frame()
  counts <- counts %>% dplyr::select(eval(shared))
  counts <- counts[, match(rownames(go_matrix), colnames(counts))]


  if (!(all(colnames(counts) == rownames(go_matrix)))) {
    stop("error during calculation, please check input format")
  }

  message("compute GO to cell matrix, might take a few secs")
  start_time <- Sys.time()
  go_mtx <- Matrix::Matrix(as.matrix(counts), sparse = TRUE) %*% Matrix::Matrix(as.matrix(go_matrix), sparse = TRUE)
  go_obj <- Seurat::CreateSeuratObject(counts = t(as.matrix(go_mtx)), meta.data = seurat_obj@meta.data)
  end_time <- Sys.time()

  message(paste0("time used: ", round(end_time - start_time, 2), " secs"))

  all_zero_terms = names(rowSums(as.matrix(go_obj@assays$RNA@layers$counts))[which(rowSums(as.matrix(go_obj@assays$RNA@layers$counts)) == 0)])
  message(paste0("removing ", length(all_zero_terms)), " all zero terms")

  go_obj <- go_obj[which(!(rownames(go_obj) %in% all_zero_terms)), ]

  message("returning GO Seurat object")
  return(go_obj)
}

#' standard seurat analysis on GO_seurat object
#' @name analyzeGOSeurat
#' @param go_seurat_obj go seurat object created by makeGOSeurat
#' @param cell_type_col column name in mera.data storing cell type classes
#' @param norm_log1p whether or not to perform data normalisation and log1p transformation, default TRUE
#' @param cluster_res resolution for Seurat FindClusters
#' @param scale.factor param for Seurat NormalizeData
#' @param nfeatures param for Seurat FindVariableFeatures
#' @param min.dist param for Seurat RunUMAP
#' @param ... additional params for all Seurat functions involved in this function
#' @return standard analyzed GO seurat object until UMAP
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' go_seurat_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#'
#' analyzeGOSeurat(go_seurat_obj = go_seurat_obj, cell_type_col = "cell_type_annotation")
#' }
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters RunUMAP Idents
#' @export
#'


analyzeGOSeurat <- function(go_seurat_obj, cell_type_col, norm_log1p=TRUE, scale.factor = 10000, nfeatures = 2000, cluster_res = 1, min.dist=0.3, ...) {
  if (!(cell_type_col %in% colnames(go_seurat_obj@meta.data))) {
    stop("cell_type_col not in annotation, please check input")
  }
  if(norm_log1p){
    message(paste0("perform normalization and log1p for ", deparse(substitute(go_seurat_obj))))
    go_seurat_obj <- Seurat::NormalizeData(go_seurat_obj, normalization.method = "LogNormalize", scale.factor = scale.factor, ...)

  }
  go_seurat_obj <- Seurat::FindVariableFeatures(go_seurat_obj, selection.method = "vst", nfeatures = nfeatures, ...)
  go_seurat_obj <- Seurat::ScaleData(object = go_seurat_obj, features = rownames(go_seurat_obj), verbose = FALSE, ...)
  go_seurat_obj <- Seurat::RunPCA(object = go_seurat_obj, npcs = 50, verbose = FALSE, ...)
  go_seurat_obj <- Seurat::FindNeighbors(go_seurat_obj, dims = 1:50, ...)
  go_seurat_obj <- Seurat::FindClusters(go_seurat_obj, resolution = cluster_res, ...)
  go_seurat_obj <- Seurat::RunUMAP(object = go_seurat_obj, reduction = "pca", min.dist = min.dist, dims = 1:50, ...)
  Seurat::Idents(go_seurat_obj) <- go_seurat_obj@meta.data[[cell_type_col]]
  return(go_seurat_obj)
}

#' get per cell type average scaled vector of GO terms
#' @name getCellTypeGO
#' @param go_seurat_obj go seurat object created by makeGOSeurat
#' @param cell_type_col column name in mera.data storing cell type classes
#' @param norm_log1p whether or not to perform data normalisation and log1p transformation, default TRUE
#' @return a table of scaled GO representation per cell type (averaged)
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' go_seurat_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#' getCellTypeGO(go_seurat_obj = go_seurat_obj, cell_type_co = "cell_type_annotation")
#' }
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData AverageExpression
#' @export
#'

getCellTypeGO <- function(go_seurat_obj, cell_type_col, norm_log1p = TRUE) {
  if (!(cell_type_col %in% colnames(go_seurat_obj@meta.data))) {
    stop("cell_type_col not in metadata, please check input")
  }

  if(norm_log1p){
    message(paste0("perform normalization and log1p for ", deparse(substitute(go_seurat_obj))))
    go_seurat_obj <- NormalizeData(go_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  }


  go_seurat_obj <- ScaleData(object = go_seurat_obj, features = rownames(go_seurat_obj), verbose = TRUE, scale.max = 10000) ## set scale.max to near-unlimited, also use all features

  tbl <- AverageExpression(object = go_seurat_obj, group.by = cell_type_col, layer = "scale.data")
  return(as.data.frame(tbl[["RNA"]]))
}


#' calculate correlation between cell types represented by scaled GO, per-species
#' @name cellTypeGOCorr
#' @param cell_type_go cell type GO table calculated via getCellTypeGO
#' @param corr_method correlation method, choose among "pearson", "kendall", "spearman", default 'pearson'
#' @return a dataframe with correlation between cell types
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' go_seurat_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#'
#' cell_type_go = getCellTypeGO(go_seurat_obj = go_seurat_obj, cell_type_co = "cell_type_annotation")
#'
#' cellTypeGOCorr(cell_type_go = cell_type_go, corr_method = "pearson")
#' }
#'
#' @importFrom stats cor
#' @export
#'


cellTypeGOCorr <- function(cell_type_go, corr_method = "pearson") {
  cell_type_go <- as.data.frame(cell_type_go)
  all_cell_types <- colnames(cell_type_go)
  corr_matrix <- data.frame()

  for (cell_type_1 in all_cell_types) {
    for (cell_type_2 in all_cell_types) {
      corr_val <- stats::cor(cell_type_go[, cell_type_1], cell_type_go[, cell_type_2], method = corr_method)

      corr_matrix[cell_type_1, cell_type_2] <- corr_val
    }
  }

  return(corr_matrix)
}

#' calculate cross-species correlation between cell types represented by scaled GO
#' @name crossSpeciesCellTypeGOCorr
#' @param species_1 name of species one
#' @param species_2 name of species two
#' @param cell_type_go_sp1 cell type GO table of species one calculated via getCellTypeGO
#' @param cell_type_go_sp2 cell type GO table of species two calculated via getCellTypeGO
#' @param corr_method correlation method, choose among "pearson", "kendall", "spearman", default 'pearson'
#' @return correlation between cell types
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' data(dme_tbl)
#' data(dme_subset)
#' mmu_go_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#' dme_go_obj = makeGOSeurat(ensembl_to_GO = dme_tbl,
#'  seurat_obj = dme_subset,
#'  feature_type = "external_gene_name")
#'
#' mmu_cell_type_go = getCellTypeGO(go_seurat_obj = mmu_go_obj, cell_type_co = "cell_type_annotation")
#' dme_cell_type_go = getCellTypeGO(go_seurat_obj = dme_go_obj, cell_type_co = "annotation")
#'
#' crossSpeciesCellTypeGOCorr(species_1 = 'mmusculus',
#'  species_2 = 'dmelanogaster',
#'  cell_type_go_sp1 = mmu_cell_type_go,
#'  cell_type_go_sp2 = dme_cell_type_go)
#' }
#' @importFrom stats cor
#' @export
#'

crossSpeciesCellTypeGOCorr <- function(species_1, species_2, cell_type_go_sp1, cell_type_go_sp2, corr_method = "pearson") {
  cell_type_go_sp1 <- as.data.frame(cell_type_go_sp1)
  cell_type_go_sp2 <- as.data.frame(cell_type_go_sp2)

  all_cell_types_sp1 <- colnames(cell_type_go_sp1)
  all_cell_types_sp2 <- colnames(cell_type_go_sp2)

  ## take intersection of GO terms and match order
  intersection <- intersect(rownames(cell_type_go_sp1), rownames(cell_type_go_sp2))
  cell_type_go_sp1 <- cell_type_go_sp1[intersection, ]
  cell_type_go_sp2 <- cell_type_go_sp2[intersection, ]

  if (!(all(rownames(cell_type_go_sp1) == rownames(cell_type_go_sp2)))) {
    stop("cross-species matching of GO terms failed, please check input")
  }

  corr_matrix <- data.frame()

  for (cell_type_1 in all_cell_types_sp1) {
    for (cell_type_2 in all_cell_types_sp2) {
      corr_val <- stats::cor(cell_type_go_sp1[, cell_type_1], cell_type_go_sp2[, cell_type_2], method = corr_method)

      corr_matrix[paste0(species_1, "_", cell_type_1), paste0(species_2, "_", cell_type_2)] <- corr_val
    }
  }

  return(corr_matrix)
}





#' plot clustered heatmap for cell type corr
#' @name plotCellTypeCorrHeatmap
#' @param corr_matrix correlation matrix from cellTypeGOCorr or crossSpeciesCellTypeGOCorr
#' @param scale scale value by column, row, or default no scaling (NA)
#' @param ... params to pass to slanter::sheatmap
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#'
#' go_seurat_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#'
#' cell_type_go = getCellTypeGO(go_seurat_obj = go_seurat_obj, cell_type_co = "cell_type_annotation")
#'
#' corr_matrix = cellTypeGOCorr(cell_type_go = cell_type_go, corr_method = "pearson")
#'
#' plotCellTypeCorrHeatmap(corr_matrix = corr_matrix, scale = "column")
#'
#' }
#' @return a sheatmap heatmap
#' @importFrom slanter sheatmap
#' @export
#'

plotCellTypeCorrHeatmap <- function(corr_matrix, scale = NA, ...) {

  # round number to the lower 0.05
  round_down_5 <- function(x){
    y = round(x, digits = 2)
    # handle machine precision when y = 0.05
    if((y %% 0.1) > 0.05 + 1e-8){
      return (round(y / 0.1, digits = 0)*0.1 - 0.05)
    } else{
      return (round(y / 0.1, digits = 0)*0.1)
    }

  }

  # round number to the upper 0.05
  round_up_5 <- function(x){
    y = round(x, digits = 2)
    # handle machine precision when y = 0.05
    if((y %% 0.1) > 0.05 + 1e-8){
      return (round(y / 0.1, digits = 0)*0.1)
    } else{
      return (round(y / 0.1, digits = 0)*0.1 + 0.05)
    }

  }


  if(scale %in% c("column", 'row')){

    heatmap <- slanter::sheatmap(corr_matrix + 0.5, scale = scale, ...)

  } else if(is.na(scale)){

    lower_bound = round_down_5(min(corr_matrix + 0.5))
    upper_bound = round_up_5(max(corr_matrix + 0.5))

    breaks = seq(lower_bound, upper_bound, by = 0.05)

    # slanter only takes positive matrix values so we have a work around
    heatmap <- slanter::sheatmap(corr_matrix + 0.5 , legend_breaks = breaks, legend_labels = round(breaks - 0.5, digits = 2), ...)

  } else{

    warning('scale = must be one of NA, row or column')
  }
  return(heatmap)
}



#' get shared up and down regulated GO terms for all pairs of cell types
#' @name getCellTypeSharedGO
#' @param species_1 name of species one
#' @param species_2 name of species two
#' @param analyzed_go_seurat_sp1 analyzed GO seurat object of species one
#' @param analyzed_go_seurat_sp2 analyzed GO seurat object of species two
#' @param cell_type_col_sp1 cell type column name for species 1 data
#' @param cell_type_col_sp2 cell type column name for species 2 data
#' @param layer_use layer to use for marker computation, default 'data' which after NormalizeData will be log1p normalized data.
#' @param p_val_threshould p value threshold for selecting DEG (p_adjust)
#' @return a list with sp1 raw, sp2 raw and shared, significant up and down regulated GO terms per cell type (pair)
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' data(dme_tbl)
#' data(dme_subset)
#'
#' mmu_go_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#' dme_go_obj = makeGOSeurat(ensembl_to_GO = dme_tbl,
#'  seurat_obj = dme_subset,
#'  feature_type = "external_gene_name")
#'
#'
#' mmu_go_obj_analyzed = analyzeGOSeurat(mmu_go_obj, "cell_type_annotation")
#' dme_go_obj_analyzed = analyzeGOSeurat(dme_go_obj, "annotation")
#'
#' getCellTypeSharedGO(species_1 = 'mmusculus',
#' species_2 = 'dmelanogaster',
#' analyzed_go_seurat_sp1 =  mmu_go_obj_analyzed,
#' analyzed_go_seurat_sp2 =  dme_go_obj_analyzed,
#' cell_type_col_sp1 = 'cell_type_annotation',
#' cell_type_col_sp2 = 'annotation',
#' layer_use = "data",
#' p_val_threshould = 0.01)
#' }
#' @importFrom Seurat FindAllMarkers Idents
#' @import limma
#' @importFrom dplyr filter mutate
#' @export
#'


getCellTypeSharedGO <- function(species_1, species_2, analyzed_go_seurat_sp1, analyzed_go_seurat_sp2, cell_type_col_sp1, cell_type_col_sp2, layer_use = "data", p_val_threshould = 0.01) {
  if (!(cell_type_col_sp1 %in% colnames(analyzed_go_seurat_sp1@meta.data))) {
    stop("cell_type_col_sp1 not in annotation, please check input")
  }

  if (!(cell_type_col_sp2 %in% colnames(analyzed_go_seurat_sp2@meta.data))) {
    stop("cell_type_col_sp2 not in annotation, please check input")
  }

  ## set idents for marker calling
  Seurat::Idents(analyzed_go_seurat_sp1) <- analyzed_go_seurat_sp1@meta.data[[cell_type_col_sp1]]
  Seurat::Idents(analyzed_go_seurat_sp2) <- analyzed_go_seurat_sp2@meta.data[[cell_type_col_sp2]]


  message(paste0("calculate cell type marker for species ", species_1, ", this will take a while"))
  sp1_markers <- FindAllMarkers(object = analyzed_go_seurat_sp1, layer = layer_use, test.use = "wilcox", verbose = TRUE, )

  message(paste0("calculate cell type marker for species ", species_2, ", this will take a while"))
  sp2_markers <- FindAllMarkers(object = analyzed_go_seurat_sp2, layer = layer_use, test.use = "wilcox", verbose = TRUE)


  sp1_cts <- levels(factor(analyzed_go_seurat_sp1@meta.data[[cell_type_col_sp1]]))
  sp2_cts <- levels(factor(analyzed_go_seurat_sp2@meta.data[[cell_type_col_sp2]]))

  message("collecting shared up regulated terms")

  shared_all <- data.frame()
  for (ct_sp1 in sp1_cts) {
    for (ct_sp2 in sp2_cts) {
      sp1_sig_up <- sp1_markers %>%
        filter(cluster == ct_sp1) %>%
        filter(avg_log2FC > 0) %>%
        filter(p_val_adj < p_val_threshould) %>%
        mutate(marker_type = "sig_up")

      sp1_sig_up_terms <- sp1_sig_up$gene

      sp2_sig_up <- sp2_markers %>%
        filter(cluster == ct_sp2) %>%
        filter(avg_log2FC > 0) %>%
        filter(p_val_adj < p_val_threshould) %>%
        mutate(marker_type = "sig_up")

      sp2_sig_up_terms <- sp2_sig_up$gene

      intersect <- intersect(sp1_sig_up_terms, sp2_sig_up_terms)

      sp1_sig_up <- sp1_sig_up %>%
        filter(gene %in% intersect) %>%
        mutate(pct_intersect = length(intersect) / length(sp1_sig_up_terms))

      sp2_sig_up <- sp2_sig_up %>%
        filter(gene %in% intersect) %>%
        mutate(pct_intersect = length(intersect) / length(sp2_sig_up_terms))

      shared <- merge(sp1_sig_up, sp2_sig_up, by = 'gene', suffixes = c("_sp1", "_sp2"))

      shared_all <- rbind(shared_all, shared)
    }
  }


  message("collecting shared down regulated terms")
  for (ct_sp1 in sp1_cts) {
    for (ct_sp2 in sp2_cts) {
      sp1_sig_down <- sp1_markers %>%
        filter(cluster == ct_sp1) %>%
        filter(avg_log2FC < 0) %>%
        filter(p_val_adj < p_val_threshould) %>%
        mutate(marker_type = "sig_down")

      sp1_sig_down_terms <- sp1_sig_down$gene

      sp2_sig_down <- sp2_markers %>%
        filter(cluster == ct_sp2) %>%
        filter(avg_log2FC < 0) %>%
        filter(p_val_adj < p_val_threshould) %>%
        mutate(marker_type = "sig_down")

      sp2_sig_down_terms <- sp2_sig_down$gene

      intersect <- intersect(sp1_sig_down_terms, sp2_sig_down_terms)

      sp1_sig_down <- sp1_sig_down %>%
        filter(gene %in% intersect) %>%
        mutate(pct_intersect = length(intersect) / length(sp1_sig_down_terms))

      sp2_sig_down <- sp2_sig_down %>%
        filter(gene %in% intersect) %>%
        mutate(pct_intersect = length(intersect) / length(sp2_sig_down_terms))

      shared <- merge(sp1_sig_down, sp2_sig_down, by = 'gene', suffixes = c("_sp1", "_sp2"))

      shared_all <- rbind(shared_all, shared)
    }
  }
  shared_all$species_1 <- species_1
  shared_all$species_2 <- species_2
  message("finish getting cell type pairs shared up and down regulated GO terms")

  results <- list(sp1_markers_raw = sp1_markers, sp2_markers_raw = sp2_markers, shared_sig_markers = shared_all)
  return(results)
}

#' plot Sankey diagram for cell type links above a certain threshould
#' @name plotCellTypeSankey
#' @param corr_matrix cell type corr matrix from crossSpeciesCellTypeGOCorr
#' @param corr_threshould minimum corr value for positively related cell types, default 0.6
#' @param ... additional params for sankeyNetwork
#' @return a Sankey plot showing related cell types
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' go_seurat_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#'
#' cell_type_go = getCellTypeGO(go_seurat_obj = go_seurat_obj, cell_type_co = "cell_type_annotation")
#' corr_matrix = cellTypeGOCorr(cell_type_go = cell_type_go, corr_method = "pearson")
#'
#' plotCellTypeSankey(corr_matrix = corr_matrix, 0.1)
#' }
#' @importFrom networkD3 sankeyNetwork
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
#'


plotCellTypeSankey <- function(corr_matrix, corr_threshould = 0.1, ...) {
  links <- corr_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "source") %>%
    tidyr::gather(key = "target", value = "value", -1) %>%
    dplyr::filter(value >= corr_threshould)

  if (nrow(links) == 0) {
    stop("no cell type pairs pass corr_threshould, please lower threshould")
  }

  nodes <- data.frame(
    name = c(as.character(links$source), as.character(links$target)) %>%
      unique()
  )

  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1

  p <- networkD3::sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "IDsource", Target = "IDtarget",
    Value = "value", NodeID = "name",
    sinksRight = FALSE, ...
  )

  return(p)
}


#' query co-up and co-down regulated GO terms from certain cell type pairs
#' @name getCellTypeSharedTerms
#' @param shared_go cell type shared GO table from getCellTypeSharedGO
#' @param cell_type_sp1 cell type from sp1 to query
#' @param cell_type_sp2 cell type from sp2 to query
#' @param return_full if return also pvals and logfc info, default FALSE
#' @param arrange_avg_log2FC arrange result by decreasing mean avg_log2FC, default TRUE
#' @return a dataframe displaying co-up or co-down regulated GO terms for the queried cell type pair
#' @examples
#' \donttest{
#' library(scGOclust)
#' library(httr)
#' httr::set_config(httr::config(ssl_verifypeer = FALSE))
#' data(mmu_tbl)
#' data(mmu_subset)
#' data(dme_tbl)
#' data(dme_subset)
#'
#' mmu_go_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
#'  seurat_obj = mmu_subset,
#'  feature_type = "external_gene_name")
#' dme_go_obj = makeGOSeurat(ensembl_to_GO = dme_tbl,
#'  seurat_obj = dme_subset,
#'  feature_type = "external_gene_name")
#'
#'
#' mmu_go_obj_analyzed = analyzeGOSeurat(mmu_go_obj, "cell_type_annotation")
#' dme_go_obj_analyzed = analyzeGOSeurat(dme_go_obj, "annotation")
#'
#' shared_go = getCellTypeSharedGO(species_1 = 'mmusculus',
#' species_2 = 'dmelanogaster',
#' analyzed_go_seurat_sp1 = mmu_go_obj_analyzed,
#' analyzed_go_seurat_sp2 = dme_go_obj_analyzed,
#' cell_type_col_sp1 = 'cell_type_annotation',
#' cell_type_col_sp2 = 'annotation',
#' layer_use = "data",
#' p_val_threshould = 0.01)
#'
#'
#' getCellTypeSharedTerms(shared_go = shared_go,
#' cell_type_sp1 = 'intestine_Enteroendocrine cell',
#' cell_type_sp2 = 'enteroendocrine cell',
#' return_full = FALSE)
#' }
#' @importFrom dplyr select filter desc
#' @importFrom magrittr %>%
#' @export
#'

getCellTypeSharedTerms <- function(shared_go, cell_type_sp1, cell_type_sp2, return_full = FALSE, arrange_avg_log2FC = TRUE) {

  ## check shared_go format
  if (!("p_val_sp1" %in% colnames(shared_go$shared_sig_markers))) {
    stop("Shared go table format incorrect, please re-compute using getCellTypeSharedGO")
  }

  if (!(cell_type_sp1 %in% levels(factor(shared_go$shared_sig_markers[["cluster_sp1"]])))) {
    message("available cell types from species 1: ")
    message(levels(factor(shared_go$shared_sig_markers[["sp1_cluster"]])))
    stop("cell type sp1 not in data, please check input")
  }

  if (!(cell_type_sp2 %in% levels(factor(shared_go$shared_sig_markers[["cluster_sp2"]])))) {
    message("available cell types from species 2: ")
    message(levels(factor(shared_go$shared_sig_markers[["sp2_cluster"]])))
    stop("cell type sp2 not in data, please check input")
  }

  if (!(return_full)) {
    tbl_ct <- shared_go$shared_sig_markers %>%
      dplyr::select(
        gene, cluster_sp1, marker_type_sp1,
        cluster_sp2, marker_type_sp2,
        avg_log2FC_sp1, avg_log2FC_sp2
      ) %>%
      dplyr::filter(cluster_sp1 == cell_type_sp1) %>%
      dplyr::filter(cluster_sp2 == cell_type_sp2)

    if(arrange_avg_log2FC){
      message("return shared terms arranged in decreasing mean avg_log2FC between species")
      tbl_ct <- tbl_ct %>% filter(avg_log2FC_sp1 > 0) %>% filter(avg_log2FC_sp2 > 0) %>%
        mutate(mean_avg_log2FC = ((avg_log2FC_sp1 + avg_log2FC_sp2) / 2)) %>%
        arrange(desc(mean_avg_log2FC))

    }

    return(tbl_ct)

  } else {

    tbl_ct <- shared_go$shared_sig_markers %>%
      dplyr::filter(cluster_sp1 == cell_type_sp1) %>%
      dplyr::filter(cluster_sp2 == cell_type_sp2)

    if(arrange_avg_log2FC){
      message("return shared terms arranged in decreasing mean avg_log2FC between species")
      tbl_ct <- tbl_ct %>% filter(avg_log2FC_sp1 > 0) %>% filter(avg_log2FC_sp2 > 0) %>%
        mutate(mean_avg_log2FC = ((avg_log2FC_sp1 + avg_log2FC_sp2) / 2)) %>%
        arrange(desc(mean_avg_log2FC))

    }

    return(tbl_ct)
  }
}
