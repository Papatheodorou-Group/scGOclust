#' get requested ensembl ID to GO mapping table
#' @name ensemblToGo
#' @param species species name matching ensembl biomaRt naming, such as hsapiens, mmusculus
#' @param GO_type GO term type, choose among 'biological_process', 'molecular_function', 'cellular_component', default 'biological_process'
#' @param GO_linkage_type GO annotation evidence codes set. Choose among 'experimental', 'phylogenetic', 'computational', 'author', 'curator', 'electronic', defaut all
#' @return a table with ensembl to GO terms mapping including requested linkage type
#' @examples
#' \dontrun{
#' ensemblToGo(species = 'mmusculus', GO_type = 'biological_process', GO_linkage_type='experimental')
#' }
#' @importFrom biomaRt useMart useDataset getBM
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
#'
#'
ensemblToGo <- function(species, GO_type = 'biological_process', GO_linkage_type = c('experimental', 'phylogenetic', 'computational', 'author', 'curator', 'electronic')) {

  ## GO source type code
  go_source = list(experimental = c("EXP", 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'),
                   phylogenetic = c("IBA", 'IBD', 'IKR', 'IRD'),
                   computational = c("ISS", 'ISO', 'ISA', 'ISM', 'IGC', 'RCA'),
                   author = c('TAS', 'NAS'),
                   curator = c('IC', 'ND'),
                   electronic = c("IEA"))


  bm <- useMart("ensembl")
  bm <- useDataset(paste0(species, "_gene_ensembl"), mart=bm)

  # Get ensembl gene ids and GO terms

  message("query biomart")


  EG2GO <- tryCatch({
    getBM(mart=bm, attributes=c('ensembl_gene_id','external_gene_name','go_id', 'name_1006', 'go_linkage_type', 'namespace_1003'))
  }, warning = function(w) {
    message("ensembl biomaRt warning:")
    print(w)
  }, error = function(e) {
    message("ensembl biomaRt error:")
    printe
  }, finally={
    message("query biomaRt successful")
  })



  # Remove blank entries
  EG2GO <- EG2GO[EG2GO$go_id != '',]

  GO_terms = EG2GO %>% dplyr::filter(namespace_1003 == GO_type) %>%
                     dplyr::filter(!(name_1006 %in% c("biological_process", 'molecular_function', 'cellular_component')))

  types_in_source = list()

  for(source_type in names(go_source)){

    type_in = levels(factor(GO_terms$go_linkage_type))[levels(factor(GO_terms$go_linkage_type)) %in% go_source[[source_type]]]
    types_in_source[[source_type]] = type_in

  }

  included_terms = sapply(GO_linkage_type, function(x) go_source[[x]])
  print("including GO link types: ")
  print(included_terms)

  use = as.character(unname(unlist(sapply(GO_linkage_type, function(x) go_source[[x]]))))

  GO_terms = GO_terms %>% dplyr::filter(go_linkage_type %in% use)

  return(GO_terms)

}




#' create a seurat object with GO terms
#' @name makeGOSeurat
#' @param ensembl_to_GO ensembl_to_go mapping table from function ensemblToGo
#' @param seurat_obj count matrix with genes to cells
#' @param feature_type feature type of count matrix, choose from ensembl_gene_id, external_gene_name, default ensembl_gene_id
#' @return a seurast object with GO terms as features
#' @examples
#' \dontrun{
#' makeGOSeurat(ensembl_to_GO, seurat_obj, feature_type = 'ensembl_gene_id')
#' }
#' @importFrom biomaRt useMart useDataset getBM
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr filter mutate arrange select
#' @export
#'
#'


makeGOSeurat <- function(ensembl_to_GO, seurat_obj, feature_type = 'ensembl_gene_id'){

  counts = as.matrix(seurat_obj@assays$RNA@counts)

  ## pivot GO to feature type table to matrix
  go_matrix = tbl %>% dplyr::mutate(placehold = 1) %>%
    tidyr::pivot_wider(id_cols = ensembl_gene_id, names_from = name_1006, values_from = placehold, values_fill = 0, values_fn = 'first')

  shared = intersect(go_matrix[[feature_type]], rownames(counts))

  go_matrix <- go_matrix %>% dplyr::filter(ensembl_gene_id %in% shared) %>% dplyr::arrange(ensembl_gene_id)
  go_matrix = go_matrix %>% tibble::column_to_rownames(feature_type)

  counts = t(counts)
  counts = counts %>% as.data.frame()
  counts = counts  %>% dplyr::select(eval(shared))
  counts = counts[, match(rownames(go_matrix), colnames(counts))]


  if(!(all(colnames(counts) == rownames(go_matrix)))) {

    stop("error during calculation, please check input format")
  }

  message("compute GO to cell matrix")
  go_mtx = as.matrix(counts) %*% as.matrix(go_matrix)
  go_obj <- CreateSeuratObject(counts = t(go_mtx), meta.data = obj@meta.data)

  message("returning GO Seurat object")
  return(go_obj)


}

#' standard seurat analysis on GO_seurat object
#' @name analyzeGOSeurat
#' @param go_seurat_obj go seurat object created by makeGOSeurat
#' @param cell_type_col column name in mera.data storing cell type classes
#' @return standard analyzed GO seurat object until findAllMarkers
#' @examples
#' \dontrun{
#'
#' }
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters RunUMAP FindAllMarkers
#' @export
#'


analyzeGOSeurat <- function(go_seurat_obj){

  go_seurat_obj <- NormalizeData(go_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  go_seurat_obj <- FindVariableFeatures(go_seurat_obj, selection.method = "vst", nfeatures = 2000)
  go_seurat_obj <- ScaleData(object = go_seurat_obj, verbose = FALSE)
  go_seurat_obj <- RunPCA(object = go_seurat_obj, npcs = 50, verbose = FALSE)
  go_seurat_obj <- FindNeighbors(go_seurat_obj, dims = 1:50)
  go_seurat_obj <- FindClusters(go_seurat_obj, resolution = 1)
  go_seurat_obj <- RunUMAP(object = go_seurat_obj, reduction = "pca", min.dist = 0.3,
                    dims = 1:50)
  go_seurat_obj <- FindAllMarkers(go_seurat_obj, test.use = "wilcox")


}

#' get per cell type average scaled vector of GO terms
#' @name getCellTypeGO
#' @param go_seurat_obj go seurat object created by makeGOSeurat
#' @param cell_type_col column name in mera.data storing cell type classes
#' @return a table of scaled GO representation per cell type (averaged)
#' @examples
#' \dontrun{
#'
#' }
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData AverageExpression
#' @export
#'

getCellTypeGO <- function(go_seurat_obj, cell_type_col){

  go_seurat_obj <- NormalizeData(go_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

  go_seurat_obj <- FindVariableFeatures(go_seurat_obj, selection.method = "vst", nfeatures = 2000)

  go_seurat_obj <- ScaleData(object = go_seurat_obj, verbose = FALSE)

  tbl <- AverageExpression(object = go_seurat_obj, group.by = cell_type_col, slot = 'scale.data')
  return(tbl[['RNA']])

}

