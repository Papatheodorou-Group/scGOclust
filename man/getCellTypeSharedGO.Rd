% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{getCellTypeSharedGO}
\alias{getCellTypeSharedGO}
\title{get shared up and down regulated GO terms for all pairs of cell types}
\usage{
getCellTypeSharedGO(
  species_1,
  species_2,
  analyzed_go_seurat_sp1,
  analyzed_go_seurat_sp2,
  cell_type_col_sp1,
  cell_type_col_sp2,
  layer_use = "data",
  p_val_threshould = 0.01
)
}
\arguments{
\item{species_1}{name of species one}

\item{species_2}{name of species two}

\item{analyzed_go_seurat_sp1}{analyzed GO seurat object of species one}

\item{analyzed_go_seurat_sp2}{analyzed GO seurat object of species two}

\item{cell_type_col_sp1}{cell type column name for species 1 data}

\item{cell_type_col_sp2}{cell type column name for species 2 data}

\item{layer_use}{layer to use for marker computation, default 'data' which after NormalizeData will be log1p normalized data.}

\item{p_val_threshould}{p value threshold for selecting DEG (p_adjust)}
}
\value{
a list with sp1 raw, sp2 raw and shared, significant up and down regulated GO terms per cell type (pair)
}
\description{
get shared up and down regulated GO terms for all pairs of cell types
}
\examples{
\donttest{
library(scGOclust)
library(httr)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
data(mmu_tbl)
data(mmu_subset)
data(dme_tbl)
data(dme_subset)

mmu_go_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
 seurat_obj = mmu_subset,
 feature_type = "external_gene_name")
dme_go_obj = makeGOSeurat(ensembl_to_GO = dme_tbl,
 seurat_obj = dme_subset,
 feature_type = "external_gene_name")


mmu_go_obj_analyzed = analyzeGOSeurat(mmu_go_obj, "cell_type_annotation")
dme_go_obj_analyzed = analyzeGOSeurat(dme_go_obj, "annotation")

getCellTypeSharedGO(species_1 = 'mmusculus',
species_2 = 'dmelanogaster',
analyzed_go_seurat_sp1 =  mmu_go_obj_analyzed,
analyzed_go_seurat_sp2 =  dme_go_obj_analyzed,
cell_type_col_sp1 = 'cell_type_annotation',
cell_type_col_sp2 = 'annotation',
layer_use = "data",
p_val_threshould = 0.01)
}
}
