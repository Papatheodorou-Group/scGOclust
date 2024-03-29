% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{getCellTypeGO}
\alias{getCellTypeGO}
\title{get per cell type average scaled vector of GO terms}
\usage{
getCellTypeGO(go_seurat_obj, cell_type_col, norm_log1p = TRUE)
}
\arguments{
\item{go_seurat_obj}{go seurat object created by makeGOSeurat}

\item{cell_type_col}{column name in mera.data storing cell type classes}

\item{norm_log1p}{whether or not to perform data normalisation and log1p transformation, default TRUE}
}
\value{
a table of scaled GO representation per cell type (averaged)
}
\description{
get per cell type average scaled vector of GO terms
}
\examples{
\donttest{
library(scGOclust)
library(httr)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
data(mmu_tbl)
data(mmu_subset)
go_seurat_obj = makeGOSeurat(ensembl_to_GO = mmu_tbl,
 seurat_obj = mmu_subset,
 feature_type = "external_gene_name")
getCellTypeGO(go_seurat_obj = go_seurat_obj, cell_type_co = "cell_type_annotation")
}
}
