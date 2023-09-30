# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song ysong@ebi.ac.uk

#' Mouse stomach and intestine scRNA-seq data, microwell-seq
#' Subset to 50 cells per cell type as an example data
#' @name mmu_subset
#' @format a `Seurat` object
#' @source <https://bis.zju.edu.cn/MCA/>
"mmu_subset"



#' Drosophila gut scRNA-seq data, 10X Chromium
#' Subset to 45 cells per cell type as an example data
#' @name dme_subset
#' @format a `Seurat` object
#' @source <https://flycellatlas.org/>
"dme_subset"


#' Mouse EMSEMBL gene and GO annotation, subset to genes present in `mmu_subset`
#' @name mmu_tbl
#' @format a `data.frame` object
#' @source <http://www.ensembl.org/>
"mmu_tbl"


#' Drosophila EMSEMBL gene and GO annotation, subset to genes present in `dme_subset`
#' @name dme_tbl
#' @format a `data.frame` object
#' @source <http://www.ensembl.org/>
"dme_tbl"
