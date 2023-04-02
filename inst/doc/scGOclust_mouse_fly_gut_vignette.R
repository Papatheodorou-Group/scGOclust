## ----install------------------------------------------------------------------
# load required libraries

library(Seurat)
library(SeuratDisk)
library(pheatmap)

if (!require("devtools")) install.packages("devtools")

#devtools::install_github("YY-SONG0718/scGOclust")

#library(scGOclust)

devtools::load_all("../")

## ----load_input---------------------------------------------------------------
# get a gene to GO BP terms mapping table
# remove electronically inferred records

mmu_tbl = ensemblToGo(species = 'mmusculus', GO_linkage_type = c('experimental', 'phylogenetic', 'computational', 'author', 'curator' ))
dme_tbl = ensemblToGo(species = 'dmelanogaster', GO_linkage_type = c('experimental', 'phylogenetic', 'computational', 'author', 'curator' ))

# load the gene expression raw count objects
mmu_obj <- LoadH5Seurat('/Users/ysong/SOFTWARE/scGOclust_data/mca_stomach_intestine_concat_counts.h5seurat')
dme_obj <- LoadH5Seurat('/Users/ysong/SOFTWARE/scGOclust_data/fca_10x_gut_cleaned_cell_types.h5seurat')



## ----build_GO_BP_profile------------------------------------------------------
## construct a Seurat object with GO BP as features

mmu_go_obj <- makeGOSeurat(ensembl_to_GO = mmu_tbl, feature_type = 'external_gene_name', seurat_obj = mmu_obj)

dme_go_obj <- makeGOSeurat(ensembl_to_GO = dme_tbl, feature_type = 'external_gene_name', seurat_obj = dme_obj)


## ----cell_type_BP-------------------------------------------------------------
# specify the column with cell type annotation in seurat_obj@meta.data

mmu_ct_go <- getCellTypeGO(go_seurat_obj = mmu_go_obj, cell_type_col = 'cell_type_annotation')
dme_ct_go <- getCellTypeGO(go_seurat_obj = dme_go_obj, cell_type_col = 'annotation')


## ----cell_type_corr_within, fig.width = 8,  fig.height = 8--------------------
# heatmap of Pearson's correlation coefficient of cell type average BP profiles within species

mmu_corr = cellTypeGOCorr(cell_type_go = mmu_ct_go, corr_method = 'pearson')
pheatmap(mmu_corr)

## ---- fig.height = 8, fig.width = 8-------------------------------------------
dme_corr = cellTypeGOCorr(cell_type_go = dme_ct_go, corr_method = 'pearson')
pheatmap(dme_corr)

## ----cell_type_corr_cross-----------------------------------------------------

# calculate Pearson's correlation coefficient of cell type average BP profiles across species

corr = crossSpeciesCellTypeGOCorr(species_1 = 'mmusculus', species_2 = 'dmelanogaster', cell_type_go_sp1 = mmu_ct_go, cell_type_go_sp2 = dme_ct_go, corr_method = 'pearson')


## ----heatmap_corr_cross, fig.width = 9,  fig.height = 10----------------------

# cross-species cell type profile heatmap

pheatmap(corr, width = 9, height = 10)

pheatmap(corr, scale = 'column', width = 9, height = 10)



## ----sheatmap_corr_cross,  fig.width = 9,  fig.height = 10--------------------

# sheatmap tries to put cells with higher values on the diagonal
# helpful when cross-species cell type similarity signal is less clear

slanter::sheatmap((corr + 0.5), width = 9, height = 10)

# scale by row or column to see relative similarity

slanter::sheatmap((corr + 0.5), scale = 'column', width = 9, height = 10)


## ----analyze_GO---------------------------------------------------------------

# analyze the cell-by-GO BP profile as a count matrix
mmu_go_analyzed = analyzeGOSeurat(go_seurat_obj = mmu_go_obj, cell_type_col = 'cell_type_annotation')

## ---- fig.width = 8, fig.height = 8-------------------------------------------
# UMAP plot of the analyzed cell-by-GO BP profile
# labeled by previously specified cell annotation column in meta.data

DimPlot(mmu_go_analyzed, label = TRUE) + NoLegend()

## -----------------------------------------------------------------------------
dme_go_analyzed = analyzeGOSeurat(go_seurat_obj = dme_go_obj, cell_type_col = 'annotation')

## ---- fig.width = 8, fig.height = 8-------------------------------------------
DimPlot(dme_go_analyzed, label = TRUE) + NoLegend()

## ----shared_go, eval = FALSE--------------------------------------------------
#  
#  ## calculation takes a few minuites due to the Wilcoxon signed rank test
#  
#  ct_shared_go = getCellTypeSharedGO(species_1 = 'mmusculus', species_2 = 'dmelanogaster', analyzed_go_seurat_sp1 = mmu_go_analyzed, analyzed_go_seurat_sp2 = dme_go_analyzed, cell_type_col_sp1 = 'cell_type_annotation', cell_type_col_sp2 = 'annotation')

## ----sessioninfo--------------------------------------------------------------

sessionInfo()


