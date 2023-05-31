## ----install and load packages------------------------------------------------
# load required libraries

library(Seurat)
library(pheatmap)
library(httr)

## if (!require("devtools")) install.packages("devtools")

## install latest from source
## for reprodubcibility we do not update dependencies
# devtools::install_github("YY-SONG0718/scGOclust", upgrade_dependencies = FALSE)

library(scGOclust)

#

## ----load_input---------------------------------------------------------------
# get a gene to GO BP terms mapping table
# remove electronically inferred records

# sometimes ensembl complains about ssh certificate has expired, this is a known issue, run this code
httr::set_config(httr::config(ssl_verifypeer = FALSE)) 

mmu_tbl = ensemblToGo(species = 'mmusculus', GO_linkage_type = c('experimental', 'phylogenetic', 'computational', 'author', 'curator' ))
dme_tbl = ensemblToGo(species = 'dmelanogaster', GO_linkage_type = c('experimental', 'phylogenetic', 'computational', 'author', 'curator' ))


## -----------------------------------------------------------------------------
# load the gene expression raw count objects
data(mmu_subset)
data(dme_subset)
ls()

## ----build_GO_BP_profile------------------------------------------------------
## construct a Seurat object with GO BP as features

mmu_go_obj <- makeGOSeurat(ensembl_to_GO = mmu_tbl, feature_type = 'external_gene_name', seurat_obj = mmu_subset)

dme_go_obj <- makeGOSeurat(ensembl_to_GO = dme_tbl, feature_type = 'external_gene_name', seurat_obj = dme_subset)


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
#  ## calculation takes a few minutes due to the Wilcoxon signed rank test
#  
#  ct_shared_go = getCellTypeSharedGO(species_1 = 'mmusculus', species_2 = 'dmelanogaster', analyzed_go_seurat_sp1 = mmu_go_analyzed, analyzed_go_seurat_sp2 = dme_go_analyzed, cell_type_col_sp1 = 'cell_type_annotation', cell_type_col_sp2 = 'annotation')
#  
#  head(ct_shared_go)

## ----shared go cell type, eval = FALSE----------------------------------------
#  
#  # query shared GO terms for specific cell type pairs
#  
#  getCellTypeSharedTerms(shared_go = ct_shared_go,
#                         cell_type_sp1 = 'intestine_Enteroendocrine cell',
#                         cell_type_sp2 = 'enteroendocrine cell',
#                         return_full = FALSE)
#  
#  
#  
#  

## -----------------------------------------------------------------------------

plotCellTypeSankey(corr_matrix = corr, corr_threshould = 0.05)

## ----sessioninfo--------------------------------------------------------------

sessionInfo()


