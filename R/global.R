# © EMBL-European Bioinformatics Institute, 2023
# Yuyao Song ysong@ebi.ac.uk

#' record some global variables: pre-defined column name in biomaRt query and markers
#' @name globalvariables
#' @importFrom utils globalVariables

utils::globalVariables(c("namespace_1003", "name_1006", "go_linkage_type", "cluster", "avg_log2FC", "p_val_adj", "gene", 'placehold'))
utils::globalVariables(c("cluster_sp1", "gene", "marker_type_sp1", "avg_log2FC_sp1",
                         "cluster_sp2", "marker_type_sp2", "avg_log2FC_sp2", "mean_avg_log2FC"))
utils::globalVariables(c('value'))
