#' Mouse gene symbols related to fatty acid metabolism pathway
#'
#' A vector of character containing 158 gene symbols related with the fatty acid metabolism pathway in Mouse.
#' @examples
#' data(fatty_acid_metabolism)
#'
#' @format A vector of character with 158 elements.
"fatty_acid_metabolism"

#' Pathways with their associated genes for Mouse
#'
#' A data frame with pathways in one column (gs_name) and their associated genes (gene_symbol).
#' @examples
#' data(pathways)
#'
#' @format A data frame with 7393 observations of 2 variables.
#' \describe{
#'   \item{gs_name}{50 pathways from Mouse.}
#'   \item{gene_symbol}{4398 gene symbols from Mouse}
#' }
"pathways"

#' Ranked list of genes from gene expression experiment.
#'
#' A vector of character ranked by logFC of each gene from starving and normal condition in Mouse experiment.
#' @examples
#' data(ranked.list)
#'
#' @format A vector of character with 24421 elements.
"ranked.list"

#' .gaf file containing GO for Homo sapiens
#'
#' A large dataset containing GO information for human in GMT format.
#' @examples
#' data(gaf_file)
#'
#' @format GAF file
#' @source http://current.geneontology.org/products/pages/downloads.html
"gaf_file"

#' .gmt file containing GO for Homo sapiens
#'
#' A large dataset containing GO information for human hallmark in GMT format.
#' @examples
#' data(gmt_file)
#'
#' @format GMT file
#' @source https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
"gmt_file"
