% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.R
\name{get_web_sets}
\alias{get_web_sets}
\title{Get sets of genes for Homo sapiens}
\usage{
get_web_sets()
}
\value{
A list containing different data frames with gene set symbols and associated genes symbols from the three databases. COlumns: gs_name, gene_symbol.
   \item{-}{[1] Reactome data frame}
   \item{-}{[2] list of data frames containing different aspects for GO terms}
       \item{-}{  [2][1] Cellulare component}
       \item{-}{  [2][2] Biological process}
       \item{-}{  [2][3] Molecular function}
   \item{-}{[3] KEGG data frane}
}
\description{
This function uses different packages to access API of databases in order to retrieve gene sets.
}
\details{
Reactome is accessed from ReactomePA package. To map gene ids to symbols, the mapIds for AnnotationDbi is used and the mapping is made on the Org.Hs.eg.db databse.
GO are accessed from the org.Hs.eg.db with the AnnotationDbi package keys function. Only GO ids are retrieved so they have to be mapped to their corresponding gene symbols with mapIds from AnnotationDbi, also done for GO terms.
KEGG is accessed with the KEGGREST package for the pathways. Than gene ids are mapped to their corresponding symbols with mapIds as before.
All the result are formated into data frames with two columns nammed as "gs_name" and "gene_symbol", All the resulting dataframes are stored in a list.
}
\examples{
res <- get_web_sets()
str(res)
}
