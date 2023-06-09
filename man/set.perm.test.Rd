% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_enrichment_statistics.R
\name{set.perm.test}
\alias{set.perm.test}
\title{Permutation test. normalized scores and Fisher test from leading-edge subset}
\usage{
set.perm.test(ranked.list, gene.set)
}
\arguments{
\item{ranked.list}{A ranked list of objects}

\item{gene.set}{A set of objects representing a process}
}
\value{
A list containing:
   ES.pvalue = ES pvalue of the set.
   NES = normalized score for the set.
   odd.ratio = odd.ratio from the fisher test.
   fp.value = pvalue of fisher test.
}
\description{
Permutation test. normalized scores and Fisher test from leading-edge subset
}
\examples{
#ranked.list <- readRDS("data/ranked_list.rds)
#gene.set <- readRDS("data/one_gene_set.rds")
#set.perm.test(ranked.list, gene.set)


}
