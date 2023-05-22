#' Maximum deviation from 0
#'
#' This function compute the maximum deviation from 0 of a vector of values.
#'
#' @param table A vector containing numerical values.
#' @return The maximum deviation from 0 among the values in the table.
#' @examples
#' max.deviation(c(seq(1:10),seq(-13)))
#'
max.deviation <- function(table){
  max <- max(table)
  min <- min(table)
  if (abs(min) >= max) return(min)
  else{return(max)}
}


#' compute the ES for a set
#'
#' The enrichment score for a given gene set is calculated using a Kolmogorov-Smirnov (KS) test-like algorithm that walks down the ranked list of genes, increasing a running-sum statistic when a gene in the set is encountered and decreasing the statistic when a gene outside the set is encountered
#'
#' @param ranked.list ranked list of objects
#' @param gene.set set of objects representing a process
#' @returns A list containing the ES of the set and it's position from the ES table
#' @examples
#' #ranked.list <- readRDS("data/ranked_list.rds")
#' #gene.set <- readRDS("data/one_gene_set.rds")
#' #compute.ES(ranked.list, gene.set)
#'
#'
compute.ES <- function(ranked.list, gene.set){
  # length of each list
  n.set <- length(gene.set)
  n.list <- length(ranked.list)
  # weights for the ks cumulative sum
  weight.in <- sqrt((n.list - n.set)/n.set)
  weight.out <- sqrt(n.set/(n.list - n.set))
  #match matrix between ranked.list and gene.set
  matrix.in <- sign(match(ranked.list, gene.set, nomatch = 0))
  matrix.out <- 1 - matrix.in
  #cumulative sum to compute the ES
  ES.sums <- cumsum(matrix.in * weight.in - matrix.out * weight.out)
  #give the max deviation from 0 so ES for the set
  ES <- max.deviation(ES.sums)
  #Give the ined of the maximum ES from the ranked list
  ES.id <- which(ES.sums == ES)
  #returns ES (max deviation from 0) and the leading-edge gene subsets
  return(list(ES = ES, ID = ES.id))
}

#' Permutation test. normalized scores and Fisher test from leading-edge subset
#'
#' @param ranked.list A ranked list of objects
#' @param gene.set A set of objects representing a process
#' @returns A list containing:
#'    ES.pvalue = ES pvalue of the set.
#'    NES = normalized score for the set.
#'    odd.ratio = odd.ratio from the fisher test.
#'    fp.value = pvalue of fisher test.
#' @examples
#' #ranked.list <- readRDS("data/ranked_list.rds)
#' #gene.set <- readRDS("data/one_gene_set.rds")
#' #set.perm.test(ranked.list, gene.set)
#'
#'
set.perm.test <- function(ranked.list, gene.set){
  ######### Permutation test ##########

  #Observed ES for a given set
  obs <- compute.ES(ranked.list, gene.set)
  # number of permutation to build then null distribution of ES
  n.perm <- 1000
  # Randomly ordered list of gene sampled from the ranked list to calculate null ES
  # Distribution.
  random.list <- lapply(1:n.perm, function(x){sample(x = ranked.list,
                                                     size = length(ranked.list),
                                                     replace = F)})
  ## Compute ES for each randomly ordered list of genes
  ES.dist <- unlist(lapply(random.list,
                            function(x){compute.ES(x,
                                                    gene.set)$ES}))
  ## Calculate the normalized ES with the z-score from the null distribution of permuted ES.
  NES <- (obs$ES - mean(ES.dist))/sd(ES.dist)
  #3 Emperical pvalue for our ES
  if(obs$ES >= 0){
    ES.pvalue <- sum(ES.dist>=obs$ES)/n.perm
  }else{
    ES.pvalue <- sum(ES.dist<=obs$ES)/n.perm
  }

  ######## Fisher test #########
  if (obs$ES > 0){
    core.genes <- ranked.list[1:obs$ID]
    not.core.genes <- ranked.list[(obs$ID + 1):length(ranked.list)]
    #in LE sets and in the giver set
    A <- sum(core.genes %in% gene.set)
    B <- length(core.genes) - A
    C <- sum(not.core.genes %in% gene.set)
    D <- length(core.genes) - C
    #contingency table
    m <- rbind(c(A,B),c(C,D))
    #associated pvalue from fisher test and the odd.ratio
    res <- fisher.test(m)
  }
  else{
    core.genes <- ranked.list[obs$ID:length(ranked.list)]
    not.core.genes <- ranked.list[1:(obs$ID - 1)]
    #in LE sets and in the giver set
    A <- sum(core.genes %in% gene.set)
    B <- length(core.genes) - A
    C <- sum(not.core.genes %in% gene.set)
    D <- length(core.genes) - C
    #contingency table
    m <- rbind(c(A,B),c(C,D))
    #associated pvalue from fisher test and the odd.ratio
    res <- fisher.test(m)
  }
  return(list( ES.pvalue = ES.pvalue,NES = NES,
               odd.ratio = res$estimate[[1]], fp.value = res$p.value))
}

#' Run GSEA and Fisher test from leading-edge subset for multiple sets of object
#'
#' @param ranked.list A ranked list of objects
#' @param pathways a data frame containing all sets names as gs_name in one column and the associated objects in another columns as gene_symbol
#' @return
#' A data frame with the following columns for each set of object :
#'    set = the set name.
#'    NES = the normalized ES for the set.
#'    p.value = The permutation test nominal pvalue for the ES.
#'    odd.ratio = The odd.ratio from the Fisher's test.
#'    fp.value = the pvalue of the Fisher's test.
#' @export
#' @examples
#' #ranked.list <- readRDS("data/ranked_list.rds)
#' #pathways <- readRDS("data/pathways.rds")
#' #multiple_set_analysis(ranked.list, pathways)
multiple_set_analysis <- function(ranked.list, pathways){
  res.sets <- lapply(unique(pathways$gs_name),function(set){
    cat("Working on :", set,"\n")
    gene.set <- pathways$gene_symbol[pathways$gs_name == set]
    res.perm <- set.perm.test(ranked.list, gene.set)
    df <- data.frame(set = set, NES = res.perm$NES, p.value = res.perm$ES.pvalue,
                     odd.ratio = res.perm$odd.ratio, fp.value = res.perm$fp.value)
    return(df)
  })

  final <- Reduce(rbind, res.sets)
  final$fp.value <- p.adjust(final$fp.value)
  final$p.value <- p.adjust(final$p.value)
  final <- final[]
  return(final)
}
