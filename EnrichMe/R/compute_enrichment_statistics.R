#' Maximum deviation from 0
#'
#' @description
#' This function compute the maximum deviation from 0 of a vector of values.
#' @details
#' To compute the maximum deviation from 0, this function first search the maximum and the minimum in the numerical vector (argument : table).
#' Then, it compares the max and absolute value of the min and returns the bigger one.
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
#' @description
#' The enrichment score for a given gene set is calculated using a Kolmogorov-Smirnov (KS) test-like algorithm.
#'
#' @details
#' This function fixes the weight for the running sum calculation as : sqrt((m.list - n.set)/n.set) for genes in the set and sqrt(n.set/(n.list - n.set)) in the other case.
#' Then, it goes through the ranked list of genes and sums the wieghts depedning if the gene is in or out the set.
#' Finally, it looks through the results of the cumulative sum and returns the maximum value deviating from 0 and it's index in the ranked list of genes.
#'
#' @param ranked.list ranked list of objects.
#' @param gene.set set of objects representing a process.
#' @returns A list containing the ES of the set and it's position in the ranked.list.
#' \item{ES}{The enrichment score for a given set.}
#' \item{ES.id}{The index of the maximum KS-cumulative sum for a given set (corresepond to the Enrichment score).}
#' @examples
#' #data(ranked.list)
#' #data(fatty_acid_metabolism)
#' #compute.ES(ranked.list, fatty_acid_metabolism)
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

#' Permutation test. normalized scores and Fisher test from leading-edge subset for one set
#'
#' @description
#' This function computes from the ES of a given set, some statistics to assess the significance of the set's ES.
#' @details
#' To this purpose, a permutation test is performed to build an empirical normal distribution of ESs by shuffling the ranked list of genes. Then, the pvalue is computed by comparing
#' this distribution to the set ES. Also, the normalized ES for the set can be assess from the normal distribution similar to Z-score.
#' The second part of this function uses the leading-edge subset genes to run a Fisher's exact test. The leading-edge subset is fixed as the top genes participating to the ES of the set and are retrieved by
#' taking all the genes before reaching the maximum ES. The two factor compared are the presence of genes from the subset in the set and the presence of genes in the set not from the subset.
#'
#' @param ranked.list A ranked list of objects.
#' @param gene.set A set of objects representing a process.
#' @param n.perm The number of permutation test. By default fixed to 1000.
#' @returns A list containing:
#'    \item{ES.pvalue}{ES pvalue of the set.}
#'    \item{NES}{Normalized score for the set.}
#'    \item{odd.ratio}{Odd.ratio from the fisher test.}
#'    \item{fp.value}{P-value of fisher test.}
#' @examples
#' #data(ranked.list)
#' #data(fatty_acid_metabolism)
#' #one_set_analysis(ranked.list,fatty_acid_metabolism)
#' @export
#'
#'
one_set_analysis <- function(ranked.list, gene.set, n.perm = 1000){
  ######### Permutation test ##########

  #Observed ES for a given set
  obs <- compute.ES(ranked.list, gene.set)
  # number of permutation to build then null distribution of ES
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
#' @description
#' This function is similar to one_set_analysis but generalized to many sets. It also correct the p-values for multiple testing.
#'
#' @details
#' The analysis run for each set is the same as the one_set_analysis function, but as many statistical tests are run for all the sets, this function
#' correct the p-values by using the build-in function p.adjust with the default method.
#'
#' @param ranked.list A ranked list of objects
#' @param pathways a data frame containing all sets names as gs_name in one column and the associated objects in another columns as gene_symbol
#' @param n.perm The number of permutation test. By default fixed to 1000.
#' @return
#' A data frame with the following columns for each set of object :
#'    \item{set}{The set name.}
#'    \item{NES}{The normalized ES for the set.}
#'    \item{p.value}{The permutation test nominal p-value for the ES.}
#'    \item{odd.ratio}{The odd.ratio from the Fisher's test.}
#'    \item{fp.value}{The p-value of the Fisher's test.}
#' @export
#' @examples
#' #data(ranked.list)
#' #data(pathways)
#' #multiple_set_analysis(ranked.list, pathways)
multiple_set_analysis <- function(ranked.list, pathways, n.perm = 1000){
  res.sets <- lapply(unique(pathways$gs_name),function(set){
    cat("Working on :", set,"\n")
    gene.set <- pathways$gene_symbol[pathways$gs_name == set]
    res.perm <- one_set_analysis(ranked.list, gene.set, n.perm)
    df <- data.frame(set = set, NES = res.perm$NES, p.value = res.perm$ES.pvalue,
                     odd.ratio = res.perm$odd.ratio, fp.value = res.perm$fp.value)
    return(df)
  })

  final <- Reduce(rbind, res.sets)
  final$fp.value <- p.adjust(final$fp.value)
  final$p.value <- p.adjust(final$p.value)
  return(final)
}
