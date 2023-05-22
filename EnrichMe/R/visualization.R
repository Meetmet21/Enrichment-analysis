#' Visualization functions for Enrichment analysis results

#' Function to plot a dot plot chart from the odd ration, the Fisher's test p.value and the normalized enrichment score of signficantly enriched èathways.

#' @param data the resulting data frame from the multiple_set_analysis function
#' @return a dotplot chart summarizing signficant results
#' @examples
#' ranked.list <- data(ranked.list)
#' pathways <- data(pathways)
#' data <- multiple_set_analysis(ranked.list, pathways)
#' dotplot(data)
#'
dotplot <- function(data){
  data <- data[data$fp.value <= 0.05,]
  ggplot2::ggplot(data = data, ggplot2::aes(x = NES, y = set,
                          color = fp.value, size = odd.ratio)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "red", high = "blue") +
    ggplot2::theme_bw() +
    ggplot2::ylab("") +
    ggplot2::xlab("NES") +
    ggplot2::ggtitle("Enrichment analysis plot")
}

