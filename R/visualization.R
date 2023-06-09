#' Visualization functions for Enrichment analysis results
#'
#' @description
#' Function to plot a dot plot chart from the odd ration, the Fisher's test p.value and the normalized enrichment score of significantly enriched pathways.
#' @details
#' Only p-values lower or equal to 0,05 are represented in the chart.
#'
#' @param data the resulting data frame from the multiple_set_analysis function.
#' @return a dotplot chart summarizing significant results.
#' @examples
#' data(ranked.list)
#' data(pathways)
#' data <- multiple_set_analysis(ranked.list, pathways)
#' dotplot(data)
#' @export
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

