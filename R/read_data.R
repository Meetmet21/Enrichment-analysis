#' Read one set from .gmt file
#'
#' @description
#' Format one set into a data frame : one set correspond to one line in .gmt format.
#'
#' @param set One line from .gmt file
#' @return A data frame with the set name and the associated gene symbol. Columns: gs_name, gene_symbol.
#' @examples
#' lines <- readLines("data/ReactomePathways.gmt")
#' read.oneset.gmt(lines[1])
#'
read.oneset.gmt <- function(set){
  vec <- unlist(strsplit(set,"\t"))
  df <- data.frame(gs_name = gsub(pattern = "HALLMARK_",replacement = "",x = vec[1]),
                   gene_symbol = vec[3:length(vec)])
  return(df)
}

#' Read whole .gmt file
#'
#' @description
#' Read all lines through .gmt file and compute a data frame for each set to finally bind them in one total data frame.
#'
#' @param path.file Path for the .gmt file
#' @return A data frame containing sets with their corresponding gene symbols. COlumns: gs_name, gene_symbol
#' @examples
#' data(gmt_file)
#' read_set_gmt(gmt_file)
#' @export
#'
read_set_gmt <- function(path.file){
  vec.sets <- as.list(readLines(path.file))
  df <- Reduce(rbind,lapply(X = vec.sets,FUN = read.oneset.gmt))
  return(df)
}

#' Get sets of genes for Homo sapiens
#'
#' @description
#' This function uses different packages to access API of databases in order to retrieve gene sets.
#' @details
#' Reactome is accessed from ReactomePA package. To map gene ids to symbols, the mapIds for AnnotationDbi is used and the mapping is made on the Org.Hs.eg.db databse.
#' GO are accessed from the org.Hs.eg.db with the AnnotationDbi package keys function. Only GO ids are retrieved so they have to be mapped to their corresponding gene symbols with mapIds from AnnotationDbi, also done for GO terms.
#' KEGG is accessed with the KEGGREST package for the pathways. Than gene ids are mapped to their corresponding symbols with mapIds as before.
#' All the result are formated into data frames with two columns nammed as "gs_name" and "gene_symbol", All the resulting dataframes are stored in a list.
#'
#' @return
#' A list containing different data frames with gene set symbols and associated genes symbols from the three databases. COlumns: gs_name, gene_symbol.
#'    \item{-}{[1] Reactome data frame}
#'    \item{-}{[2] list of data frames containing different aspects for GO terms}
#'        \item{-}{  [2][1] Cellulare component}
#'        \item{-}{  [2][2] Biological process}
#'        \item{-}{  [2][3] Molecular function}
#'    \item{-}{[3] KEGG data frane}
#' @examples
#' res <- get_web_sets()
#' str(res)
#' @export
#'
get_web_sets <- function(){
  ## REACTOME **
  raw.data.reactome <- ReactomePA::gson_Reactome(organism = "human")
  data.reactome <- as.data.frame(raw.data.reactome@gsid2gene)
  colnames(data.reactome) <- c("gs_name","gene_symbol")
  data.reactome$gene_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=data.reactome$gene_symbol,
                                          keytype="ENTREZID", column="SYMBOL")
  data.reactome$gs_name <- unlist(as.list(reactome.db::reactomePATHID2NAME)[data.reactome$gs_name])

  ## Gene ontology
  go_ids <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "GO")
  gene_sy <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, go_ids,
                               c("GO","SYMBOL"),"GO")
  go_name <- AnnotationDbi::select(GO.db::GO.db, gene_sy$GO,
                                   c("GOID","TERM"),"GOID")
  data.go <- data.frame(gs_name = go_name[,2],
                   gene_symbol = gene_sy$SYMBOL,
                   aspect = gene_sy$ONTOLOGY)

  ## KEGG
  path <- KEGGREST::keggLink("pathway","hsa")
  path.gene.id <- sub("hsa:","", names(path))
  path.name.id <-  gsub("path:","",as.vector(path))
  match.name <- gsub(" - Homo sapiens (human)","",
                     KEGGREST::keggList("pathway","hsa"),
                     fixed = TRUE)

  path.gene.symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                            keys=path.gene.id,
                                            keytype="ENTREZID", column="SYMBOL")
  path.name <- unlist(lapply(path.name.id, function(id){
    match.name[[id]]
  }))
  data.kegg <- data.frame(gs_name = path.name,
                          gene_symbol = path.gene.symbol)

  #Final result
  return(list(Reactome = data.reactome,
              GO = list(CC = data.go[data.go$aspect=="CC",c(1,2)],
                        BP = data.go[data.go$aspect=="BP",c(1,2)],
                        MF = data.go[data.go$aspect=="MF",c(1,2)]),
              KEGG = data.kegg))
}

#' Read .gaf file
#'
#' @description
#' This function reads .gaf file and format it to be used by analysis functions for a given GO aspect.
#'
#' @details
#' This function also maps GO ids to their corresponding terms with the AnnotationDbi package via mapIds function.
#'
#' @param gaf a gaf file
#' @param aspect one of the GO aspect :"C" : Cellular component "P" : Biological process "M" : Molecular function
#' @return a data frame containing sets with their corresponding gene symbols. Columns: gs_name, gene_symbol.
#' @examples
#' data(gaf_file)
#' read_set_gaf(gaf_file, "C")
#' read_set_gaf(gaf_file, "M")
#' read_set_gaf(gaf_file, "P")
#' @export
read_set_gaf <- function(gaf, aspect){
  #Reads lines from .gaf file
  lines <- readLines(gaf)
  #Selects lines containing informaiton
  lines <- lines[which(!startsWith(lines,"!"))]
  #Splits each lines to obtain each columns for a given gene
  vecs <- strsplit(lines, "\t", fixed = TRUE)
  #Selects lines corresponding to the aspect
  res <- lapply(vecs, function(vec){
    if(vec[9] == toupper(aspect)){
      list(vec[5],vec[3])
    }
  })
  #Remove NULL values
  res <- Filter(Negate(is.null), res)
  #Retrieves gene set names as GO ID
  gs_name <- unlist(lapply(res,function(x){x[1]}))
  #Retrieves gene symbol from each line
  gene_symbol <- unlist(lapply(res,function(x){x[2]}))
  #Maps GO IDS to their corresponding terms
  terms <- AnnotationDbi::select(GO.db::GO.db,
                                 keys = gs_name,
                                 columns = c("GOID","TERM"),
                                 keytype = "GOID")
  #Outputs result in a data frame with non repetitive values
  unique(data.frame(gs_name = terms[,2], gene_symbol = gene_symbol))
}

