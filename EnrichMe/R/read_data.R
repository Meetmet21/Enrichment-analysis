#' Read one set from .gmt file
#'
#' Format one set into a data frame : one set correspond to one line in .gmt format.
#'
#' @param set One line from .gmt file
#' @return A data frame with the set name and the associated gene symbols
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
#' Read all lines through .gmt file and compute a data frame for each set to finally bind them in one total data frame.
#'
#' @param path.file Path for the .gmt file
#' @return A data frame containing sets with their corresponding gene symbols.
#' @examples
#' gmt <- data(gmt_file)
#' read_set_gmt(gmt)
#' @export
#'
read_set_gmt <- function(path.file){
  vec.sets <- as.list(readLines(path.file))
  df <- Reduce(rbind,lapply(X = vec.sets,FUN = read.oneset.gmt))
  return(df)
}

#' Get sets of genes for Homo sapiens from Reactome, Gene ontology and KEGG from the web
#'
#' This function uses different packages to access API of databases in order to retrieve gene sets.
#' Reactome is accessed from ReactomePA package.
#' GO terms are downloaded from the gene ontology database for different aspects.
#' KEGG is accessed with the KEGGREST package
#'
#' @return
#' A list containing different data frames with gene set symbols and associated genes symbols from the three databases.
#'    \item{[1] Reactome data frame}
#'    \item{[2] list of data frames containing different aspects for GO terms}
#'        \item{  [2][1] Cellulare component}
#'        \item{  [2][2] Biological process}
#'        \item{  [2][3] Molecular function}
#'    \item{[3] KEGG data frane}
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
#' This function reads .gaf file and format it to be used analysis functions for a given GO aspect.
#'
#' @param gaf a gaf file
#' @param aspect one of the GO aspect :"C" : Cellulare component "P" : Biological process "M" : Molecular function
#' @return a data frame containing sets with their corresponding gene symbols.
#' @examples
#' gaf <- data(gaf_file)
#' CC <- read_set_gaf(gaf, "C")
#' MF <- read_set_gaf(gaf, "M")
#' BP <- read_set_gaf(gaf, "P")
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

