############################
# Function get_annotations #
############################

#' This function annotates a column of transcripts or gene IDs (ENSEMBL) with information of the Biomart.
#' If transcript IDs are provided, they are also annotated with information of the genes to which they belong.
#' 
#' The Gene information added include:
#' - ENSEMBL ID
#' - HGNC Symbol and description
#' - Start and End position
#' - Gene length
#'
#' @param ensembl_ids The column of transcripts to be used as input.
#' @param mode To specify the IDs provided, between "transcripts" or "genes". Default = genes.
#' @param version This function can use the version 103 or the current version of the Biomart. Default = Current.
#' @param format The output will be a table named tx2gene in .csv or .xlsx formats. Default = csv.
#' @export

get_annotations <- function(ensembl_ids, mode = "genes", version = "", format = "csv") {

  require("biomaRt")
  require("ensembldb")
  
  if(version == "103"){
    ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl",
                      host = "feb2021.archive.ensembl.org")
  } else {
    ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl",
                      host = "useast.ensembl.org")
  }
  
  if(mode == "transcripts"){
    
    df <- data.frame(transcriptID = ensembl_ids)
    filename <- "tx2gene"
  
    # There are more annotations available in the biomaRt, check listAttributes(mart = ensembl)
    genemap <- getBM(attributes = c("ensembl_transcript_id_version",
                                    "ensembl_gene_id",
                                    "hgnc_symbol",
                                    "start_position", "end_position",
                                    "description"),
                     filters = "ensembl_transcript_id_version",
                     values = df$transcriptID,
                     mart = ensembl)
  
    idx <- match(df$transcriptID, genemap$ensembl_transcript_id_version)
  
    df$geneID <- genemap$ensembl_gene_id[idx]
    df$symbol <- genemap$hgnc_symbol[idx]
    df$gene_start <- genemap$start_position[idx]
    df$gene_end <- genemap$end_position[idx]
    df$description <- genemap$description[idx]
    df$gene_length <- df$gene_end - tx2gene$gene_start + 1
  
  } else {
    
    df <- data.frame(geneID = ensembl_ids)
    filename <- "gene_annotations"
    
    # There are more annotations available in the biomaRt, check listAttributes(mart = ensembl)
    genemap <- getBM(attributes = c("ensembl_gene_id",
                                    "hgnc_symbol",
                                    "start_position", "end_position",
                                    "description"),
                     filters = "ensembl_gene_id",
                     values = df$geneID,
                     mart = ensembl)
    
    idx <- match(df$geneID, genemap$ensembl_gene_id)
    
    df$symbol <- genemap$hgnc_symbol[idx]
    df$gene_start <- genemap$start_position[idx]
    df$gene_end <- genemap$end_position[idx]
    df$description <- genemap$description[idx]
    df$gene_length <- df$gene_end - df$gene_start + 1
    
  }
  
  if(format == "xlsx"){
    require("openxlsx")
    write.xlsx(df, file = paste0(filename, ".xlsx"), colNames = T, rowNames = F, append = F)
  } else {
    write.csv(df, rowNames = F, file = paste0(filename, ".csv"))
  }
  
  return(df)

}
