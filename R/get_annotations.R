############################
# Function get_annotations #
############################

#' This function annotates a column of transcripts IDs (ENSEMBL) with information of the genes to which they belong.
#'
#' The Gene information added include:
#' - ENSEMBL ID
#' - HGNC Symbol and description
#' - Start and End position
#' - Gene length
#'
#' @param tx_ensembl_ids The column of transcripts to be used as input.
#' @param version This function can use the version 103 or the current version of the Biomart. Default = 103.
#' @param format The output will be a table named tx2gene in .csv or .xlsx formats. Default = csv.
#' @export

get_annotations <- function(tx_ensembl_ids, version = "103", format = "csv") {

  require("biomaRt")
  require("ensembldb")

  tx2gene <- data.frame(transcriptIDs = tx_ensembl_ids)

  if(version == "103"){
    ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl",
                      host = "feb2021.archive.ensembl.org")
  } else {
    ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl",
                      host = "useast.ensembl.org")
  }

  # More annotations available in the biomaRt:
  #listAttributes(mart = ensembl)

  genemap <- getBM(attributes = c("ensembl_transcript_id_version",
                                  "ensembl_gene_id",
                                  "hgnc_symbol",
                                  "start_position", "end_position",
                                  "description"),
                   filters = "ensembl_transcript_id_version",
                   values = tx2gene$transcriptIDs,
                   mart = ensembl)

  idx <- match(tx2gene$transcriptIDs, genemap$ensembl_transcript_id_version)

  tx2gene$geneID <- genemap$ensembl_gene_id[idx]
  tx2gene$symbol <- genemap$hgnc_symbol[idx]
  tx2gene$gene_start <- genemap$start_position[idx]
  tx2gene$gene_end <- genemap$end_position[idx]
  tx2gene$description <- genemap$description[idx]
  tx2gene$gene_length <- tx2gene$gene_end - tx2gene$gene_start + 1

  # If you want to annotate the PRKACA-DNAJB1 chimera (if quantified),
  # use these lines and be careful to not leave NA values in geneID
  # tx2gene$geneID[tx2gene$transcriptIDs == "PRKACA-DNAJB1.1"] <- "PRKACA-DNAJB1"
  # tx2gene$symbol[tx2gene$transcriptIDs == "PRKACA-DNAJB1.1"] <- "PRKACA-DNAJB1"
  # tx2gene$description[tx2gene$transcriptIDs == "PRKACA-DNAJB1.1"] <- "Chimeric protein connecting exons 2-10 of PRKACA and exon 1 of DNAJB1"

  if(format == "xlsx"){
    install.packages("openxlsx")
    write.xlsx(tx2gene, file = "tx2gene.xlsx", col.names = T, row.names = F, append = F)
  } else {
    write.csv(tx2gene, row.names = F, file = "tx2gene.csv")
  }

  return(tx2gene)
}
