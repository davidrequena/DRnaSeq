########################
# Add gene annotations #
########################

#' A function to add annotations to a table of gene counts.
#'
#' @param object A table of gene counts (rows: genes, columns:samples).
#' @param key The variable name (column) in common between the object and the reference, which contains the indentifiers. Default: "geneID".
#' @param reference A reference table with the annotations. Default: geneID.details.
#' @param variables Annotations (columns) from the reference table to add, i.e. gene symbol, gene_length, description.
#' @export

add_annotations <- function(object, reference = geneID.details, key = "geneID", variables = NULL){

  df <- as.data.frame(object)
  df[,key] <- rownames(df)
  index <- match(df[,key], geneID.details[,key])

  # Add a new column for each variable
  for(i in 1:length(variables)){
    df[,variables[i]] <- geneID.details[index, variables[i]]
  }

  return(df)
}
