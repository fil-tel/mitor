#' Adjust Gene Coordinates Based on Gaps in a Multiple Sequence Alignment
#'
#' This function adjusts the coordinates of genes from the revised Cambridge Reference Sequence (rCRS)
#' to account for the positions of gaps ("-") in a reference sequence within a multiple sequence alignment (MSA).
#'
#' @param msa A multiple sequence alignment object, typically a `DNAStringSet`
#'            from the `Biostrings` package. The reference sequence in the MSA should be identified
#'            as `NC_012920`.
#'
#' @return A data frame containing the adjusted gene coordinates:
#' \describe{
#'   \item{\code{seqnames}}{Character. The origin of the sequence: "Homo sapiens".}
#'   \item{\code{start}}{Integer. The starting coordinate of the gene in the MSA.}
#'   \item{\code{end}}{Integer. The ending coordinate of the gene in the MSA.}
#'   \item{\code{width}}{Integer. The length of the gene in base pairs.}
#'   \item{\code{strand}}{Character. The strand orientation of the gene: \code{"+"} or \code{"-"}.}
#'   \item{\code{type}}{Character. \code{"gene"} for all.}
#'   \item{\code{gene}}{Character. The gene identifier (e.g., \code{TRNF}, \code{RNR1}).}
#'   \item{\code{nomenclature}}{Character. Nomenclature annotations including the official symbol, full name, and source of the information.}
#'   \item{\code{db_xref}}{Character. Database cross-references for the gene.}
#'   \item{\code{loctype}}{Character. Location type, \code{"normal"} for all.}
#'   \item{\code{gene_synonym}}{Character. Synonyms for the gene, if any.}
#'   \item{\code{gene_id}}{Character. The unique gene identifier.}
#' }
#'
#' @details
#' The function identifies gaps ("-") in the reference sequence \code{NC_012920} within the provided MSA
#' and adjusts the start and end positions of the genes accordingly. This ensures that the gene coordinates
#' remain accurate relative to the alignment.
#'
#' @note
#' The input MSA must have the reference sequence named \code{NC_012920} for the function to correctly
#' identify and adjust the gaps. Before running the MSA, it is essential to include the revised Cambridge reference sequence
#' \code{rCRS} already provided in the package.
#'
#'
#'
#' @export
genes_coord <- function(msa){
  genes <- rCRS_genes_df
  # identify locations of the gaps in the reference to adjust the indices of the genes
  gaps <- unlist(gregexpr("-",as.character(msa$`NC_012920`)))
  # adjust the indices according to the position of the gaps in the reference
  for (i in gaps) {
    genes$start <- genes$start+as.integer(genes$start>=i)
    genes$end <- genes$end+as.integer(genes$end>=i)
  }
  genes
}


