#' Cambridge Reference Sequence Gene Coordinates
#'
#' This dataset contains the coordinates of genes in the revised Cambridge reference sequence of the human mitochondrial genome.
#' It includes details such as gene locations, strand orientation, type, and nomenclature annotations.
#'
#' @format A data frame with 37 rows and 12 columns:
#' \describe{
#'   \item{\code{seqnames}}{Character. The origin of the sequence: "Homo sapiens".}
#'   \item{\code{start}}{Integer. The starting coordinate of the gene on the sequence.}
#'   \item{\code{end}}{Integer. The ending coordinate of the gene on the sequence.}
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
#' This dataset provides detailed annotations of the genes present in the Cambridge reference sequence, which represents
#' the human mitochondrial genome.
#'
#' @source [https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/)
#'
#'
"rCRS_genes_df"

#' Cambridge Reference Sequence - Human Mitochondrial Genome
#'
#' This object contains the full DNA sequence of the revised Cambridge Reference Sequence (rCRS), representing the human mitochondrial genome.
#'
#' @format A \code{DNAStringSet} object of length 1:
#' \describe{
#'   \item{\code{width}}{Integer. The total length of the DNA sequence, 16,569 base pairs.}
#'   \item{\code{seq}}{Character. The complete nucleotide sequence of the human mitochondrial genome.}
#'   \item{\code{names}}{Character. The accession number and identifier of the sequence, \code{"NC_012920"}.}
#' }
#'
#' @details
#' It serves as the reference for sequence alignment, annotation, and comparative studies in mitochondrial genetics.
#' For more details check this [website](https://www.mitomap.org/foswiki/bin/view/MITOMAP/HumanMitoSeq)
#'
#' @source [https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/)
#'
#'
"rCRS"

