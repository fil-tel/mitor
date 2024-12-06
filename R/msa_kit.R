#' Adjust Gene Coordinates Based on Gaps in a Multiple Sequence Alignment
#'
#' This function adjusts the coordinates of genes from the revised Cambridge Reference Sequence (rCRS)
#' to account for the positions of gaps ("-") in a reference sequence within a multiple sequence alignment (MSA).
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
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

#' Save Genes from a mtDNA Multiple Sequence Alignment in FASTA files
#'
#' This function extracts the 37 mitochondrial genes (protein-coding genes, tRNAs, and rRNAs)
#' from a multiple sequence alignment (MSA) and saves them as FASTA files organized by type.
#' Each FASTA file will correspond to the subset of the MSA containing a specific gene.
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
#'            from the `Biostrings` package.
#' @param coord_df Data frame. A data frame created using the function \code{genes_coord}.
#' @param dir String. Name of the directory where to save the extracted genes.
#'
#'
#' @export
#'
save_genes <- function(msa, coord_df, dir="data"){

  cds_path <- paste0(dir, "/genes/cds")
  trna_path <- paste0(dir,"/genes/trna")
  rrna_path <- paste0(dir,"/genes/rrna")

  unlink(cds_path, recursive = TRUE, force = TRUE)
  dir.create(cds_path, recursive = TRUE)
  unlink(trna_path, recursive = TRUE, force = TRUE)
  dir.create(trna_path, recursive = TRUE)
  unlink(rrna_path, recursive = TRUE, force = TRUE)
  dir.create(rrna_path, recursive = TRUE)

  for (i in rownames(coord_df)) {
    # save in fasta the subsets of the alignment containing only the genes
    if (startsWith(coord_df[i,]$gene_id, "R")) {
      Biostrings::writeXStringSet(subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end), paste0(rrna_path, "/", coord_df[i,]$gene_id ,'.fa'))
    }else if (startsWith(coord_df[i,]$gene_id, "T")) {
      Biostrings::writeXStringSet(subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end), paste0(trna_path, "/", coord_df[i,]$gene_id ,'.fa'))
    }else{
      Biostrings::writeXStringSet(subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end), paste0(cds_path, "/", coord_df[i,]$gene_id ,'.fa'))
    }
  }
}

#' Extract rRNA Genes from a Mitochondrial DNA Multiple Sequence Alignment
#'
#' This function extracts ribosomal RNA (rRNA) genes from a multiple sequence alignment (MSA)
#' based on the coordinates provided in a data frame. The extracted rRNA sequences are returned
#' as a named `DNAStringSetList`.
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
#'            from the `Biostrings` package.
#' @param coord_df Data frame. A data frame created using the function \code{genes_coord}.
#'
#' @return A `DNAStringSetList` object containing the extracted rRNA gene sequences,
#'         with gene identifiers as names.
#' @export
#'
extract_rrna <- function(msa, coord_df){
  rrna_list <- Biostrings::DNAStringSetList()
  rrna_names <- c()
  for (i in rownames(coord_df)) {
    if (startsWith(coord_df[i,]$gene_id, "R")) {
      rrna_list <- append(rrna_list, Biostrings::DNAStringSetList(subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end)))
      rrna_names <- c(rrna_names, coord_df[i,]$gene_id)
    }
  }
  names(rrna_list) <- rrna_names
  rrna_list
}

#' Extract tRNA Genes from a Mitochondrial DNA Multiple Sequence Alignment
#'
#' This function extracts transfer RNA (tRNA) genes from a multiple sequence alignment (MSA)
#' based on the coordinates provided in a data frame. The extracted tRNA sequences are returned
#' as a named `DNAStringSetList`.
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
#'            from the `Biostrings` package.
#' @param coord_df Data frame. A data frame created using the function \code{genes_coord}.
#'
#' @return A `DNAStringSetList` object containing the extracted tRNA gene sequences,
#'         with gene identifiers as names.
#' @export
#'
extract_trna <- function(msa, coord_df){
  trna_list <- Biostrings::DNAStringSetList()
  trna_names <- c()
  for (i in rownames(coord_df)) {
    if (startsWith(coord_df[i,]$gene_id, "T")) {
      trna_list <- append(trna_list, Biostrings::DNAStringSetList(subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end)))
      trna_names <- c(trna_names, coord_df[i,]$gene_id)
    }
  }
  names(trna_list) <- trna_names
  trna_list
}

#' Extract CDS Genes from a Mitochondrial DNA Multiple Sequence Alignment
#'
#' This function extracts coding sequences (CDS) genes from a multiple sequence alignment (MSA)
#' based on the coordinates provided in a data frame. The extracted CDS are returned
#' as a named `DNAStringSetList`.
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
#'            from the `Biostrings` package.
#' @param coord_df Data frame. A data frame created using the function \code{genes_coord}.
#'
#' @return A `DNAStringSetList` object containing the extracted CDS gene sequences,
#'         with gene identifiers as names.
#' @export
#'
extract_cds <- function(msa, coord_df){
  cds_list <- Biostrings::DNAStringSetList()
  cds_names <- c()
  for (i in rownames(coord_df)) {
    if (!startsWith(coord_df[i,]$gene_id, "T")&!startsWith(coord_df[i,]$gene_id, "R")) {
      cds_list <- append(cds_list, Biostrings::DNAStringSetList(subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end)))
      cds_names <- c(cds_names, coord_df[i,]$gene_id)
    }
  }
  names(cds_list) <- cds_names
  cds_list
}




