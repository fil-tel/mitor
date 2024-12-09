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
      Biostrings::writeXStringSet(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end), paste0(rrna_path, "/", coord_df[i,]$gene_id ,'.fa'))
    }else if (startsWith(coord_df[i,]$gene_id, "T")) {
      Biostrings::writeXStringSet(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end), paste0(trna_path, "/", coord_df[i,]$gene_id ,'.fa'))
    }else{
      # ND6 is the only gene on the - strand, since later I will translate them I reversed it now
      if(coord_df[i,]$strand == "+"){
        Biostrings::writeXStringSet(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end), paste0(cds_path, "/", coord_df[i,]$gene_id ,'.fa'))
      }
      else{
        Biostrings::writeXStringSet(Biostrings::reverseComplement(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end), paste0(cds_path, "/", coord_df[i,]$gene_id ,'.fa')))
      }
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
      rrna_list <- append(rrna_list, Biostrings::DNAStringSetList(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end)))
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
      trna_list <- append(trna_list, Biostrings::DNAStringSetList(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end)))
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
      # ND6 is the only gene on the - strand, since later I will translate them I reversed it now
      if(coord_df[i,]$strand == "+"){
        cds_list <- append(cds_list, Biostrings::DNAStringSetList(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end)))
      }
      else{
        cds_list <- append(cds_list, Biostrings::DNAStringSetList(Biostrings::reverseComplement(Biostrings::subseq(msa, start = coord_df[i,]$start, end = coord_df[i,]$end))))
      }

      cds_names <- c(cds_names, coord_df[i,]$gene_id)
    }
  }
  names(cds_list) <- cds_names
  cds_list
}


#' Translate MSA Corresponding to a CDS
#'
#' This function translates a given MSA of mitochondrial CDS into the corresponding AA sequences.
#'
#' @param msa A multiple sequence alignment object of mitochondrial CDS, a `DNAStringSet`
#'            from the `Biostrings` package. Typically a DNAStringSet object contained in the list
#'            generated by the function \code{extract_cds}.
#'
#' @return  A `DNAStringSet` from the `Biostrings` package. The result will not necessarily be a MSA, since gaps are removed in the translation.
#' @export
#'
#' @note
#' CDS MSA containing gaps are not reccomended for this function.
translation <- function(msa){
  # ungapped msa, the presence of gaps messes up the translation
  ungapped <- Biostrings::DNAStringSet(gsub("-", "", msa))
  widths <- Biostrings::width(ungapped)
  n_codons <- numeric(length(widhts))
  # number of codons excluding stop codon, check if integer because of polyA stop codon
  for (i in seq(length(n_codons))) {
    if(widths[i]%%3==0){
      n_codons[i] <- widths[i]/3-1
    }
    else{
      n_codons[i] <- as.integer(widths[i]/3)
    }
  }
  # extract cds without stop codon
  prot2trans <- Biostrings::subseq(ungapped, start = 1, end = n_codons*3)
  # translate
  Biostrings::translate(prot2trans, genetic.code = getGeneticCode("2"), if.fuzzy.codon = "X")
}


#' Check Presence of Gaps in a MSA
#'
#' @param msa A multiple sequence alignment object of mitochondrial CDS, a `DNAStringSet`
#'            from the `Biostrings` package.
#'
#' @return Boolean. TRUE if gaps are presents, FALSE if not.
#' @export
#'
#' @examples
#'
#' aln <- Biostrings::DNAStringSet(c( "ATGGCATCTACTTTGTATGACTATTGCAAGTGCCCATGGGTGACATCTGTAAGAAAGATGGGGATAAGCGCTGTAAGCTT", "ATGGCATCTACTTCGTATGACTATTGCAGAGTGCCCATG--------------GAAGACGGGGATAAGCGCTGTAAGCTT", "ATGGCATCTACTTTGTATGACTATTGCAGAGTGCCCATGGGTGACATCTGTAGAAAGATGGGGATAAGCGCTGTAAGCTT"))
#' # Check for gaps
#' check_gaps(aln)
#'
check_gaps <- function(msa){
  !sum(nchar(gsub("-", "", msa)))==sum(nchar(msa))
}


#' Identify Sequence Variations in Multiple Sequence Alignment (MSA)
#'
#' This function detects and returns variations between a target sequence and a reference sequence (typically rCRS)
#' within a given multiple sequence alignment (MSA). Variations are annotated following the
#' nomenclature described at <https://www.hgmd.cf.ac.uk/docs/mut_nom.html>. Note that different
#' nomenclature systems may appear in other studies.
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
#'            from the `Biostrings` package. Containing the target and the reference.
#' @param target String. Name of the target sequence in contained in the MSA.
#' @param ref String. Name of the reference sequence in the MSA.
#'            Defaults to `"NC_012920"`.
#'
#' @return A character vector containing the identified variations between the target and reference
#'         sequences. Variations are described using the following formats:
#'         - Substitutions: `<position><reference_base>><target_base>` (e.g., "123A>G").
#'         - Insertions: `<start_position>_<end_position>ins<bases>` (e.g., "123_124insAT").
#'         - Deletions: `<start_position>_<end_position>del<bases>` (e.g., "123_125delGTC").
#'
#' @details
#' Positions are adjusted to match the reference sequence, accounting for gaps introduced
#'          during the alignment process. Ambiguous bases (e.g., "N") are ignored.
#'
#' @export
#'
#' @examples
#' aln <- Biostrings::DNAStringSet(c( "ATGGCATCTACTTTGTATGACTATTGCAAGTGCCCATGGGTGACATCTGTAAGAAAGATGGGGATAAGCGCTGTAAGCTT", "ATGGCATCTACTTCGTATGACTATTGCAGAGTGCCCATG--------------GAAGACGGGGATAAGCGCTGTAAGCTT", "ATGGCATCTACTTTGTATGACTATTGCAGAGTGCCCATGGGTGACATCTGTAGAAAGATGGGGATAAGCGCTGTAAGCTT"))
#' names(aln) <- paste0("seq", 1:3)
#' find_variants(aln, target = "seq2", ref = "seq1")
find_variants <- function(msa, target, ref="NC_012920"){
  # ref seq goes in front
  if(!target %in% names(msa)) stop("The target sequence is not present in the MSA. Check the argument you passed.")
  if(!ref %in% names(msa)) stop("The reference sequence is not present in the MSA. Check the argument you passed.")
  seqs <- as.matrix(msa[c(ref,target)])
  variants_list <- apply(seqs, 2, unique)
  gaps <- unlist(gregexpr("-",as.character(msa[[ref]])))
  mutations <- c()
  i <- 1
  # variants_list
  while (i<=length(variants_list)) {
    # if in the position there is only 1 character (no difference) skip the iteration
    if(length(variants_list[[i]])==1 |  "N" %in% variants_list[[i]]){
      i <- i+1
      next
    }
    # Insertion case: if there is a gap in the reference
    else if(variants_list[[i]][1]=="-"){
      # position after the beginning of the gap
      if(i+1<=length(variants_list)){
        c <- i+1
      }
      # if the insertion is longer than one base we need to check how long it is
      while (variants_list[[c]][1]=="-" & c+1<=length(variants_list) & length(variants_list[[c]])==2) {
        c <- c+1
      }
      # adjust the insertion indices according to the position in the ref seq not in the msa, as always
      ins_start <- i-sum(gaps<i)-1
      ins_end <- ins_start+1
      # extract the bases inserted in the target seq
      ins_bases <- paste0(sapply(variants_list[i:(c-1)], function(x) x[2]), collapse = "")
      mutations <- c(mutations, paste0(ins_start, "_", ins_end, "ins", ins_bases))
      # this is needed in case we have a gap in the last position, otherwise if don't do this we get stuck in the loop
      if(!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1>length(variants_list)){
        i <- c+1
      }
      else{
        i <- c
      }
    }
    # Deletion case: if there is a gap in the target
    else if(variants_list[[i]][2]=="-"){
      # position after the beginning of the gap
      if(i+1<=length(variants_list)){
        c <- i+1
      }
      # if the deletion is longer than one base we need to check how long it is
      # the other conditions are to avoid to go out of bound when if at the end of the alignment
      while (!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1<=length(variants_list) & length(variants_list[[c]])==2) {
        c <- c+1
      }
      # adjust the insertion indices according to the position in the ref seq not in the msa, as always
      del_start <- i-sum(gaps<i)
      del_end <- del_start+(c-i)-1
      # one base del
      if(del_start==del_end){
        del_bases <- variants_list[[i]][1]
        mutations <- c(mutations, paste0(del_start, "del", del_bases))
      }
      # more than one base deletion
      else{
        # extract the bases inserted in the target seq
        del_bases <- paste0(sapply(variants_list[i:(c-1)], function(x) x[1]), collapse = "")
        mutations <- c(mutations, paste0(del_start, "_", del_end, "del", del_bases))
      }
      # this is needed in case we have a gap in the last position, otherwise if don't do this we get stuck in the loop
      if(!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1>length(variants_list)){
        i <- c+1
      }
      else{
        i <- c
      }
    }
    # subsituition case
    else{
      # adjust index according to the msa
      sub_i <- i-sum(gaps<i)
      mutations <- c(mutations, paste0(sub_i, variants_list[[i]][1], ">", variants_list[[i]][2]))
      i <- i+1
    }
  }
  mutations
}

#' Find Non Conserved Columns in a Multiple Sequence Alignment (MSA)
#'
#' This function finds and extracts the non conserved columns in a MSA.
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
#'            from the `Biostrings` package.
#' @param type String. Type of MSA, amino acids ("AA") or nucleotides ("DNA").
#'
#' @return Matrix. The function returns a matrix that looks like this. The names of the columns coresponds to the
#' position in the MSA, while the rownames to the name of the sequence.
#'
#'      14  29  30  31  32  33  36  37  38  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  59
#' seq1 "T" "A" "G" "T" "G" "C" "A" "T" "G" "G" "T" "G" "A" "C" "A" "T" "C" "T" "G" "T" "A" "A" "G" "A" "T"
#' seq2 "C" "G" "A" "G" "T" "G" "C" "A" "T" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "G" "C"
#' seq3 "T" "G" "A" "G" "T" "G" "C" "A" "T" "G" "G" "T" "G" "A" "C" "A" "T" "C" "T" "G" "T" "A" "G" "A" "T"
#'
#' @export
#'
#' @examples
#' aln <- Biostrings::DNAStringSet(c( "ATGGCATCTACTTTGTATGACTATTGCAAGTGCCCATGGGTGACATCTGTAAGAAAGATGGGGATAAGCGCTGTAAGCTT", "ATGGCATCTACTTCGTATGACTATTGCAGAGTGCCCATG--------------GAAGACGGGGATAAGCGCTGTAAGCTT", "ATGGCATCTACTTTGTATGACTATTGCAGAGTGCCCATGGGTGACATCTGTAGAAAGATGGGGATAAGCGCTGTAAGCTT"))
#' names(aln) <- paste0("seq", 1:3)
#' find_var_pos(aln, type="DNA")
#'
find_var_pos <- function(msa, type="AA"){
  msa_matrix <- as.matrix(msa)

  # find positions of non conserved columns
  if(type=="AA"){
    mismatch_positions <- which(apply(msa_matrix, 2, function(col) length(unique(col)) > 1 & !(any(unique(col)=="X") & length(unique(col))==2)))
  }
  else if(type=="DNA"){
    mismatch_positions <- which(apply(msa_matrix, 2, function(col) length(unique(col)) > 1 & !(any(unique(col)=="N") & length(unique(col))==2)))
  }
  else{stop("Argument type not available!")}

  # return something only if the msa is not totally conserved, otherwise is NULL
  if(!identical(mismatch_positions, integer(0))){
    # extract only non conserved columns
    non_cons_cols <- msa_matrix[, mismatch_positions, drop=FALSE]
    colnames(non_cons_cols) <- mismatch_positions
    non_cons_cols
  }

}
