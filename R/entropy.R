#' Calculate Shannon Entropy
#'
#' This function apply the formula for the Shannon entropy for a vector of probabilities (!!!not absolute frequencies!!!).
#' \deqn{H = -\sum_{i=1}^n p_i \log_2(p_i)}
#' where:
#' - \eqn{H} is the Shannon entropy.
#' - \eqn{p_i} is the proportion of occurrences of the \eqn{i}-th category (e.g., nucleotide or amino acid frequency).
#' - The logarithm is base 2.
#'
#' Shannon entropy is highest when all categories are equally likely and is lowest (zero) when one category dominates completely.
#'
#' @param prob Vector. Vector of probabilities (ex. frequencies of nucleotides in a column of a MSA)
#'
#' @return Double. Shannon entropy value.
#' @export
#'
calc_entropy <- function(prob){
  # consider only non zero frequencies
  prob <- prob[prob>0]
  abs(sum(prob*log2(prob)))
}

#' Find Shannon Entropy for each position of a MSA
#'
#' @param msa A multiple sequence alignment object, a `DNAStringSet`
#'            from the `Biostrings` package.
#'
#' @return A vector with the entropy value for each position of the MSA.
#' @export
#'
#'
position_entropy_dna <- function(msa){
  # find consensus matrix as counts and calc the entropy for dna seq
  cons_matrix <- Biostrings::consensusMatrix(msa)
  freq <- cons_matrix[c("A", "C", "G", "T", "-"),]
  freq <- sweep(freq, 2, colSums(freq), "/")
  entropy <- apply(freq, 2, calc_entropy)
  entropy
}

#' Find Shannon Entropy for each position of a MSA of Proteins
#'
#' @param msa A multiple sequence alignment object, a `AAStringSet`
#'            from the `Biostrings` package.
#'
#' @return A vector with the entropy value for each position of the MSA.
#' @export
#'
#'
position_entropy_AA <- function(msa){
  # find consensus matrix as counts and calc the entropy for proteins seq
  cons_matrix <- Biostrings::consensusMatrix(msa)
  # consider only the 20 amino acid, no ambiguities
  freq <- cons_matrix[1:20,]
  freq <- sweep(freq, 2, colSums(freq), "/")
  entropy <- apply(freq, 2, calc_entropy)
  entropy
}
