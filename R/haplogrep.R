#' Assign the Haplogroups of a set of human mtDNA Sequences
#'
#' `classify_haplogroup` classify mitochondrial DNA haplogroups using [Haplogrep 3](https://haplogrep.readthedocs.io/en/latest/installation/)
#'
#'
#' @param sequences \code{DNAStringSet}. Set of human mitochondrial DNA sequences to classify.
#'
#' @return A table.
#' \describe{
#'   \item{SampleID}{The ID of the sequence.}
#'   \item{Haplogroup}{The predicted haplogroup.}
#'   \item{Quality}{The quality score of the prediction.}
#'   \item{Range}{The range of the sequence used for prediction.}
#' }
#'
#' Output example:
#' \preformatted{
#'       SampleID   Haplogroup Rank Quality             Range
#' 1   DQ856316.1        HV1b2    1  1.0000           1-16569
#' 2   EU603401.1        K1a4a    1  0.9272           1-16569
#' 3   KX681447.1          H7e    1  0.9610           1-16569
#' 4   KX702230.1           H7    1  0.8705           1-16569
#' }
#' @export
#' @details
#' `classify_haplogroup` invokes \emph{haplogrep3}, which has to be previously installed on the user's machine.
#' For a step-by-step guide on installing and setting up haplogrep3, visit [github]().
#'
#' The function might take 30-40 min for large datasets (~60000 sequences).
#'
#' For more information regarding the software and the results visit \url{https://haplogrep.readthedocs.io/en/latest/installation/}
#'
#'
classify_haplogroup <- function(sequences) {
  if (Sys.which("haplogrep3") == "")
    stop(
      "The haplogrep binary couldn't be found. \n Did you add it to your PATH variable? \n
  Visit https://github.com/fil-tel/mitor for a step-by-step guide.                                                                   Visit githublink for a step-by-step guide."
    )

  # if sequences have IUPAC ambiguities not accepted by haplogrep they will be removed
  if (any(grepl("[URYSWKMBDHV]", sequences))) {
    cat(
      "Warning! Some of the sequences contain IUPAC ambiguities not accepted by haplogrep3.\n"
    )
    cont <-
      readline('Do you wish to continue? [Y/N] \n If yes these sequences will be removed.')
    if (cont != 'Y')
      stop('Aborted by user')
    sequences <- sequences[!grepl("[URYSWKMBDHV]", sequences)]
  }
  # temporary files for input and output
  input_file <- tempfile(fileext = ".fasta")
  output_file <- tempfile(fileext = ".tsv")
  Biostrings::writeXStringSet(sequences, filepath = input_file)

  # run haplogrep
  system(
    paste(
      "haplogrep3 classify --in",
      input_file,
      "--tree phylotree-rcrs@17.2" ,
      "--out",
      output_file
    ),
    wait = TRUE
  )

  result <- read.table(output_file, header = TRUE)

  # remove temporary files
  unlink(input_file)
  unlink(output_file)

  return(result)
}

