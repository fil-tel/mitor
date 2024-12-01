#' fetch_seq
#'
#' @param query NCBI query to download
#' @param dir_path path where to save the files
#' @param n number of sequences to be download
#'
#' @return nothing
#' @export
#'
#' @examples ceee
fetch_seq <- function(query, dir_path, n=NULL){

  search <- entrez_search(db = "nuccore",term = query, use_history = TRUE)

  unlink(dir_path, recursive = TRUE, force = TRUE)
  dir.create(dir_path, recursive = TRUE)

  if (is.null(n)) {
    n=search$count
  }

  else if(search$count<n){
    n=search$count
  }

  # this condition it's needed in the case where we wanted to download a lot of sequences,
  # it's much more reasonable and faster to download them in chunks of fasta file due to limitations of the server,
  # we won't download the genebank file though

  if (n>=5000) {
    for( seq_start in seq(1,n,5000)){
      max=5000
      if(n-(seq_start-1)<5000){
        max=n-(seq_start-1)
      }
      recs <- entrez_fetch(db="nuccore", web_history=search$web_history,
                           rettype="fasta", retmax=max, retstart=seq_start-1)
      write(recs, file=paste0(path, seq_start, "-", seq_start+max-1, ".fa"))
      cat(seq_start+max-1, "sequences downloaded\r")
    }
  }

  else{
    for (i in seq(n)) {
      cat(sprintf("Processing sequence [%d/%d]\n", i, n))

      # fetch a GenBank file in a string format
      gb_txt <- entrez_fetch(db = "nuccore", web_history = search$web_history,
                             rettype = "gb", retmax = 1, retstart = i-1)

      cat("  - Fetched GenBank file from the internet\n")

      # extract accession # from the file for saving the data
      accession_id <-
        gb_txt %>%
        strsplit("\n") %>%
        { grep("ACCESSION", .[[1]], value = TRUE, ignore.case = TRUE) } %>%
        strsplit(" +") %>%
        { .[[1]][2] } # "ACCESSION <accession ID> <potentially other stuff>

      cat(sprintf("  - Extracted accession code %s\n", accession_id))

      gb_file <- file.path(dir_path, paste0(accession_id, ".gb"))
      fa_file <- file.path(dir_path, paste0(accession_id, ".fa"))

      cat(sprintf("  - Saving GenBank file to %s\n", gb_file))

      # write the GenBank to a file
      write(gb_txt, file = gb_file)

      cat(sprintf("  - Converting GenBank to FASTA file %s\n", fa_file))

      # convert to FASTA
      tryCatch(
        gb2fasta(gb_file, fa_file),
        error = function(e) {
          if (e$message == "argument of length 0") {
            cat("  - Empty FASTA string, skipping this record\n")
            unlink(gb_file)
            unlink(fa_file)
          } else {
            stop("Unknown error!", .call = FALSE)
          }
        }
      )
    }
  }
}
