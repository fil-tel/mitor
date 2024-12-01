#' Fetching Sequences From The NCBI Database
#'
#' `fetch_seq()` fetches sequences from the specified NCBI database and save them in a specified directory on your machine.
#' The sequences will be downloaded and stored in a FASTA files with maximum 5000 sequences in each. If the argument `gb` is passed the genbank file will also be
#' downloaded.
#'
#' @param query String. NCBI query to search.
#'
#'  A documentation for how NCBI querys are structured is available [here](https://www.ncbi.nlm.nih.gov/books/NBK44863/).
#'  It is helpful to first manually try different querys on the [NCBI websites](https://www.ncbi.nlm.nih.gov/) to understand which suites better for your goal.
#' @param dir_path String. Path of the directory where to save the files.
#'
#'  It is essential to specify whether the directory is new or it already exists using the argument `create`.
#' @param n Integer. Number of sequences to be downloaded. DEFAULT: all the records corresponding to the query.
#'
#'  * If n is smaller than the number of records corresponding to the query, the first n records will be downloaded.
#'  * If n is greater than the number of records corresponding to the query, all the records will be downloaded.
#' @param create Boolean that specify whether to create a new directory (TRUE) or to use an already existence one (FALSE).
#'
#'  IMPORTANT: If TRUE, but the directory already exist, it will deleted with all its content and recreated.
#' @param db String. NCBI database from where the sequences have to be downloaded.
#'
#'  The available databases can be found [here](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html).
#' @param gb Boolean. If TRUE, the the genbank file will be download together with the FASTA.
#' @param filename String. Name you want to give to the file.
#'
#' @examples
#' ref_query <- "NC_012920"
#' ref_path <- "data/mtdna/ref"
#' fetch_seq(ref_query, ref_path, filename = "ref", gb = TRUE)
#'
#'
#' @export
#'
fetch_seq <- function(query, dir_path, filename, create=TRUE, n=NULL, db = "nuccore", gb = FALSE){

  if(create && file.exists(dir_path)){
    cat("Warning! The directory alredy exists. All its content will be deleted.", '\n')
    cont <- readline('Do you wish to continue? [y/n] ')
    if(cont != 'y') stop('Aborted by user')
  }

  search <- entrez_search(db = db,term = query, use_history = TRUE)

  if (create) {
    unlink(dir_path, recursive = TRUE, force = TRUE)
    dir.create(dir_path, recursive = TRUE)
  }


  if (is.null(n)) {
    n=search$count
  }

  else if(search$count<n){
    n=search$count
  }

  for (seq_start in seq(1, n, 5000)) {
    max = 5000
    if (n - (seq_start - 1) < 5000) {
      max = n - (seq_start - 1)
    }
    recs_fa <-
      entrez_fetch(
        db = db,
        web_history = search$web_history,
        rettype = "fasta",
        retmax = max,
        retstart = seq_start - 1
      )
    write(recs_fa, file = paste0(dir_path,  "/", filename, ".fa"), append = TRUE)
    if(gb){
      recs_gb <-
        entrez_fetch(
          db = db,
          web_history = search$web_history,
          rettype = "gb",
          retmax = max,
          retstart = seq_start - 1
        )
      write(recs_gb, file = paste0(dir_path, "/", filename, ".gb"), append = TRUE)
    }
    cat(seq_start + max - 1, "sequences downloaded\r")
  }

}
