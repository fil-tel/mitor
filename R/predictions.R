process_curl <- function(input_file) {
  curl_content <- readLines(input_file, warn = FALSE)
  url <- gsub("'", "", unlist(strsplit(curl_content, split = " "))[2])

  headers_pattern <- "-H '([^']*)'"
  headers <- regmatches(curl_content, gregexpr(headers_pattern, curl_content))[[1]]
  headers <- gsub("-H '", "", headers)
  headers <- gsub("'$", "", headers)

  post_data <-  gsub("'", "",unlist(strsplit(curl_content, split = "--data-raw "))[2])

  generic <- stringr::str_extract(post_data, "(?<=generic).*")

  headers_dic <- c()
  for (header in headers) {
    split_header <- strsplit(header, ": ", fixed = TRUE)[[1]]
    if (length(split_header) == 2) {
      key <- trimws(split_header[1])
      value <- trimws(split_header[2])
      headers_dic[key] <- value
    }
  }

  return(list(url = url, headers = headers_dic, generic = generic))
}


download_zip <- function(dir, name, download_code, list){
  url <- paste0("https://alphafoldserver.com/_/v2/folds/", download_code, "/fold.zip")

  headers <- httr::add_headers(
    `User-Agent` = as.character(list["User-Agent"]),
    `Accept` = "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    `Accept-Language` = "en-US,en;q=0.5",
    `Accept-Encoding` = "gzip, deflate, br, zstd",
    `Alt-Used` = "alphafoldserver.com",
    `Connection` = "keep-alive",
    `Cookie` = as.character(list["Cookie"]),
    `Upgrade-Insecure-Requests` = "1",
    `Sec-Fetch-Dest` = "document",
    `Sec-Fetch-Mode` = "navigate",
    `Sec-Fetch-Site` = "none",
    `Sec-Fetch-User` = "?1",
    `Priority` = "u=0, i",
    `TE` = "trailers"
  )

  status=500
  while(status!=200){
  Sys.sleep(40)
  response <- httr::GET(url, headers)
  status=httr::status_code(response)
  # print(status)
  if (status == 200) {
    file_path <- paste0(dir, "/" , name, ".zip")
    writeBin(content(response, "raw"), file_path)
    cat(paste0("The file ", name, ".zip has bin downloaded in ", file_path))
  }
  }
}

#' Predicting Protein Structure With AlphaFold3 Server
#'
#' @param seq String. Amino acids sequence of the protein, only canonical AAs are accepted. AAString object are accepted.
#' @param name String. Name of the zip file containing the prediction to be created.
#' @param dir String. Path of the directory where to save the zip file with the prediction. If it doesn't already exist, it will be created.
#'
#' @return The function
#' @export
#'
predict_af3 <- function(seq = NULL, name = NULL, dir = NULL){

  seq <- as.character(seq)
  if(!grepl(paste0("^[", paste0(Biostrings::AA_ALPHABET[1:20], collapse = ""), "]+$"), seq, ignore.case = TRUE)) stop("The AA sequence contains character that are not accepted.")
  if(is.null(seq)) stop("You need to pass a protein sequence!")
  if(is.null(name)) stop("You need to pass a name for yout job!")
  if(is.null(dir)) dir = getwd()

  path <- system.file("extdata", "curl.txt", package = "mitor")
  if(identical(readLines(path), character(0))){
    stop("You didn't provide the your cookies to access the AlphaFold3 server. Follow the instructions here:
         linktogthub")
  }
  else{
    list_post <- process_curl(path)
  }

  post_data <- paste0("f.req=%5B%5B%5B%22kMnDgb%22%2C%22%5B%5B%5C%22", name, "%5C%22%2C%5B%5B%5Bnull%2C%5C%22", toupper(seq), "%5C%22%2Cnull%2Cnull%2Cnull%2Cnull%2C1%5D%5D%5D%5D%2C1%2C%5B2093574080%5D%5D%22%2Cnull%2C%22generic ",list_post$generic)

  response <- httr::POST(
    list_post$url,
    httr::add_headers(.headers = list_post$headers),
    body = post_data,
    encode = "raw"
  )

  if(response$status_code==200) cat("The call to the server went well! Wait sometime for the download!\n")
  else stop(paste0("OPS! Something went wrong!", "Error:", status_code(response), "\n"))

  download_code <- unlist(strsplit(rawToChar(response$content), "\\\\\\\""))[2]
  cat(download_code)
  dir.create(dir, recursive = TRUE)
  future::plan(future::multisession)
  future::future(download_zip(dir = dir, name = name, download_code = download_code, list = list_post$headers), packages = "httr")
}

sequence <- paste0(sample(Biostrings::AA_ALPHABET[1:20], size = 120, replace = TRUE), collapse = "")
predict_af3(seq = sequence, name = "bobob", dir = "boia/cazz")



