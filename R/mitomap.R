
# Function to split the tables returned by
# the search allele in mitomap into different dataframe
split_df <- function(df){
  rows2split <- which(df[,1]=="")-1
  rows2split <- c(rows2split, length(df[,1]))
  df_list <- list()
  c <- 1
  for (i in seq(length(rows2split))) {
    df_list[[df[c,1]]] <- df[c:rows2split[i], ]
    c <- rows2split[i]+2
  }

  for (j in seq_along(df_list)){
    tmp <- df_list[[j]]
    df_list[[j]] <- tmp[-c(1,2), ]
    colnames(df_list[[j]]) <- tmp[2,]
    df_list[[j]] <- df_list[[j]][, !is.na(colnames(df_list[[j]]))]
  }
  # I will discard the Somatic Variants dataframe
  df_list[names(df_list)!="MITOMAP: mtDNA Somatic Variants"]
}


#' Search Positions Information in the Mitomap Database
#'
#' search_allele allows the user to programmatically retrieve information regarding human mitochondrial varaints at a given position as shown in the
#' MITOMAP database. It is equivalent to manually searching it in [webserver](https://www.mitomap.org/allelesearch.html),
#' but with the advantage of returning results in a structured data frame format, making it easier to analyze.
#'
#' @param pos Vector. Numeric vector of mitochondrial positions (max. 100).
#'
#' @return A list of data frames. Each data frame corresponds to the variants in Coding Region Sequence Variants,
#' Control Region Sequence Variants, Reported Mitochondrial DNA Base Substitution Diseases.
#' @export
#'
#' @examples
#' \dontrun{
#'   pos <- 13800:13900
#'  example <- search_allele(pos)
#'  test
#'   }
search_allele <- function(pos=NULL){

  if(is.null(pos)) stop("Error! Provide at least 1 position!")
  if(length(pos)>101) stop("Attention! The number of positions is greater than 100.")
  if(max(pos)>16569 | min(pos)<1) stop("Attention! Position cannot be greater than the length of the rCRS (16569) and less than 1.")

  pstn <- paste0(pos, collapse = "+")

  # URL for the search allele db of mitomap
  url <- "https://www.mitomap.org/cgi-bin/search_allele"
  # headers for the request (taken from my call from my laptop, should anyway work always)
  # usually the server assigns also a temporary cookie, apparently, as for the keys it is not
  # really important
  headers <- c(
    "User-Agent" = "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:134.0) Gecko/20100101 Firefox/134.0",
    "Accept" = "*/*",
    "Accept-Language" = "en-US,en;q=0.5",
    "Accept-Encoding" = "gzip, deflate, br, zstd",
    "Content-Type" = "application/x-www-form-urlencoded; charset=UTF-8",
    "X-Requested-With" = "XMLHttpRequest",
    "Origin" = "https://www.mitomap.org",
    "Connection" = "keep-alive",
    "Referer" = "https://www.mitomap.org/allelesearch.html",
    "Sec-Fetch-Dest" = "empty",
    "Sec-Fetch-Mode" = "cors",
    "Sec-Fetch-Site" = "same-origin",
    "Priority" = "u=0"
  )

  # keys to access the server, generally these are assigned by the server once you connect the first time
  # I realized that even two simple keys as the one below work
  form_data <- list(
    validation_key = "22222222222222222222222222222222",
    key = "1111111111111",
    pstns = pstn
  )
  # classic POST request
  response <- httr::POST(
    url,
    httr::add_headers(.headers = headers),
    body = form_data,
    encode = "form"
  )

  content <- httr::content(response, as = "text", encoding = "UTF-8")

  tab <- rvest::html_table(rvest::read_html(content))
  if(length(tab)==0 | any(unique(colnames(tab[[1]]))=="MITOMAP: mtDNA Somatic Variants")) stop("No reported variants are reported in the MITOMAP databse
                                                                    for the given positions.")

  # modifying the html because scraping the table from it
  # it is a bit problematic
  tmp <- as.data.frame(rvest::html_table(rvest::read_html(content))[[1]])
  tmp2 <- gsub(">\\d+\\.\\d+%", "", content)
  tmp3 <- gsub("\\(count", "[count", tmp2)
  tmp4 <- gsub("<br>\\(", "</td><td>", tmp3)
  def <- as.data.frame(rvest::html_table(rvest::read_html(tmp4))[[1]])

  if("Conservation" %in% tmp[2,]){
    cons_col <- tmp[ ,which(tmp[2,]=="Conservation")]
    def[[which(def[2,]=="Conservation")]] <- cons_col
  }

  splitted <- split_df(def)

  # format values
  for(i in 1:length(splitted)){
    splitted[[i]][["Nucleotide Position"]] <- as.numeric(splitted[[i]][["Nucleotide Position"]])
    if("Conservation" %in% colnames(splitted[[i]])){
      splitted[[i]][["Conservation"]] <- as.numeric(gsub("%", "", splitted[[i]][["Conservation"]]))*0.01
    }
    if("AA Change♦)" %in% colnames(splitted[[i]])){
      splitted[[i]][["AA Change♦)"]] <- gsub(")", "", splitted[[i]][["AA Change♦)"]])
    }
    if("Mitomap Frequency‡in 61134 FL Seqs" %in% colnames(splitted[[i]])){
      splitted[[i]][["Mitomap Frequency‡in 61134 FL Seqs"]] <- as.numeric(gsub("%", "", splitted[[i]][["Mitomap Frequency‡in 61134 FL Seqs"]]))/61134
    }
    if("HelixFrequency‡‡‡in 195983 seqs" %in% colnames(splitted[[i]])){
      splitted[[i]][["HelixFrequency‡‡‡in 195983 seqs"]] <- as.numeric(gsub("%", "", splitted[[i]][["HelixFrequency‡‡‡in 195983 seqs"]]))/195983
    }
    if("gnomAD 3.1Frequency‡‡[count/total)" %in% colnames(splitted[[i]])){
      splitted[[i]][["gnomAD 3.1Frequency‡‡[count/total)"]] <- sapply(splitted[[i]][["gnomAD 3.1Frequency‡‡[count/total)"]], function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      })
    }
  }

  splitted

}



# list <- split_df(test)

# strsplit(list$`MITOMAP: mtDNA Coding Region Sequence Variants`$`Nucleotide†(AA Change♦)`, split = "\\(tRNA\\)")
#
# list$`MITOMAP: mtDNA Coding Region Sequence Variants` %>%
#   separate(`Nucleotide†(AA Change♦)`, into = c("Nucleotide", "type"), sep = "\\(", extra = "merge", fill = "right")
#
#
# ca <- httr::GET(url = "https://www.mitomap.org/mitomaster/conservation.cgi?loc=1390&locus=33&pos=1&res=S&all=0&frameshift=0&_=1738069039120")
# s <- content(ca, as = "text")
#
# ex <- fromJSON(s)
# ex



