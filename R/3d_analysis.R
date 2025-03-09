

eucl_distance <- function(x1, x2) {
  sqrt(sum((x1 - x2) ^ 2))
}

k_t <- function(variants, complex = NULL) {
  if (is.null(complex))
    stop("Specify the complex of interest. I, III, IV, V.")
  df <- switch (
    complex,
    I = complexI_coord,
    III = complexIII_coord,
    IV = complexIV_coord,
    V = complexV_coord
  )

  coords <- merge(variants, df, by = c("Protein", "Pos"))



}

#' Get 3D Coordinates of Protein Variants
#'
#' @param prot_list
#' @param t
#'
#' @return Data frame containing the 3D coordinates of the positions of non-synonymous variants
#' @export
#'
get_coords <- function(prot_list, t=NULL){
  prot_df <- data.frame(Protein = character(),
                        Protein_position = character(),
                        Pattern = character(),
                        stringsAsFactors = FALSE)
  if(is.null(t)){
    ls <- lapply(prot_list, function(x) colnames(find_var_pos(x)))
  }
  else{
    ls <- lapply(prot_list, function(x) find_var_pos(x, t = t))
  }

  for (name in names(ls)) {
    for (pos in ls[[name]]) {
      prot_df <- rbind(prot_df, cbind(data.frame("Protein" = paste0("MT-", name), "Pos"=pos)))
    }
  }
  merge(prot_df, rbind(complexI_coord, complexIII_coord, complexIV_coord, complexV_coord), by=c("Protein", "Pos"))
}

# get_coords(protein_list, t=0.99)
