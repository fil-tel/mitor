

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
