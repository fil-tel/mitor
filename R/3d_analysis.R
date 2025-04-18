
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
#' @param prot_list List. Protein list generated by translating the output of *extract_cds*.
#' @param t Double. Frequency threshold, select positions in the MSA which its most conserved amino acid/nucleotide has a frequency of at most *t*.
#' @param target Character. Optional, if specified, the function returns the variants between the reference and the target sequence.
#' @param ref Character. Reference sequence.
#'
#' @return Data frame containing the 3D coordinates of the positions of non-synonymous variants. It can be
#' used by the function *map_variants* to visualize the sites in 3D
#' @export
#'
get_coords <-
  function(prot_list,
           t = NULL,
           target = NULL,
           ref = "NC_012920") {
    prot_df <- data.frame(
      Protein = character(),
      Protein_position = character(),
      stringsAsFactors = FALSE
    )
    if (is.null(t) &
        is.null(target)) {
      ls <- lapply(prot_list, function(x)
        colnames(find_var_pos(x)))
      for (name in names(ls)) {
        for (pos in ls[[name]]) {
          prot_df <-
            rbind(prot_df, cbind(data.frame(
              "Protein" = paste0("MT-", name), "Pos" = pos
            )))
        }
      }
    }
    else if (is.null(target)) {
      ls <- lapply(prot_list, function(x)
        find_var_pos(x, t = t))
      for (name in names(ls)) {
        for (pos in ls[[name]]) {
          prot_df <-
            rbind(prot_df, cbind(data.frame(
              "Protein" = paste0("MT-", name), "Pos" = pos
            )))
        }
      }
    }
    else{
      ls <- lapply(prot_list, function(x)
        find_variants_AA(x, target = target, ref = ref))
      for (name in names(ls)) {
        df <- ls[[name]]
        if (!is.null(df)) {
          # print(nrow(df))
          prot_df <-
            rbind(prot_df, cbind(data.frame("Protein" = rep(
              paste0("MT-", name), nrow(df)
            )), df))
        }
      }
      names(prot_df)[names(prot_df) == "Position"] <- "Pos"
      prot_df$Ref <- seqinr::aaa(prot_df$Ref)
      prot_df$Alt <- seqinr::aaa(prot_df$Alt)
      prot_df$Gr_score <- grantham::grantham_distance(prot_df$Ref, prot_df$Alt)$d
    }

    # print(prot_df)

    if (nrow(prot_df) == 0)
      stop("There are no variants in your MSA.")

    merge(
      prot_df,
      rbind(
        complexI_coord,
        complexIII_coord,
        complexIV_coord,
        complexV_coord
      ),
      by = c("Protein", "Pos")
    )
  }

#' Map Non-Synonymous Variants To Protein Structure
#'
#' `map_variants()` maps the location of the variants on the 3D structure of the complexes using PyMol.
#' The variants location are highlighted in red.
#'
#' @param variants Data frame. Data frame containing the location of the variants, generated by the function *get_coords*
#' @param complex String. Roman number corresponding to the electron transport chain complexes of "I", "III", "IV", "V"
#'
#' @export
#'
#'
map_variants <- function(variants, complex = NULL) {
  if (Sys.which("pymol") == "")
    stop(
      "The pymol binary couldn't be found. \n Did you add it to your PATH variable? \n
  Visit https://github.com/fil-tel/mitor for a step-by-step guide.                                                                   Visit githublink for a step-by-step guide."
    )
  if (is.null(complex))
    stop("You have to specify which Electron Transport Chain complex you want to visualize!\n")
  df <- switch (
    complex,
    I = merge(complexI_coord, variants),
    III = merge(complexIII_coord, variants),
    IV = merge(complexIV_coord, variants),
    V = merge(complexV_coord, variants)
  )

  struc <- switch (
    complex,
    I = "5xtd.pdb",
    III = "5xte.pdb",
    IV = "5z62.pdb",
    V = "8h9t.pdb"
  )
  input_file <- tempfile(fileext = ".txt")
  write.csv(df, input_file, row.names = FALSE, quote = FALSE)
  script <-
    system.file("extdata", "pymol4mitor.py", package = "mitor")
  struc_path <- system.file("extdata/pdb", struc, package = "mitor")

  future::plan(future::multisession)
  future::future(system(
    paste0("pymol ", script, " ", input_file, " ", struc_path, " ", complex)
  ))
}
