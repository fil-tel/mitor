.onAttach <- function(libname, pkgname) {
  ## check for Pymol and Haplogrep libraries
  ## a warning if it's not present
  cmds  <- c("pymol", "haplogrep3")
  path_check_pymol <- Sys.which(cmds[1]) != ""
  path_check_haplo <- Sys.which(cmds[2]) != ""
  if (!path_check_pymol) {
    packageStartupMessage(
      "The PyMol binary required by some of the mitor's functions are not in your $PATH variable.\n",
      "This will not prevent you from using most of the mitor's functions, but you won't be able to map the variant sites on the 3D structure of complexes.\n",
      "\nOne way to solve this is to modify the $PATH variable in your ~/.Renviron\n",
      "file so that it points to the directory containing the PyMol binary.\n\n",
      "Current contents of the $PATH:\n",
      Sys.getenv("PATH", "\n")
    )
  }
  if (!path_check_haplo) {
    packageStartupMessage(
      "The haplogrep binary required by some of the mitor's functions are not in your $PATH variable.\n",
      "This will not prevent you from using most of the mitor's functions, but you won't be able to classify haplogroups using haplogrep.\n",
      "\nOne way to solve this is to modify the $PATH variable in your ~/.Renviron\n",
      "file so that it points to the directory containing the haplogrep binary.\n\n",
      "Current contents of the $PATH:\n",
      Sys.getenv("PATH", "\n")
    )
  }
}
