% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msa_kit.R
\name{find_variants_AA}
\alias{find_variants_AA}
\title{Identify Sequence Variations in Multiple Sequence Alignment (MSA)}
\usage{
find_variants_AA(msa, target, ref = "NC_012920")
}
\arguments{
\item{msa}{A multiple sequence alignment object, a \code{DNAStringSet}
from the \code{Biostrings} package. Containing the target and the reference.}

\item{target}{String. Name of the target sequence in contained in the MSA.}

\item{ref}{String. Name of the reference sequence in the MSA.
Defaults to \code{"NC_012920"}.}
}
\value{
A data frame containing the identified variations between the target and reference
sequences. Variations are described using the following formats:
- Substitutions: \verb{<position><reference_base>><target_base>} (e.g., "123A>G").
- Insertions: \verb{<start_position>_<end_position>ins<bases>} (e.g., "123_124insAT").
- Deletions: \verb{<start_position>_<end_position>del<bases>} (e.g., "123_125delGTC").

\if{html}{\out{<div class="sourceCode">}}\preformatted{   Returns *NULL* if no variants are detected.
}\if{html}{\out{</div>}}
}
\description{
This function detects and returns variations between a target sequence and a reference sequence (typically rCRS)
within a given multiple sequence alignment (MSA). Variations are annotated following the
nomenclature described at \url{https://www.hgmd.cf.ac.uk/docs/mut_nom.html}. Note that different
nomenclature systems may appear in other studies.
}
\details{
Positions are adjusted to match the reference sequence, accounting for gaps introduced
during the alignment process. Ambiguous bases (e.g., "N") are ignored.
}
