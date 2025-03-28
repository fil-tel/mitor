% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy.R
\name{position_entropy_dna}
\alias{position_entropy_dna}
\title{Find Shannon Entropy for each position of a MSA}
\usage{
position_entropy_dna(msa)
}
\arguments{
\item{msa}{A multiple sequence alignment object, a \code{DNAStringSet}
from the \code{Biostrings} package.}
}
\value{
A vector with the entropy value for each position of the MSA.
}
\description{
Find Shannon Entropy for each position of a MSA
}
