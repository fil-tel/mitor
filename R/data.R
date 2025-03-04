#' Cambridge Reference Sequence Gene Coordinates
#'
#' This dataset contains the coordinates of genes in the revised Cambridge reference sequence of the human mitochondrial genome.
#' It includes details such as gene locations, strand orientation, type, and nomenclature annotations.
#'
#' @format A data frame with 37 rows and 12 columns:
#' \describe{
#'   \item{\code{seqnames}}{Character. The origin of the sequence: "Homo sapiens".}
#'   \item{\code{start}}{Integer. The starting coordinate of the gene on the sequence.}
#'   \item{\code{end}}{Integer. The ending coordinate of the gene on the sequence.}
#'   \item{\code{width}}{Integer. The length of the gene in base pairs.}
#'   \item{\code{strand}}{Character. The strand orientation of the gene: \code{"+"} or \code{"-"}.}
#'   \item{\code{type}}{Character. \code{"gene"} for all.}
#'   \item{\code{gene}}{Character. The gene identifier (e.g., \code{TRNF}, \code{RNR1}).}
#'   \item{\code{nomenclature}}{Character. Nomenclature annotations including the official symbol, full name, and source of the information.}
#'   \item{\code{db_xref}}{Character. Database cross-references for the gene.}
#'   \item{\code{loctype}}{Character. Location type, \code{"normal"} for all.}
#'   \item{\code{gene_synonym}}{Character. Synonyms for the gene, if any.}
#'   \item{\code{gene_id}}{Character. The unique gene identifier.}
#' }
#'
#' @details
#' This dataset provides detailed annotations of the genes present in the Cambridge reference sequence, which represents
#' the human mitochondrial genome.
#'
#' @source [https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/)
#'
#'
"rCRS_genes_df"

#' Cambridge Reference Sequence - Human Mitochondrial Genome
#'
#' This object contains the full DNA sequence of the revised Cambridge Reference Sequence (rCRS), representing the human mitochondrial genome.
#'
#' @format A \code{DNAStringSet} object of length 1:
#' \describe{
#'   \item{\code{width}}{Integer. The total length of the DNA sequence, 16,569 base pairs.}
#'   \item{\code{seq}}{Character. The complete nucleotide sequence of the human mitochondrial genome.}
#'   \item{\code{names}}{Character. The accession number and identifier of the sequence, \code{"NC_012920"}.}
#' }
#'
#' @details
#' It serves as the reference for sequence alignment, annotation, and comparative studies in mitochondrial genetics.
#' For more details check this [website](https://www.mitomap.org/foswiki/bin/view/MITOMAP/HumanMitoSeq)
#'
#' @source [https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1/)
#'
#'
"rCRS"

#' gnomAD v3.1 Mitochondrial DNA Variants
#'
#' This dataframe contains the Mitochondrial DNA (mtDNA) variants for [gnomaAD](https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/).
#'
#' @format A \code{data frame} of dimension 10850X99:
#' \describe{
#'   \item{\code{CHROM}}{Character. The chromsome, "chrM" for all the entries.}
#'   \item{\code{POS}}{Integer. The position of the variant in the rCRS.}
#'   \item{\code{ID}}{Character. refSNP id for the variant. Not always available.}
#'    \item{\code{AF_hom}}{Double. Frequency of the variant in homoplasmic state in the database. }
#'    \item{\code{Allele}}{Character. Allele state of the alternative variant. }
#'        \item{\code{Consequence}}{Character. Type of variant: "intergenic_variant", "non_coding_transcript_exon_variant", "start_lost",
#'     "stop_gained&start_lost",     "synonymous_variant" ,
#'     "missense_variant",
#'     "frameshift_variant",
#'     "stop_gained",
#'     "frameshift_variant&stop_lost" ,
#'      "stop_lost",
#'      "stop_retained_variant",
#'      "coding_sequence_variant&3_prime_UTR_variant",
#'     "inframe_deletion",
#'      "incomplete_terminal_codon_variant&coding_sequence_variant". }
#'    \item{\code{Protein_position}}{Numeric. Position of the amino acid in the protein.}
#'    \item{\code{AA_ref}}{Character. Amino acid (one letter code) in the rCRS.}
#'    \item{\code{AA_alt}}{Character. Amino acid (one letter code) in the variant.}
#' }
#'
#' @details
#' This data frame contains information of the human mitochondrial variants as stored in the [gnomaAD](https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/).
#' Since the file on the database is a vcf file, this data frame allows the user to navigate through them in a more
#' user friendly way. Moreover, the data frame presents only the variants that passed the quality filter pipeline, as described in their [website](https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/).
#'
#' @source [gnomAD v3.1 Mitochondrial DNA Variants](https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/)
#'
#'
"gnomAD_df"


#' AA Coordinates of the Mitochondiral Encoded Subuints of Complex I
#'
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex I (MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6).
#'
#' @format A \code{data frame}.
#' \describe{
#'   \item{\code{Protein}}{Character. The protein name.}
#'   \item{\code{ChainID}}{Character. The ID of the chain in the pdb file.}
#'   \item{\code{Pos}}{Numeric. The position of the residue in the chain.}
#'    \item{\code{x}}{Numeric. x coordinate. }
#'        \item{\code{y}}{Numeric. y coordinate. }
#'            \item{\code{z}}{Numeric. z coordinate. }
#'    \item{\code{Location}}{Character. Location of the residue in the complex (intermembrane, membrane, matrix.}
#' }
#'
#' @details
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex I (MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6). The structure use to extract the coordinates it's available
#' [here](https://opm.phar.umich.edu/proteins/3532). The scripts used to extract the coordinates are available in the
#' directory \code{inst/scripts/coordinates/} available in the directory where the package \code{mitor} was installed.
#'
#'
"complexI_coord"

#' AA Coordinates of the Mitochondiral Encoded Subuints of Complex III
#'
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex III (MT-CYTB).
#'
#' @format A \code{data frame}.
#' \describe{
#'   \item{\code{Protein}}{Character. The protein name.}
#'   \item{\code{ChainID}}{Character. The ID of the chain in the pdb file.}
#'   \item{\code{Pos}}{Numeric. The position of the residue in the chain.}
#'    \item{\code{x}}{Numeric. x coordinate. }
#'        \item{\code{y}}{Numeric. y coordinate. }
#'            \item{\code{z}}{Numeric. z coordinate. }
#'    \item{\code{Location}}{Character. Location of the residue in the complex (intermembrane, membrane, matrix.}
#' }
#'
#' @details
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex III (MT-CYTB). The structure use to extract the coordinates it's available
#' [here](https://opm.phar.umich.edu/proteins/3531). The scripts used to extract the coordinates are available in the
#' directory \code{inst/scripts/coordinates/} available in the directory where the package \code{mitor} was installed.
#'
#'
"complexIII_coord"


#' AA Coordinates of the Mitochondiral Encoded Subuints of Complex IV
#'
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex IV (MT-COX1, MT-COX2, MT-COX3).
#'
#' @format A \code{data frame}.
#' \describe{
#'   \item{\code{Protein}}{Character. The protein name.}
#'   \item{\code{ChainID}}{Character. The ID of the chain in the pdb file.}
#'   \item{\code{Pos}}{Numeric. The position of the residue in the chain.}
#'    \item{\code{x}}{Numeric. x coordinate. }
#'        \item{\code{y}}{Numeric. y coordinate. }
#'            \item{\code{z}}{Numeric. z coordinate. }
#'    \item{\code{Location}}{Character. Location of the residue in the complex (intermembrane, membrane, matrix.}
#' }
#'
#' @details
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex IV (MT-COX1, MT-COX2, MT-COX3). The structure use to extract the coordinates it's available
#' [here](https://opm.phar.umich.edu/proteins/4308). The scripts used to extract the coordinates are available in the
#' directory \code{inst/scripts/coordinates/} available in the directory where the package \code{mitor} was installed.
#'
#'
"complexIV_coord"

#' AA Coordinates of the Mitochondiral Encoded Subuints of Complex V
#'
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex V (MT-ATP6, MT-ATP8).
#'
#' @format A \code{data frame}.
#' \describe{
#'   \item{\code{Protein}}{Character. The protein name.}
#'   \item{\code{ChainID}}{Character. The ID of the chain in the pdb file.}
#'   \item{\code{Pos}}{Numeric. The position of the residue in the chain.}
#'    \item{\code{x}}{Numeric. x coordinate. }
#'        \item{\code{y}}{Numeric. y coordinate. }
#'            \item{\code{z}}{Numeric. z coordinate. }
#'    \item{\code{Location}}{Character. Location of the residue in the complex (intermembrane, membrane, matrix.}
#' }
#'
#' @details
#' This data frame contains the coordinates of the carbon \eqn{\alpha} of the amino acids of the mitochondrial encoded subunits of
#' Complex V (MT-ATP6, MT-ATP8). The structure use to extract the coordinates it's available
#' [here](https://opm.phar.umich.edu/proteins/10340). The scripts used to extract the coordinates are available in the
#' directory \code{inst/scripts/coordinates/} available in the directory where the package \code{mitor} was installed.
#'
#'
"complexV_coord"
