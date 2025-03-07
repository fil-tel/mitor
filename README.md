
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mitor <img src="man/figures/logo.png" align="right" height="138" alt="logo"/>

<!-- badges: start -->
<!-- badges: end -->

## Overview

*mitor* is a kit that provides a set of tools for the handling of human
mitochondrial DNA (mtDNA) sequences and the identification and mapping
of non-synonymous variants in the 3D structure of the electron transport
chain (ETC) complexes in R.

Originally, *mitor* was developed as a pipeline for characterizing,
analyzing, and 3D visualizing the mutational landscape of archaic and
modern human mtDNA sequences. As we realized its potential for broader
applications and scenarios, we decided to develop it into a package, in
the hope that it will be useful to future users.

## Main features

Here is a summary of *mitor*’s main features. The R package allows you
to:

- **Fetch mitochondrial sequences** from the [NCBI nucleotide
  database](https://www.ncbi.nlm.nih.gov/nucleotide/) and download them
  on your local machine.
- **Haplogroup classification** of mtDNA sequences using
  [Haplogrep](https://haplogrep.i-med.ac.at/).
- **Handle a multiple sequence alignment (MSA) of human mtDNA
  sequences**, identify and extract the mitochondrial genes in the MSA,
  translate the coding sequences into amino acid sequences.
- **Identify variable positions** in a MSA of mtDNA sequences.
- **Mapping mitochondrial protein variable sites** for a given MSA onto
  the 3D structure of the ETC complexes and visualizing them using
  [PyMol](https://www.pymol.org/).
- **Query the [MITOMAP](https://www.mitomap.org/allelesearch.html)
  allele search database** to obtain information regarding variants in
  the human population at given positions.
- **Predict the structure of an amino acid sequence** using the
  [AlphaFold Server](https://alphafoldserver.com/). **Note:** this
  feature is currently under maintenance because some bugs were
  encountered.

## Installation

*mitor* can be installed directly from [GitHub](https://github.com/)
using the following commands:

``` r
# install.packages("devtools")
devtools::install_github("fil-tel/mitor")
```

**Note:** to ensure that *mitor* works correctly it is recommended to
manually install *Biostrings* from
[Bioconductor](https://www.bioconductor.org/) and *genbankr* from
[Github](https://github.com/gmbecker/genbankr) using the following
commands:

``` r
# Installation of Biostrings
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
```

``` r
# Installation of genbankr
# install.packages("devtools")
devtools::install_github("gmbecker/genbankr")
```

## External software

*mitor* has several functions that relies on two external software:
[PyMol](https://www.pymol.org/) and
[Haplogrep](https://haplogrep.i-med.ac.at/). These **does not** come
with the package and they have to be installed separetly by the user.
You can install them following the instructions listed in their
websites:

- [PyMol](https://pymol.org/dokuwiki/doku.php?id=installation) -
  required for 3D visualization of mtDNA variant sites.
- [Haplogrep](https://haplogrep.readthedocs.io/en/latest/installation/) -
  required for haplogroup classification.

**Note that** if the software’s binaries are not found in the R’s PATH
variable, you need to manually add them adding these two lines to your
*.Renviron* file.

``` r
## Add haplogrep binary's directory to PATH variable
PATH="${PATH}:/path/to/haplogrep/"
## Add pymol binary's directory to PATH variable
PATH="${PATH}:/path/to/pymol"
```

The *.Renviorn* file is usually found in *~/*.

## AlphaFold Server

One of *mitor*’s features is predicting the 3D structure of a given
amino acid sequence using the [AlphaFold
Server](https://alphafoldserver.com/). In order to “unlock” this feature
it is necessary to follow this step-by-step guide:

- **Create an account.** [AlphaFold
  Server](https://alphafoldserver.com/) requires you to create an
  account (you can log in using your Google account).
- **Main page**. Once you created your account and logged in you will
  have something like this:

<figure>
<img src="man/figures/alphafold_main.png"
alt="AlphaFold Server main page." />
<figcaption aria-hidden="true">AlphaFold Server main page.</figcaption>
</figure>

- **Copy request as cURL.** When submitting a job to the server, the
  server needs to know who is submitting the request. For this reason,
  AlphaFold Server is assigning you a token that identifies your account
  and some other information. As *mitor* needs to have access to this
  data to be able to submit jobs to the server, what you need to do now
  is to copy a call as cURL. Once on the main page, open the Web
  Developer Tools right-clicking anywhere in the page and selecting
  *Inspect(Q)*:

<figure>
<img src="man/figures/inspect.png" alt="Web Developer Tools" />
<figcaption aria-hidden="true">Web Developer Tools</figcaption>
</figure>

- Now a section like this will open, select the Network option:

<figure>
<img src="man/figures/network.png" alt="Network" />
<figcaption aria-hidden="true">Network</figcaption>
</figure>

- At this point, select one of the POST requests, right-click on it and
  select *Copy Value*, and then *Copy as cURL*:

<figure>
<img src="man/figures/curl.png" alt="Copy as cURL" />
<figcaption aria-hidden="true">Copy as cURL</figcaption>
</figure>

- Now you are almost done, go to the directory where you installed
  *mitor* (find.package(“mitor”) in R), and paste the cURL request in
  the *curl.txt* file that you find in *mitor/extdata/*.

- **You are now done!**

## Operating System Support

*mitor* has only been tested on Linux (Ubuntu version 22.04). Some
features will definitely not work on other Operating Systems.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(mitor)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.pngpressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
