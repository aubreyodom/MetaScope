% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_refseq.R
\name{download_refseq}
\alias{download_refseq}
\title{Download RefSeq genome libraries}
\usage{
download_refseq(
  taxon,
  reference = TRUE,
  representative = FALSE,
  compress = TRUE,
  patho_out = FALSE
)
}
\arguments{
\item{taxon}{Select one taxon to download. The taxon name should be a
recognized NCBI scientific name, with no grammatical or capitalization
inconsistencies. All available taxonomies are visible by accessing the
\code{taxonomy_table} object included in the package.}

\item{reference}{Download only RefSeq reference genomes?
Defaults to \code{TRUE}.
Automatically set to \code{TRUE} if \code{representative = TRUE}.}

\item{representative}{Download only RefSeq representative genomes?
Defaults to \code{FALSE}. If \code{TRUE}, reference is automatically
set at \code{TRUE}.}

\item{compress}{Compress the output .fasta file? Defaults to \code{TRUE}.}

\item{patho_out}{Create duplicate outpute files compatible with PathoScope?
Defaults to \code{FALSE}.}
}
\value{
Returns a .fasta or .fasta.gz file of the desired RefSeq genomes.
This file is named after the kindom selectd and saved to the current
directory (e.g. 'bacteria.fasta.gz'). Currently, this function also returns
a .fasta file formatted for PathoScope as well
(e.g. bacteria.pathoscope.fasta.gz') if \code{path_out = TRUE}.
}
\description{
This function will automatically download RefSeq genome libraries in a
.fasta format from the specified taxon. The function will first
download the summary report at:
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/**kingdom**/assembly_summary.txt,
and then use this file to download the genome(s) and combine them in a single
compressed or uncompressed .fasta file.
}
\examples{
#### Download RefSeq genomes

## Download all RefSeq reference bacterial superkingdom genomes
download_refseq('bacteria', reference = TRUE, representative = FALSE)

## Download all RefSeq representative mononegavirales genomes
download_refseq('mononegavirales', representative = TRUE)

## Download all RefSeq morbillivirus genomes
download_refseq('morbillivirus', reference = FALSE)

## Download all RefSeq bacilli reference genomes, uncompressed
download_refseq('Bacilli', reference = TRUE,
                representative = FALSE, compress = FALSE)

## Download RefSeq Escherichia coli IAI1 strain 
download_refseq('Escherichia coli IAI1', reference = FALSE, compress = FALSE)

}
