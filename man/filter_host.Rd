% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_host.R
\name{filter_host}
\alias{filter_host}
\title{Align reads against one or more filter libraries and subsequently
remove mapped reads}
\usage{
filter_host(
  reads_bam,
  libs,
  lib_dir = NULL,
  output = paste(tools::file_path_sans_ext(reads_bam), "filtered", "bam", sep = "."),
  settings = align_details
)
}
\arguments{
\item{reads_bam}{The name of a merged, sorted .bam file that has previously
been aligned to a reference library. Likely, the output from running an
instance of \code{align_target()}.}

\item{libs}{A vector of character strings giving the basenames of the
Subread index files for alignment. These should be located in the same
directory and the current directory if lib_dir=NULL.}

\item{lib_dir}{path to the library index files (all libraries should be here)}

\item{output}{The desired name of the output .bam file. Default is
the basename of \code{reads_bam} + \code{.filtered.bam}.}

\item{settings}{A named \code{list} specifying alignment parameters for
the \code{Rsubread::align()} function, which is
called inside \code{align_target()}.
Elements should include type, nthreads, maxMismatches, nsubreads,
phredOffset, unique, and nBestLocations. Descriptions of these parameters
are available under \code{?Rsubread::align}. Defaults to the global
\code{align_details} object.}
}
\value{
The name of a filtered, sorted .bam file written to the user's
current working directory.
}
\description{
After a sample is aligned to a target library with \code{align_target()},
we may use \code{filter_host()} to remove unwelcome host contamination using
filter reference libraries. This function takes as input the name
of the .bam file produced via \code{align_target()}, and produces a
sorted .bam file with any reads that match the filter libraries removed.
This resulting .bam file may be used upstream for further analysis.
}
\examples{
#### Filter reads from bam file that align to any of the filter libraries

## Assuming a bam file has been created previously with align_target()

## Create object with path to the example filter library
refPath <- system.file("extdata","filter.fasta", package = "MetaScope")

## Copy the example filter library to the current directory
file.copy(from = refPath, to = file.path(".", "filter.fasta"))

## Make subread index of filter library
mk_subread_index('filter.fasta')

## Create object with path to the previously aligned bam file
readPath <- system.file("extdata", "subread_target.bam", package = "MetaScope")

## Filter bam file 
filter_host(readPath, libs = "filter")

}
