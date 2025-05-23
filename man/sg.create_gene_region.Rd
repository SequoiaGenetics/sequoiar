% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_gene_region.R
\name{sg.create_gene_region}
\alias{sg.create_gene_region}
\title{Create a gene region list from RACER annotations}
\usage{
sg.create_gene_region(gene_name_chr, build = 38)
}
\arguments{
\item{gene_name_chr}{A character string. The gene name to look up (must match RACER annotation).}

\item{build}{An integer, either `38` (default) or `37`, specifying the genome build to use.}
}
\value{
A named list containing `gene_name`, `chr`, `lower` (start), and `upper` (end) if the gene is found uniquely;
otherwise returns `NULL`.
}
\description{
Given a gene name and genome build, this function retrieves the chromosome,
start, and end coordinates of the specified protein-coding gene using
annotations from the `RACER` package.
}
