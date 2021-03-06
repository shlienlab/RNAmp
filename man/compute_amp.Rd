% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_amp.R
\name{compute_amp}
\alias{compute_amp}
\title{Compute sample amplification levels using the new meta-amp method}
\usage{
compute_amp(rnamp.obj)
}
\arguments{
\item{rnamp.obj}{RNAmp data object from \code{process.variants}}
}
\description{
\code{compute_amp} calculates a sample's amplification level.
}
\details{
After importing data with \code{import.data}, and processing with
\code{process.variants}, this function will estimate a tumour sample's
transcriptional output level.
}
