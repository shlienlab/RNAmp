% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equations.R
\name{amp.germline}
\alias{amp.germline}
\title{Calculate transcriptional amplification value in LOH germline variants.}
\usage{
amp.germline(vaf.rna, purity, tum.cn, snp.cn)
}
\arguments{
\item{vaf.rna}{Variant allele fraction in RNA}

\item{purity}{Tumour purity. Numeric [0-1].}

\item{tum.cn}{Local total DNA copy number in the tumour cell fraction.}

\item{snp.cn}{SNP copy number. Should be equal to 0 or the major CN.}
}
\description{
\code{amp.germline} will use RNA variant fractions, along with local tumour
copy number and purity information to calculate the amplification value for a
given germline variant in an LOH region.
}
\details{
This algorithm uses tumour purity, local copy number, SNP copy number, and
the dna and rna fractions of a given germline variant to infer the level of
transcriptional amplification taking place. This will only work for germline
variants in regions of LOH.
}
