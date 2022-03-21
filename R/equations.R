#' Calculate allelic ratio
#'
#' \code{allelic.ratio} Given DNA and RNA variant fractions, this function will
#' calculate the allelic ratio for a given variant.
#'
#' The allelic ratio calculation is similar to allele specific expression, but
#' can also be applied to somatic mutations and copy number changed regions for
#' both germline and somatic variants. When the RNA fraction has increased, the
#' allelic ratio is given by the equation: \deqn{(vaf.rna - vaf.dna) / (1 -
#' vaf.dna)} and when the RNA fraction has decreased by: \deqn{(vaf.rna -
#' vaf.dna) / vaf.dna}
#'
#' This formula provides an indirect measure of transcriptional output, but has
#' the advantage of not requiring copy number, purity, or ploidy information.
#' @param vaf.dna Variant allele fraction in DNA
#' @param vaf.rna Variant allele fraction in RNA
#' @export
allelic.ratio <- function(vaf.dna, vaf.rna) {
    vaf.diff = vaf.rna - vaf.dna
    ifelse(vaf.rna == vaf.dna, 0, ifelse(vaf.rna > vaf.dna, ((vaf.diff/(1 - vaf.dna))), ifelse(vaf.rna < 
        vaf.dna, (vaf.diff/vaf.dna), NA)))
}

#' Calculate transcriptional amplification value in somatic variants.
#'
#' \code{amp.somatic} will use DNA and RNA variant fractions, along with local
#' tumour copy number and purity information to calculate the amplification
#' value for a given somatic variant.
#'
#' This algorithm uses tumour purity, local copy number, and the dna and rna
#' fractions of a given somatic variant to infer the level of transcriptional
#' amplification taking place.
#' @inheritParams allelic.ratio
#' @param norm.cn Local total DNA copy number in the normal cell fraction.
#'   (Default: 2).
#' @param tum.cn Local total DNA copy number in the tumour cell fraction.
#' @param purity Tumour purity. Numeric [0-1].
#' @export
amp.somatic <- function(vaf.dna, vaf.rna, norm.cn = 2, tum.cn, purity) {
    amp = ((vaf.rna * norm.cn * (1 - purity))/(((vaf.dna * norm.cn * (1 - purity)) - (purity * tum.cn * (vaf.rna - 
        vaf.dna)))))
    return(amp)
}

# alternate formula not used calc_amp_alt <- function(tum.cn, mut.cn, norm.cn, vaf.rna, purity) { amp <-
# (-norm.cn * vaf.rna + norm.cn * purity * vaf.rna)/(purity * (-mut.cn + tum.cn * vaf.rna)) return(amp) }

#' Calculate transcriptional amplification value in LOH germline variants.
#'
#' \code{amp.germline} will use RNA variant fractions, along with local tumour
#' copy number and purity information to calculate the amplification value for a
#' given germline variant in an LOH region.
#'
#' This algorithm uses tumour purity, local copy number, SNP copy number, and
#' the dna and rna fractions of a given germline variant to infer the level of
#' transcriptional amplification taking place. This will only work for germline
#' variants in regions of LOH.
#' @inheritParams amp.somatic
#' @param snp.cn SNP copy number. Should be equal to 0 or the major CN.
#' @export
#' 
amp.germline <- function(vaf.rna, purity, tum.cn, snp.cn) {
    amp <- ( (1 - purity) + 2 * vaf.rna * (purity - 1)  ) / ( purity * ( tum.cn * vaf.rna - snp.cn ) )
    return(amp)
}
amp.germline.orig <- function(vaf.rna, norm.cn = 2, purity, tum.cn, snp.cn) {
    # Need to copy and adjust assuming ewual gene dosage when LOH occurs in the tumour and recompute over the
    # data
    amp <- ((2 * purity * vaf.rna - 2 * vaf.rna - purity + 1)/(purity * (tum.cn * vaf.rna - snp.cn)))
    return(amp)
}



#' Calculate transcriptional amplification value.
#'
#' \code{calc.amp} is a wrapper for \code{amp.somatic} and \code{amp.germline}
#' used to calculate amplification values across all informative variants.
#' @inheritParams amp.somatic
#' @inheritParams amp.germline
#' @param ID Indicate whether the variant is somatic, \emph{0} or germline,
#'   \emph{1}.
#' @export
calc.amp <- function(vaf.dna, vaf.rna, norm.cn = 2, tum.cn, snp.cn = NULL, purity, ID) {
    amp <- ifelse(ID == "0" | ID == "SOMATIC", amp.somatic(vaf.dna, vaf.rna, norm.cn, tum.cn, purity), amp.germline.orig(vaf.rna = vaf.rna, norm.cn = norm.cn, purity = purity, tum.cn = tum.cn, snp.cn = snp.cn))
    return(amp)
}

#' Calculate tumour RNA content.
#'
#' \code{rna.content} estimates the amount of tumour derived RNA in a given
#' primary tumour sample based on its purity and amplification level.
#'
#' @inheritParams amp.somatic
#' @param amp Mean amplification level for a tumour.
#' @export
rna.content <- function(amp, purity) {
    tr = purity/(((1 - purity)/amp) + purity)
    return(tr)
}

#' Refactor data frames
#'
#' \code{refactor_df} Refactorizes data frames
#'
#' @export
refactor_df <- function(df) {
    fct.idx <- lapply(X = df[], FUN = is.factor) %>% unlist %>% which
    df[fct.idx] <- lapply(X = df[fct.idx], factor)
    return(df)
}



### Function: Std Dev of Posterior Beta
sd_of_posterior <- function(m, n, N, Y) {
    a <- Y + (n * m)
    b <- N - Y + (n * (1 - m))
    sigma_posterior <- sqrt((a * b)/(((a + b)^2) * (a + b + 1)))
    return(sigma_posterior)
}
