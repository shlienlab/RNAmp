# Subclonal functions
# Measure clones vs subclones as described by ABSOLUTE
# As in absolute (not exactly but the CCF CI method):
#' Absolute simplified clonal method
#'
#' \code{absolute.method} computes the cancer cell fraction, confidence intervals, and clonal or subclonal status of a given mutation
#' 
#' Using total and alternative read counts, sample purity, locus total copy number, this method computes the posterior probability distribution over given cancer cell fraction values (from 0.01 to 1). Mutations were defined as clonal if the 95% confidence interval overlapped 1, and subclonal otherwise
#'
#' @param alt_reads Number of alternative reads for the mutation
#' @param tot_reads Number of total reads for the mutation
#' @param purity Estimated sample purity
#' @param total_cn Total integer copy number at the mutation's locus
#' @param ccf_grid Grid of ccf values to compute the distribtion (default: 0.01 - 1)
#' @export
absolute.method <- function(alt_reads, tot_reads, purity, total_cn, ccf_grid = seq(0, 1, by=0.01)) {
  # For testing
  # ccf_grid = seq(0, 1, by=0.01)
  # alt_reads = 0
  # tot_reads = 60
  # purity = 0.5
  # total_cn = 2
  ###
  
  vaf.dna = alt_reads / tot_reads
  
  f = (purity * ccf_grid)/(2 * (1 - purity) + purity * total_cn)
  ll_grid = dbinom(alt_reads, tot_reads, f, log=TRUE)
  # not sure what these two lines do -- converts into probabilty grid -- posterior probability?
  log_add = ll_grid[which.max(ll_grid)] + log(sum(exp(ll_grid - ll_grid[which.max(ll_grid)])))
  pr_grid = exp(ll_grid  - log_add)
  ###
  ccf_hat = ccf_grid[which.max(pr_grid)]
  ecdf = cumsum(pr_grid)
  ccf_ci95 = approx(x=ecdf, y=ccf_grid, xout=c(0.025, 0.975))$y
  
  clonal = ifelse(ccf_ci95[2] > 0.95, "clonal", "subclonal")
  
  res = c(ccf_hat, ccf_ci95, clonal)
  names(res) = c("ccf", "ccf_95CI_low", "ccf_95CI_high", "clonal")
  return((res))
}

#' Measure subclones
#'
#' \code{measure.subclones} will take processed rnamp object and apply  a simplified method from ABSOLUTE to measure clonal and subclonal mutations
#'

#' @param rnamp.obj RNAmp data object from \code{import.data}
#' @export
measure.subclones <- function(rnamp.obj) {
  somatic.df <- rnamp.obj$variants.processed[rnamp.obj$variants.processed$ID == "SOMATIC" & rnamp.obj$variants.processed$dna_cov == "pass" & rnamp.obj$variants.processed$TUMCN > 0,]
  purity = rnamp.obj$purity
  
  clonal <- "NA"
  
  if (nrow(somatic.df) > 0) {
    clon.dat <- suppressWarnings(t(with(somatic.df, (mapply(absolute.method, alt_reads = altCount.dna, tot_reads = totalCount.dna, purity = purity, total_cn = TUMCN)))))
    
    clon.dat <- data.frame(clon.dat, stringsAsFactors = FALSE)
    clon.dat[1:3] <- lapply(clon.dat[1:3], as.numeric)
    clon.dat[4] <- lapply(clon.dat[4], as.factor)
    
    clonal <- cbind(somatic.df, clon.dat)
    rnamp.obj$variants.amp <- merge(rnamp.obj$variants.amp, clonal, all.x = TRUE)
  }
  
  # also add clonal information to amp used variant df
  rnamp.obj$clonal <- clonal
  return(rnamp.obj)
}
