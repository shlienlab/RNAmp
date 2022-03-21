#' Compute sample amplification levels using the new meta-amp method
#' 
#' \code{compute_amp} calculates a sample's amplification level.
#' 
#' After importing data with \code{import.data}, and processing with
#' \code{process.variants}, this function will estimate a tumour sample's
#' transcriptional output level.
#' @param rnamp.obj RNAmp data object from \code{process.variants}
#' @export
#' @import dplyr
compute_amp <- function(rnamp.obj) {
  # pull out data from the rnamp obj
  var.df <- rnamp.obj$variants.processed
  p <- rnamp.obj$purity
  sample_name <- rnamp.obj$sample_name
  
  # calculate vaf difference
  var.df$vaf_diff <- var.df$vaf.rna - var.df$vaf.dna

  # Somatic part
  # Calculate mutation clonality
  clonal <- "NA"
  somatic.df <- var.df[var.df$ID == "SOMATIC" & var.df$TUMCN > 0,]
  
  if (nrow(somatic.df) > 0) {
    clon.dat <- suppressWarnings(t(with(somatic.df, (mapply(absolute.method, alt_reads = altCount.dna, tot_reads = totalCount.dna, purity = 1, total_cn = TUMCN)))))
    
    clon.dat <- data.frame(clon.dat, stringsAsFactors = FALSE)
    clon.dat[1:3] <- lapply(clon.dat[1:3], as.numeric)
    clon.dat[4] <- lapply(clon.dat[4], as.factor)
    
    clonal <- cbind(somatic.df, clon.dat)
    # Remerge onto the variants data frame
    var.df <- merge(var.df, clonal, all.x = TRUE) %>% refactor_df()
    
    var.df.f1.somatic <- dplyr::filter(var.df, sexChr == "no", ID == "SOMATIC", cov == "pass", CLASSIFICATION %in% c("nonsynonymous SNV", "synonymous SNV", "Missense_Mutation", "Silent"), clonal %in%  c("clonal"))
    
    # Meta somatic variant approach
    amp.somatic <- var.df.f1.somatic %>% summarise(var.count = n(), vaf.dna.adj.glob.mean = mean(vaf.dna), vaf.rna.glob.mean = mean(vaf.rna), mean.tumcn = mean(TUMCN), mean.mutcn = mean(MUTCN), mean.normcn = mean(NORMCN)) %>%
      mutate(mean.mutcn.recomputed = ((vaf.dna.adj.glob.mean / p) * ((p * mean.tumcn) + mean.normcn * (1 - p)))) %>% 
      mutate(amp.calc.mean = amp.somatic(vaf.dna = vaf.dna.adj.glob.mean, vaf.rna = vaf.rna.glob.mean, norm.cn = mean.normcn, tum.cn = mean.tumcn, purity = p)) %>%
      as.data.frame()
  } else {
    cat("No somatic mutations in this sample...\n")
    var.df.f1.somatic <- "NA"
    amp.somatic <- data.frame(matrix(ncol =11, nrow = 1))
    # Intialize empty dataframe to handle errors  
    colnames(amp.somatic) <- c("var.count", "vaf.dna.adj.glob.mean", 
                                         "vaf.rna.glob.mean", "mean.tumcn", "mean.mutcn", 
                                         "mean.normcn", "mean.mutcn.recomputed", "amp.calc.mean")
    
    amp.somatic$var.count = 0
    }
 
  
  # LOH SNPs
  variants.f1.loh.snp <- filter(var.df, sexChr == "no", cov == "pass", loh.snp == "yes", loh.snp.type %in% c("LOH_alt", "LOH_ref"), SNPCN_mismatch == "SNPCN match", CLASSIFICATION %in% c("Missense_Mutation", "Silent"), refCount.rna >=0, altCount.rna >=0)
  
  # Use MAJCN instead of snpcn b/c vaf.rna is flipped for ref lost
  amp.loh <- variants.f1.loh.snp %>% 
    summarise(var.count = n(), vaf.dna.adj.glob.mean = mean(vaf.dna), vaf.rna.glob.mean = mean(vaf.rna), mean.tumcn = mean(TUMCN), mean.snpcn = mean(MAJCN), mean.mutcn = mean(MUTCN), mean.normcn = mean(NORMCN)) %>%
    mutate(amp.calc.mean = amp.germline(vaf.rna = vaf.rna.glob.mean, tum.cn = mean.tumcn, purity = p, snp.cn = mean.snpcn)) %>% 
    as.data.frame()
  
  # Combine the loh and somatic amps
  amp.somatic$sample_name <- rnamp.obj$sample_name

  amp.loh$sample_name <- rnamp.obj$sample_name
  
  amp.comb <- merge(amp.somatic, amp.loh, by = c("sample_name"), suffixes = c(".somatic", ".loh"))
  
  amp.comb$prop.loh.var <- with(amp.comb, var.count.loh / (var.count.somatic + var.count.loh))
  amp.comb$prop.som.var <- with(amp.comb, var.count.somatic / (var.count.somatic + var.count.loh))
  
  # Get total meta var count
  amp.comb$var.count.total <- with(amp.comb, var.count.somatic + var.count.loh)
  
  # Convert all the negative or zero amp results to NAs
  amp.comb$amp.calc.mean.somatic <- ifelse(amp.comb$amp.calc.mean.somatic <=0, NA, amp.comb$amp.calc.mean.somatic)

  amp.comb$amp.calc.mean.loh <- ifelse(amp.comb$amp.calc.mean.loh <=0, NA, amp.comb$amp.calc.mean.loh)

  
  # Add a flag for samples passing certain # of variants
  amp.comb$loh.count.flag <- flag_coverage(depth = amp.comb$var.count.loh, threshold = 25, test = "greater eq")
  
  amp.comb$som.count.flag <- flag_coverage(depth = amp.comb$var.count.somatic, threshold = 25, test = "greater eq")
  
  amp.comb$count.flag <- with(amp.comb, factor(ifelse(loh.count.flag == "pass" & som.count.flag == "pass", "both", ifelse(loh.count.flag == "fail" & som.count.flag == "pass", "somatic", ifelse(loh.count.flag == "pass" & som.count.flag == "fail", "loh", "fail")))))
  

  # Stepwise weighted calculation depending on passing both, or either flag
  amp.comb$amp.calc.mean.weighted <- with(amp.comb, ifelse(count.flag == "both", (prop.loh.var * amp.calc.mean.loh) + (prop.som.var * amp.calc.mean.somatic), ifelse(count.flag == "loh", amp.calc.mean.loh, ifelse(count.flag == "somatic", amp.calc.mean.somatic, NA))))
  
  
  # Recover the samples that are 'both' but NA because one or the other fails
  amp.comb$amp.calc.mean.weighted <- with(amp.comb, ifelse(count.flag == "both" & is.na(amp.calc.mean.weighted) & !is.na(amp.calc.mean.loh), amp.calc.mean.loh, ifelse(count.flag == "both" & is.na(amp.calc.mean.weighted) & !is.na(amp.calc.mean.somatic), amp.calc.mean.somatic, amp.calc.mean.weighted)))
  

  # Create a 'Total_mean' variable for the overall RNA output value for the sample
  amp.comb$Total_mean <- amp.comb$amp.calc.mean.weighted
  

  
  # Simplify output
  sample.amp <- amp.comb[,c('sample_name', 'var.count.somatic', 'var.count.loh', 'var.count.total', 'som.count.flag', 'loh.count.flag', 'count.flag')]
  sample.amp$fold_change <- amp.comb$Total_mean
  
  # Add to the RNAmp object and return
  rnamp.obj$sample.amp <- sample.amp
  
  # Include detailed output
  rnamp.obj$sample.amp.dat <- amp.comb
  
  # Log all used variants in the calculations
  rnamp.obj$sample.amp.variants <- rbind(variants.f1.loh.snp, var.df.f1.somatic)

  return(rnamp.obj)
}
