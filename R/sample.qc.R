#' Get sample level QC statistics
#'
#' \code{sample.qc} will perform QC for a given sample.
#'
#' After importing data with \code{import.data}, and processing with
#' \code{process.variants}, this function will perform sample level QC assessing
#' mean and median depth of coverage per variant, number of variants per sample,
#' number of informative variants and their coverage. See the \strong{Metrics}
#' section for further details.
#'
#' A variant is considered as \strong{covered} if it has sufficient read depth
#' in both DNA and RNA as defined by \code{dna.depth} and \code{rna.depth}
#' thresholds in \code{process.variants}.
#'
#' A variant is considered \strong{informative} if it is useful for
#' transcriptional output measurements. These include somatic variants, and
#' heterozygous germline variants found in LOH regions.
#' @section Metrics:
#' \strong{n_var}: Total number of variants \cr
#' \strong{n_var_cov}: Total number of variants covered . \cr
#' \strong{prop_var_cov}: Proportion of variants covered \cr
#' \strong{n_var_inf}: Number of informative variants \cr
#' \strong{n_var_inf_cov}: Number of covered informative variants \cr
#' \strong{prop_inf_cov}: Proportion of informative variants covered \cr
#' \strong{prop_cov_inf}: Proportion of covered variants that are informative \cr
#' \strong{n_som_var}: Number of somatic variants \cr
#' \strong{n_som_var_cov}: Number of covered somatic variants \cr
#' \strong{prop_som_var_cov}: Proportion of somatic variants covered \cr
#' \strong{n_som_var_inf}: Number of informative somatic variants covered \cr
#' \strong{n_som_var_inf_cov}: Number of informative somatic variants covered \cr
#' \strong{prop_som_var_inf_cov}: Proportion of informative somatic variants covered \cr
#' \strong{n_germ}: Number of germline variants \cr
#' \strong{n_germ_cov}: Number of covered germline variants \cr
#' \strong{prop_germ_cov}: Proportion of germline variants covered \cr
#' \strong{n_germ_inf}: Number of informative germline variants \cr
#' \strong{n_germ_inf_cov}: Number of informative germline variants covered \cr
#' \strong{prop_germ_inf_cov}: Proportion of informative germline variants covered \cr
#' \strong{n_germ_loh_alt}: Number of germline LOH alt variants \cr
#' \strong{n_germ_loh_alt_cov}: Number of germline LOH alt variants covered \cr
#' \strong{prop_germ_loh_alt_cov}: Proportion of germline LOH alt variants covered \cr
#' \strong{n_germ_loh_ref}: Number of germline LOH ref variants  \cr
#' \strong{n_germ_loh_ref_cov}: Number of germline LOH ref variants covered \cr
#' \strong{prop_germ_loh_ref_cov}: Number of germline LOH ref variants covered \cr
#' \strong{alt_ref_ratio}: LOH alt to LOH ref germline variant ratio \cr
#' \strong{alt_ref_cov_ratio}: LOH alt to LOH ref germline variant coverage ratio \cr
#' \strong{n_germ_dip}: Number of germline diploid variants \cr
#' \strong{n_germ_dip_cov}: Number of germline diploid variants covered \cr
#' \strong{prop_germ_dip_cov}: Proportion of germline diploid variants covered \cr
#' \strong{mean_dna_cov}: Mean variant coverage in DNA \cr
#' \strong{med_dna_cov}: Median variant coverage in DNA \cr
#' \strong{prop_dna_cov}: Proprotion of variants covered in DNA \cr
#' \strong{n_var_dna_cov}: Number of variants covered in DNA \cr
#' \strong{mean_rna_cov}: Mean variant coverage in RNA \cr
#' \strong{med_rna_cov}: Median variant coverage in RNA \cr
#' \strong{prop_rna_cov}: Proprotion of variants covered in RNA \cr
#' \strong{n_var_rna_cov}: Number of variants covered in RNA \cr
#' @export
sample.qc <- function(rnamp.obj) {
    
    # pull out the processed variants
    var.df <- rnamp.obj$variants.processed
    sample_name <- rnamp.obj$sample_name
    
    ### Overall variant coverage statistics Total coverage
    n_var = nrow(var.df)
    n_var_cov = nrow(dplyr::filter(var.df, cov == "pass"))
    prop_var_cov = n_var_cov/n_var
    
    # Informative variant coverage
    n_var_inf = nrow(dplyr::filter(var.df, inform == "yes"))
    n_var_inf_cov = nrow(dplyr::filter(var.df, inform == "yes", cov == "pass"))
    prop_inf_cov = n_var_inf_cov/n_var_inf  #prop var.df that are cov.
    prop_cov_inf = n_var_inf_cov/n_var_cov  #prop. cov. var.df that are inform
    
    # Somatic variant coverage
    n_som_var = nrow(dplyr::filter(var.df, ID == "SOMATIC"))
    n_som_var_cov = nrow(dplyr::filter(var.df, ID == "SOMATIC", cov == "pass"))
    prop_som_var_cov = n_som_var_cov/n_som_var
    
    n_som_var_inf = nrow(dplyr::filter(var.df, ID == "SOMATIC", inform == "yes"))
    n_som_var_inf_cov = nrow(dplyr::filter(var.df, ID == "SOMATIC", inform == "yes", cov == "pass"))
    prop_som_var_inf_cov = n_som_var_inf_cov/n_som_var_inf
    
    
    # Germline variant coverage
    n_germ = nrow(dplyr::filter(var.df, ID == "GERMLINE"))
    n_germ_cov = nrow(dplyr::filter(var.df, ID == "GERMLINE", cov == "pass"))
    prop_germ_cov = n_germ_cov/n_germ
    
    # Informative germline variant coverage
    n_germ_inf = nrow(dplyr::filter(var.df, ID == "GERMLINE", inform == "yes"))
    n_germ_inf_cov = nrow(dplyr::filter(var.df, ID == "GERMLINE", inform == "yes", cov == "pass"))
    prop_germ_inf_cov = n_germ_inf_cov/n_germ_inf
    
    # SNP type coverage breakdown
    n_germ_loh_alt <- nrow(dplyr::filter(var.df, ID == "GERMLINE", variant.type == "SNP_LOH_alt"))
    n_germ_loh_alt_cov <- nrow(dplyr::filter(var.df, ID == "GERMLINE", variant.type == "SNP_LOH_alt", cov == 
        "pass"))
    prop_germ_loh_alt_cov <- n_germ_loh_alt_cov/n_germ_loh_alt
    
    n_germ_loh_ref <- nrow(dplyr::filter(var.df, ID == "GERMLINE", variant.type == "SNP_LOH_ref"))
    n_germ_loh_ref_cov <- nrow(dplyr::filter(var.df, ID == "GERMLINE", variant.type == "SNP_LOH_ref", cov == 
        "pass"))
    prop_germ_loh_ref_cov <- n_germ_loh_ref_cov/n_germ_loh_ref
    
    alt_ref_ratio <- n_germ_loh_alt/n_germ_loh_ref
    alt_ref_cov_ratio <- n_germ_loh_alt_cov/n_germ_loh_ref_cov
    
    n_germ_dip <- nrow(dplyr::filter(var.df, ID == "GERMLINE", ID == "GERMLINE", diploid == "yes"))
    n_germ_dip_cov <- nrow(dplyr::filter(var.df, ID == "GERMLINE", ID == "GERMLINE", diploid == "yes", cov == 
        "pass"))
    prop_germ_dip_cov <- n_germ_dip_cov/n_germ_dip
    
    # DNA coverage statistics
    mean_dna_cov = mean(dplyr::filter(var.df, totalCount.dna > 0)$totalCount.dna)  #mean coverage
    med_dna_cov = median(dplyr::filter(var.df, totalCount.dna > 0)$totalCount.dna)  #median coverage
    prop_dna_cov = nrow(dplyr::filter(var.df, dna_cov == "pass"))/nrow(var.df)  #proportion covered
    n_var_dna_cov = nrow(dplyr::filter(var.df, dna_cov == "pass"))
    
    # RNA coverage statistics
    mean_rna_cov = mean(dplyr::filter(var.df, totalCount.rna > 0)$totalCount.rna)
    med_rna_cov = median(dplyr::filter(var.df, totalCount.rna > 0)$totalCount.rna)
    prop_rna_cov = nrow(dplyr::filter(var.df, rna_cov == "pass"))/nrow(var.df)
    n_var_rna_cov = nrow(dplyr::filter(var.df, rna_cov == "pass"))
    
    # Create statistics summary
    qc <- cbind.data.frame(rep(sample_name, length.out = 1), round(cbind(n_var, n_var_cov, prop_var_cov, 
        n_var_inf, n_var_inf_cov, prop_inf_cov, prop_cov_inf, n_som_var, n_som_var_cov, prop_som_var_cov, 
        n_som_var_inf, n_som_var_inf_cov, prop_som_var_inf_cov, n_germ, n_germ_cov, prop_germ_cov, n_germ_inf, 
        n_germ_inf_cov, prop_germ_inf_cov, n_germ_loh_alt, n_germ_loh_alt_cov, prop_germ_loh_alt_cov, n_germ_loh_ref, 
        n_germ_loh_ref_cov, prop_germ_loh_ref_cov, alt_ref_ratio, alt_ref_cov_ratio, n_germ_dip, n_germ_dip_cov, 
        prop_germ_dip_cov, mean_dna_cov, med_dna_cov, prop_dna_cov, n_var_dna_cov, mean_rna_cov, med_rna_cov, 
        prop_rna_cov, n_var_rna_cov), 3))
    
    rnamp.obj$qc <- qc
    return(rnamp.obj)
}
