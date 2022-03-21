#' Process RNAmp variant data
#'
#' \code{process.variants} will clean, annotate, and process variant data to
#' make it suitable for downstream processing.
#'
#' After importing data with \code{import.data}, this function will \enumerate{
#' \item CLEAN: Set missing counts to zero, remove MT variants and reorder by
#' chr:pos \item MERGE: Copy number and variant data \item COMPUTE: VAF in DNA
#' and RNA, mutation and SNP copy numbers and allelic ratio and amp \item FLAG:
#' Low coverage DNA and RNA sites, diploid and copy neutral regions, false
#' positive heterozygous SNPs, imprinted regions, and LOH alleles \item
#' ANNOTATE: Informative variants for transcriptional output estimation}
#' Processed variants are attached to the RNAmp data object as
#' 'variants.processed'.
#' @param rnamp.obj RNAmp data object from \code{import.data}
#' @param dna.depth Depth filter for DNA (Default: 8)
#' @param rna.depth Depth filter for RNA (Default: 30)
#' @export
process.variants <- function(rnamp.obj, dna.depth = 8, rna.depth = 30) {
    
    # make the var.df df available within our function
    var.df <- rnamp.obj$variants
    p <- rnamp.obj$purity
    
    cat("Setting missing counts to 0\n")
    var.df[is.na(var.df$refCount.dna), "refCount.dna"] <- 0
    var.df[is.na(var.df$altCount.dna), "altCount.dna"] <- 0
    var.df[is.na(var.df$totalCount.dna), "totalCount.dna"] <- 0
    
    var.df[is.na(var.df$refCount.rna), "refCount.rna"] <- 0
    var.df[is.na(var.df$altCount.rna), "altCount.rna"] <- 0
    var.df[is.na(var.df$totalCount.rna), "totalCount.rna"] <- 0
    
    # filter out mitochondrial chromosome
    var.df <- dplyr::filter(var.df, !grepl("M", CHROM))
    
    # remove chr prefix
    var.df$CHROM = factor(gsub(pattern = "chr", replacement = "", x = var.df$CHROM), levels = c("1", "2", 
        "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
        "21", "22", "X", "Y"))
    
    # reorder data by chromosome and position
    var.df <- var.df[order(var.df$CHROM, var.df$POS), ]
    
    # pre-filter variants with 0 coverage in DNA
    var.df <- dplyr::filter(var.df, totalCount.dna != 0)
    
    # should segment merging happen here? probably better than in the perl script in case people come from
    # different pipelines
    
    # calculate VAF in DNA and RNA
    cat("Calculing VAF in DNA and RNA\n")
    var.df$vaf.dna <- with(var.df, altCount.dna/totalCount.dna)
    var.df$vaf.rna <- with(var.df, altCount.rna/totalCount.rna)
    
    # save unadjusted VAF.dna and RNA values
    var.df$vaf.dna.unadj <- var.df$vaf.dna
    var.df$vaf.rna.unadj <- var.df$vaf.rna
    
    # flag poor coverage regions in DNA
    cat("Flag poorly covered variants in DNA (<", dna.depth, " read support)\n")
    var.df$dna_cov <- with(var.df, flag_coverage(totalCount.dna, dna.depth))
    
    # flag poor coverage regions in RNA
    cat("Flag poorly covered variants in RNA (<", rna.depth, " read support)\n")
    var.df$rna_cov <- with(var.df, flag_coverage(totalCount.rna, rna.depth))
    
    # flag for both DNA and RNA coverage
    var.df$cov <- with(var.df, factor(ifelse(dna_cov == "pass" & rna_cov == "pass", 0, 1), levels = c(0, 
        1), labels = c("pass", "fail")))
    
    # label sex chromosomes
    cat("Labelling sex chromosomes\n")
    var.df$sexChr <- with(var.df, factor(ifelse(grepl("X|Y", CHROM), 0, 1), levels = c(0, 1), labels = c("yes", 
        "no")))
    
    # label diploid genome regions
    cat("Labelling diploid genome regions\n")
    var.df$diploid <- with(var.df, factor(ifelse(MINCN == 1 & MAJCN == 1, 0, 1), levels = c(0, 1), labels = c("yes", 
        "no")))
    
    # label copy neutral genome regions
    cat("Labelling copy neutral genome regions\n")
    var.df$copy_neutral <- with(var.df, factor(ifelse(MINCN == MAJCN, 0, 1), levels = c(0, 1), labels = c("yes", 
        "no")))
    
    # flag false positive het snps
    cat("Flagging false positive heterozygous SNPs\n")
    var.df$het.fp <- with(var.df, factor(ifelse(ID == "GERMLINE" & (vaf.dna < 0.025 | vaf.dna > 0.975), 0, 
        1), levels = c(0, 1), labels = c("yes", "no")))
    
    # label putative imprinted SNPs -- unexpressed shows as NA?
    cat("Labelling imprinted germline SNPs\n")
    var.df$imprinted <- with(var.df, factor(ifelse((ID == "GERMLINE" & (vaf.rna < 0.025 | vaf.rna > 0.975)), 
        0, 1), levels = c(0, 1), labels = c("yes", "no")))
    
    # flag LOH var.df only for autosomes -- ignoring sites that have no genetic material
    cat("Flagging LOH var.df on autosomal chromosomes\n")
    var.df$loh <- with(var.df, factor(ifelse(MINCN == 0 & MAJCN != 0 & sexChr != "yes", 0, 1), levels = c(0, 
        1), labels = c("yes", "no")))
    
    # flag LOH SNPs
    cat("Flagging SNPs in LOH regions\n")
    var.df$loh.snp <- with(var.df, factor(ifelse(loh == "yes" & ID == "GERMLINE", 0, 1), levels = c(0, 1), 
        labels = c("yes", "no")))
    
    # measure SNP CN
    cat("Measuring tumour SNP CN\n")
    var.df$SNPCN <- with(var.df, ifelse(ID == "GERMLINE", ((vaf.dna * ((p * (MINCN + MAJCN)) + (2 * (1 - 
        p))) - (1 - p))/p), NA))
    var.df[which(var.df$SNPCN < 0), "SNPCN"] <- 0  #if SNPCN is negative make it zero
    
    # round SNP CN to the nearest modal value
    var.df$SNPCN.modal <- round(var.df$SNPCN)
    
    # Label mismatches in SNPCN
    var.df$SNPCN_mismatch <- with(var.df, factor(ifelse(SNPCN.modal == MINCN | SNPCN.modal == MAJCN, "SNPCN match", "SNPCN mismatch")))
    
    # Total copy number
    var.df$TUMCN = with(var.df, MAJCN + MINCN)
    
    # measure mutation copy number
    var.df$MUTCN = with(var.df, ifelse(ID == "SOMATIC", (vaf.dna/p) * ((p * TUMCN) + NORMCN * (1 - p)), NA))
    var.df$MUTCN.modal <- round(var.df$MUTCN)
    
    # identify LOH snp type
    var.df$loh.snp.type <- with(var.df, factor(ifelse(loh.snp == "yes" & (SNPCN.modal == MINCN), 0, ifelse(loh.snp == 
        "yes" & (SNPCN.modal == MAJCN), 1, ifelse(loh.snp == "yes", 2, ifelse(loh.snp == "no" & ID == "GERMLINE", 
        3, 4)))), levels = c(0, 1, 2, 3, 4), labels = c("LOH_alt", "LOH_ref", "LOH_unknown", "non_LOH", "somatic")))
    # summary(var.df$loh.snp.type)
    
    # classify variants
    var.df$variant.type <- with(var.df, factor(ifelse(loh.snp == "yes" & (SNPCN.modal == MINCN), 0, ifelse(loh.snp == 
        "yes" & (SNPCN.modal == MAJCN), 1, ifelse(loh.snp == "yes", 2, ifelse(loh.snp == "no" & ID == "GERMLINE" & 
        copy_neutral == "yes", 3, ifelse(loh.snp == "no" & ID == "GERMLINE" & copy_neutral == "no", 4, ifelse(ID == 
        "SOMATIC", 5, 6)))))), levels = c(0, 1, 2, 3, 4, 5), labels = c("SNP_LOH_alt", "SNP_LOH_ref", "SNP_LOH_unknown", 
        "SNP_copy_neutral", "SNP_copy_aberrant", "somatic")))
    # summary(var.df$variant.type)
    
    # flag informative variants for transcriptional output
    cat("Flagging informative variants\n")
    var.df$inform <- with(var.df, factor(ifelse((variant.type == "SNP_LOH_alt" | variant.type == "SNP_LOH_ref" | 
        variant.type == "somatic") & (CLASSIFICATION %in% c("nonsynonymous SNV", "synonymous SNV", "Missense_Mutation", 
        "Silent") & (het.fp == "no" & imprinted != "yes")), 0, 1), levels = c(0, 1), labels = c("yes", "no")))
    # with(var.df, table(variant.type, inform))
    
    # adjust vaf.dna and vaf,rna for LOH alt SNPs
    var.df[var.df$loh.snp.type == "LOH_alt", "vaf.dna"] <- 1 - var.df[var.df$loh.snp.type == "LOH_alt", "vaf.dna"]
    var.df[var.df$loh.snp.type == "LOH_alt", "vaf.rna"] <- 1 - var.df[var.df$loh.snp.type == "LOH_alt", "vaf.rna"]
    
    rnamp.obj$variants.processed <- var.df
    return(rnamp.obj)
}
