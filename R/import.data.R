#' Import RNAmp raw data
#'
#' \code{import.data} Imports variant allelic read counts along with segmented
#' copy number data into an RNAmp object.
#'
#' This function will import data DNA and RNA allele count data and merge with copy number and annotation inormation or processing by RNAmp
#' @param DNA_ASE_Table,RNA_ASE_Table Allelic depth count tables produced for DNA
#'  and RNA sequencing from RNAmp.pl.
#' @param Variant_Table_file Variant table file produced from RNamp.pl script.
#' @param Segments_File Optional: Segments file.\cr Tab seperated file in the
#'  format: \emph{Chromosome, Start Position, End Position, Major CN, Minor CN}.
#' @param purity Optional: Sample purity
#' @param ploidy Optional: Sample ploidy
#' @param sample_name Optional: Name of the sample.
#' @export
import.data <- function(DNA_ASE_Table, RNA_ASE_Table, Variant_Table_File, Segments_File = NULL, purity = NULL, 
    ploidy = NULL, sample_name = NULL) {
    # define ASEReadCounter Column classes for ease
    colCls = c(position = "numeric", refCount = "numeric", altCount = "numeric", totalCount = "numeric", 
        lowMAPQDepth = "numeric", lowBaseQDepth = "numeric", rawDepth = "numeric", otherBases = "numeric", 
        improperPairs = "numeric")
    
    if (is.null(purity)) {
        cat("WARNING: Purity estimate not provided. Amp algorithm will not work.")
    }
    
    # read in allele counts
    cat("Reading DNA Allele Counts\n")
    dnaASE <- read.table(file = DNA_ASE_Table, header = TRUE, sep = "\t", colClasses = colCls)
    cat("Reading RNA Allele Counts\n")
    rnaASE <- read.table(file = RNA_ASE_Table, header = TRUE, sep = "\t", colClasses = colCls)
    
    # read in variant data
    cat("Reading Variant Table\n")
    variant_sites <- read.table(file = Variant_Table_File, header = TRUE, sep = "\t", quote = "\"")
    
    # check for existence of copy number segments
    segs = NULL
    if (!is.null(Segments_File) & !is.na(Segments_File) & Segments_File != "NA") {
        cat("Reading Allele Specific Copy Number Segments\n")
        segs <- read.table(file = Segments_File, header = TRUE, sep = "\t")
    } else {
        cat("No segments file providing. Assuming diploid sample.\n")
        segs = NULL
    }
    
    # split apart the ref and alt alleles for each site
    variant_sites$refAllele <- as.factor(substr(variant_sites$AA, 1, 1))
    variant_sites$altAllele <- as.factor(substr(variant_sites$AA, 2, 2))
    
    # merge allele counts with all initial sites
    variants.dna <- merge(x = variant_sites[, -4], y = dnaASE[, c(-4, -5)], by.x = c("CHROM", "POS", "ID"), 
        by.y = c("contig", "position", "variantID"), all = TRUE)
    
    variants.rna <- merge(x = variant_sites[, -4], y = rnaASE[, c(-4, -5)], by.x = c("CHROM", "POS", "ID"), 
        by.y = c("contig", "position", "variantID"), all = TRUE)
    
    # create unmelted table
    cat("Merging DNA and RNA counts\n")
    variants <- merge(x = variants.dna, y = variants.rna, all = FALSE, suffixes = c(".dna", ".rna"), by = c("CHROM", 
        "POS", "ID", "refAllele", "altAllele", "MAJCN", "MINCN", "NORMCN", "CLASSIFICATION"))
    
    # Add sample name to variants table
    variants <- cbind.data.frame(sample_name = rep(sample_name, length.out = nrow(variants)), variants)
    
    # create melted table (stack RNA on DNA)
    variants.dna$type <- as.factor("DNA")
    variants.rna$type <- as.factor("RNA")
    # variants.melted <- rbind(x = variants.dna, y = variants.rna, make.row.names = FALSE)
    
    return(list(sample_name = sample_name, variants = variants, segments = segs, purity = purity, ploidy = ploidy))
}
