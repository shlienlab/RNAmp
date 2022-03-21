# Load the RNAmp library
library(rnampR)

outdir = "output/"

# Data import
rnamp.obj <- import.data(DNA_ASE_Table = "input/DNA.rtable", RNA_ASE_Table = "input/RNA.rtable", Variant_Table_File = "input/variants.txt", Segments_File = "input/segments.txt", purity = 0.42, ploidy = 3.21, sample_name = "TCGA-05-4250-01A")

# Processing
rnamp.obj <- process.variants(rnamp.obj = rnamp.obj)

# Sample QC
rnamp.obj <- sample.qc(rnamp.obj = rnamp.obj)


# Calculate RNA output
rnamp.obj <- compute_amp(rnamp.obj)

# Save output data
cat("Saving data to", outdir, "\n")

# RNAMP Rdata file with all processed information
save(rnamp.obj, file = paste(outdir, "/", rnamp.obj$sample_name, ".rnamp.Rdata", sep = ""))

# Save sample summary output
write.table(rnamp.obj$sample.amp, file = paste0(outdir, "/", rnamp.obj$sample_name, ".summary.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)

# Save RNAmp variant list
write.table(rnamp.obj$sample.amp.variants, file = paste0(outdir, "/", rnamp.obj$sample_name, ".rnamp.variants.txt"), quote = FALSE, sep = "\t", row.names = F, col.names = T)
