# RNAmp
RNA Amplification in Mixed Primary Tumors

RNAmp calculates cancer cell specific transcriptional output using DNA and RNA sequencing data from primary tumors.

## Installation
```R
devtools::install_github('shlienlab/rnamp')
```

## Example
An example of running RNAmp is provided in `ExampleData/rnamp_example.R`

## Input Data
Example input data can be found in `ExampleData/inputs`. These include:
- `DNA.rtable` and `RNA.rtable` - DNA and RNA read counts for somatic and germline variant data produced by GATK’s ASEReadCounter tool
- `variant.txt` – A list of somatic and/or germline variants and their focal copy number
- `segments.txt`- Sample allele specific copy number data

## Outputs
Example input data can be found in `ExampleData/inputs`. These include:
- The `summary.txt` file: This file is the main output of RNAmp, and reports the estimated fold change in RNA output for the cancer cell fraction. In addition it also includes the number of somatic and LOH SNP variants used in the RNAmp calculation for the sample, and flags for whether the sample met minimum variant count thresholds for each variant type.
- The raw Rdata file (`rnamp.Rdata`): Containing all data produced from RNAmp
- The raw variant data (`rnamp.variants.txt`)

### RNAmp Example Run
```R
# Load the RNAmp library
library(rnamp)

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
```
