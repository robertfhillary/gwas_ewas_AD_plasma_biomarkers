library(devtools)

library(TwoSampleMR)
library(MRInstruments)

snps <- read.csv("GWAS.csv")

list = unique(snps$Somamer)

## Set up loop of outcome data - Hong et al Biomarkers as an example

loop2 = list.files("CSF_biomarkers_Hong/", ".txt")
loop2 = loop2[-grep("IV", loop2)]

for(j in loop2){ 
output <- list()

for(i in list){ 

## Protein as exposure 

print(i)

## Open exposure - Protein
file = paste0("/MR_Inputs/",i, ".txt")


exposure_dat <- read_exposure_data( 
filename = file, 
sep = "\t", 
snp_col = 'rsid',
beta_col = 'BETA',
se_col = 'SE',
effect_allele_col = 'A1',
phenotype_col = '',
units_col = '', 
other_allele_col = 'A2',
eaf_col = 'AF1', 
samplesize_col = 'N', 
ncase_col = '',
ncontrol_col = '',
gene_col = '',
pval_col = 'P' 
)


## Clump Exposures 

exposure_dat <- clump_data(exposure_dat)

## Open Outcome  


file2 =  paste0("/CSF_biomarkers_Hong/", j)


outcome_dat <- read_outcome_data(
filename = file2,
sep = ",",
snp_col = 'SNP',
beta_col = 'BETA',
se_col = 'SE',
effect_allele_col = 'A1',
phenotype_col = '',
units_col = '',
other_allele_col = '',
eaf_col = 'FREQ1',
samplesize_col = '',
ncase_col = '',
ncontrol_col = '',
gene_col = '',
pval_col = 'P'
)


## Check to see if SNPs match 
if(length(which(exposure_dat$SNP %in% outcome_dat$SNP)) == 0) { 
NULL 
} else { 


## Harmonisation of Alleles 

dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)

  if(length(dat$mr_keep) == length(which(dat$mr_keep %in% "FALSE"))){ 
NULL

} else { 
## Run MR 

mr_results <- mr(dat) 
mr_results$Somamer <- i

output[[i]] <- mr_results

}
} 

}


output1 <- do.call("rbind", output)


write.csv(output1, paste0("Proteins_on_", j, "_Hong.csv"), row.names = F)

} 
