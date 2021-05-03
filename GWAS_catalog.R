## Read in GWAS Hits to be queried 

gwas = read.csv("GWAS.csv")

## Set up list of SNPs to be queried 
list = as.array(as.character(gwas$Marker))

## Load catalog 
library(gwasrapidd)

## query the SNPs 
results = get_associations(variant_id = list)

## As it a GenomicRanges object, we need to process it to a dataframe

## this releases the S4 GenomicRanges object
unclass(results)

## This takes out the studies
associations = slot(results, "associations")
## This takes out the SNPs we queried
ids = slot(results, "risk_alleles")

## Check the order of associations are same in both 
table(associations$association_id == ids$association_id)

## Include SNP information in association file 
associations$variant_id <- ids$variant_id
associations$risk_allele <- ids$risk_allele
associations$genome_wide <- ids$genome_wide


## Convert association IDs to study and trait info

traits= association_to_trait(associations$association_id, verbose = T)
traits = as.data.frame(traits)

together = merge(traits, associations, by= "association_id")

## Get studies from EFO IDs

studies = get_traits(efo_id = traits$efo_id)
studies1 = slot(studies, "traits")
studies2 = as.data.frame(studies1)

final = merge(studies2, together, by = "efo_id")
final$uri <- NULL
final = final[which(final$pvalue < 5e-8),]
final2 = final[-which(duplicated(final$association_id)),]

ad = c("Alzheimer's disease", "t-tau measurement", "p-tau measurement", "memory performance", "amyloid-beta measurement", "cerebrospinal fluid biomarker measurement", "cerebral amyloid deposition measurement", "cerebral amyloid deposition measurement", "p-tau:beta-amyloid 1-42 ratio measurement","cognitive impairment measurement",  "family history of Alzheimer???s disease","Mental deterioration")
final3 = final2[final2$trait %in% ad,]

final2$AD_relevant <- "-"
final2[which(final2$association_id %in% final3$association_id),"AD_relevant"] <- "Yes"

final_gwas <- final2[,c("variant_id", "risk_allele", "trait", "or_per_copy_number", "beta_number", "beta_direction", "pvalue", "AD_relevant")]
names(final_gwas) <- c("SNP", "Risk Allele", "Trait", "OR per allele", "Beta per allele", "Beta Direction", "P Value", "AD Relevant Association")

gwas1 <- gwas[,c("Marker","Gene.of.Hit", "SOMAmer", "Gene.of.SOMAmer")]
names(gwas1) <- c("SNP", "Gene of Hit", "SOMAmer", "Gene of SOMAmer")

final_gwas <- merge(gwas1, final_gwas, by = "SNP")

final_gwas <- final_gwas[rev(order(final_gwas$`AD Relevant Association`)), ]

write.csv(final_gwas, "GWAS_Catalog.csv", row.names = F)
