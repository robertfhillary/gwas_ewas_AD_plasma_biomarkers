## Moloc script 

library(moloc)
library(phenoscanner)
Sys.setlocale("LC_ALL", "C")
snps = read.csv("GWAS.csv")
## Read in eQTL data 

eqtls = readRDS("eQTLgen_data.rds")

eqtls$BETA <- eqtls$Zscore/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))
eqtls$SE <- 1/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))


## Read in GWAS hits with colocalised eQTLs 

eqtl_coloc <- read.csv("eQTLs_lookup.csv") 

## Read in mQTLs 

mqtls  <- read.csv("mqtls_lookup.csv") 

## Find overlap - these are SNPs which have eQTL association (supported by coloc evidence) and mQTLs 

eqtl_coloc <- eqtl_coloc[which(eqtl_coloc$SNP %in% mqtls$SNP) ,]

## Trim to those with coloc PP.H4 > 0.75

eqtl_coloc <- eqtl_coloc[eqtl_coloc$PP.H4 > 0.75, ]

nrow(eqtl_coloc) ## 2 


## Moloc loop

output <- matrix(nrow = nrow(eqtl_coloc), ncol = 7)
output <- as.data.frame(output)
names(output) <- c("SNP", "Somamer", "Gene", "bestSNP", "nsnps", "PP_eqtl", "PP_mqtl")

for(i in 1:length(eqtl_coloc$SNP)){ 

## get snp
snp = eqtl_coloc$SNP[i]
## get chr of snp 
chr= snps[which(snps$Marker %in% snp) , "Chr_of_Hit"]
## get pos of snp 
position = snps[which(snps$Marker %in% snp) , "Position_of_Hit"]
## get somamer  
gene = "PCSK7"
somamer = mqtls[mqtls$SNP %in% snp, "SOMAmer"]

 ##### EQTL DATASET 
## Subset eqtl file to the chr of the queried SNP 
tmp_eqtl <- eqtls[which(eqtls$SNPChr %in% chr), ]
## Subset eqtl file to region +/- 200 kb from queried SNP 
tmp_eqtl <- tmp_eqtl[which(tmp_eqtl$SNPPos >= (position - 2e5)),]
tmp_eqtl <- tmp_eqtl[which(tmp_eqtl$SNPPos <= (position + 2e5)),]
tmp_eqtl <- tmp_eqtl[tmp_eqtl$GeneSymbol %in% gene, ]


### PQTL DATASET - GWAS 
## Go to GWAS outputs and extract Somamer 

GWAS = fread(paste0("/Outputs/", somamer, ".txt"))
GWAS = as.data.frame(GWAS)

print(somamer)

## Subset GWAS file to chr + position 
GWAS <- GWAS[which(GWAS$CHR %in% chr), ]
## Subset eqtl file to region +/- 200 kb from queried SNP 
GWAS <- GWAS[which(GWAS$POS >= (position - 2e5)),]
GWAS <- GWAS[which(GWAS$POS <= (position + 2e5)),]
## Get MAF - coloc needs this 
GWAS$AF1 <- as.numeric(GWAS$AF1)
GWAS$MAF <- ifelse(GWAS$AF1 > 0.5, 1-GWAS$AF1, GWAS$AF1)


### MQTL DATASET - PHENOSCANNER 

cpg = mqtls[which(mqtls$SOMAmer %in% somamer),"CpG"]
GWAS_mqtl = fread(paste0("/Outputs/", cpg, ".txt"))
GWAS_mqtl = as.data.frame(GWAS_mqtl)

print(cpg)

## Subset GWAS file to chr + position 
GWAS_mqtl <- GWAS_mqtl[which(GWAS_mqtl$CHR %in% chr), ]
## Subset eqtl file to region +/- 200 kb from queried SNP 
GWAS_mqtl <- GWAS_mqtl[which(GWAS_mqtl$POS >= (position - 2e5)),]
GWAS_mqtl <- GWAS_mqtl[which(GWAS_mqtl$POS <= (position + 2e5)),]
## Get MAF - coloc needs this 
GWAS_mqtl$AF1 <- as.numeric(GWAS_mqtl$AF1)
GWAS_mqtl$MAF <- ifelse(GWAS_mqtl$AF1 > 0.5, 1-GWAS_mqtl$AF1, GWAS_mqtl$AF1)


## Clean up datasets 
## Missing MAF
if(length(which(is.na(tmp_eqtl$MAF))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(is.na(tmp_eqtl$MAF)),]
} 


if(length(which(is.na(GWAS$MAF))) >= 1) { 
  GWAS <- GWAS[-which(is.na(GWAS$MAF)),]
} 

if(length(which(is.na(GWAS_mqtl$MAF))) >= 1) { 
  GWAS_mqtl <- GWAS_mqtl[-which(is.na(GWAS_mqtl$MAF)),]
} 

## Duplicated SNPs
if(length(which(duplicated(tmp_eqtl$SNP))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(duplicated(tmp_eqtl$SNP)),]
} 

if(length(which(duplicated(GWAS$rsid))) >= 1) { 
  GWAS <- GWAS[-which(duplicated(GWAS$rsid)),]
} 

if(length(which(duplicated(GWAS_mqtl$rsid))) >= 1) { 
  GWAS_mqtl <- GWAS_mqtl[-which(duplicated(GWAS_mqtl$rsid)),]
} 

## Missing beta

if(length(which(is.na(tmp_eqtl$BETA))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(is.na(tmp_eqtl$BETA)),]
} 


if(length(which(is.na(GWAS$BETA))) >= 1) { 
  GWAS <- GWAS[-which(is.na(GWAS$BETA)),]
} 

if(length(which(is.na(as.numeric(GWAS_mqtl$BETA)))) >= 1) { 
  GWAS_mqtl <- GWAS_mqtl[-which(is.na(as.numeric(GWAS_mqtl$BETA))),]
} 


## Prepare final datasets 

tmp_eqtl <- tmp_eqtl[,c("SNP", "BETA", "SE", "MAF", "NrSamples")]
names(tmp_eqtl) <- c("SNP", "BETA", "SE", "MAF", "N") 
tmp_eqtl$BETA <- as.numeric(tmp_eqtl$BETA) 
tmp_eqtl$SE <- as.numeric(tmp_eqtl$SE) 
tmp_eqtl$MAF <- as.numeric(tmp_eqtl$MAF) 
tmp_eqtl$N <- as.numeric(tmp_eqtl$N) 


GWAS <- GWAS[,c("rsid", "BETA", "SE", "MAF", "N")]
names(GWAS) <- c("SNP", "BETA", "SE", "MAF", "N") 
GWAS$BETA <- as.numeric(GWAS$BETA) 
GWAS$SE <- as.numeric(GWAS$SE) 
GWAS$MAF <- as.numeric(GWAS$MAF) 
GWAS$N <- as.numeric(GWAS$N) 

GWAS_mqtl <- GWAS_mqtl[,c("rsid", "BETA", "SE", "MAF", "N")]
names(GWAS_mqtl) <- c("SNP", "BETA", "SE", "MAF", "N") 
GWAS_mqtl$BETA <- as.numeric(GWAS_mqtl$BETA) 
GWAS_mqtl$SE <- as.numeric(GWAS_mqtl$SE) 
GWAS_mqtl$MAF <- as.numeric(GWAS_mqtl$MAF) 
GWAS_mqtl$N <- as.numeric(GWAS_mqtl$N) 

list_of_data <- list(GWAS, tmp_eqtl, GWAS_mqtl)

moloc <- moloc_test(list_of_data, prior_var = c(0.01, 0.1, 0.5), priors = c(1e-04, 1e-06,
  1e-07), save.SNP.info = FALSE)

## Store output 
ind = i
output[ind, 1] <- snp
output[ind, 2] <- somamer
output[ind, 3] <- somamer_gene
output[ind, 4] <- moloc$best_snp[[2]][1]
output[ind, 5] <- moloc$nsnps
output[ind, 6] <- moloc$priors_lkl_ppa[[4]][[10]] + moloc$priors_lkl_ppa[[4]][[11]] + moloc$priors_lkl_ppa[[4]][[14]]
output[ind, 7] <- moloc$priors_lkl_ppa[[4]][[8]] + moloc$priors_lkl_ppa[[4]][[12]] + moloc$priors_lkl_ppa[[4]][[14]]


} 

## Tidy and save output 
output[,c(6,7)] <- signif(output[,c(6,7)],2)


write.csv("moloc_results.csv",  row.names = F)
