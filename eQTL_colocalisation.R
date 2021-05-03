
## eQTLGen Analyses 

## Read in eQTL Gen data 

eqtls = readRDS("eQTLgen_data.rds")

eqtls$BETA <- eqtls$Zscore/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))
eqtls$SE <- 1/sqrt((2*eqtls$MAF)*(1-eqtls$MAF)*(eqtls$NrSamples + eqtls$Zscore^2))


## Read in BayesR significant SNPs 
Sys.setlocale("LC_ALL", "C")
snps = read.csv("GWAS.csv")

## Subset to Cis 
snps = snps[snps$Type %in% "Cis", ]
## Subset to sentinel SNPs - top assoication per gene 

sentinel <- list()

for(i in unique(snps$Gene_of_Somamer)){ 
print(i)

tmp = snps[which(snps$Gene_of_Somamer %in% i),]

if(range(tmp$PIP)[1] == range(tmp$PIP)[2]) {
sentinel[[i]] <- tmp[which.max(abs(tmp$Beta)),] 

snps[which(snps$Gene_of_Somamer %in% i),]
} else { 

sentinel[[i]] <- tmp[which.max(abs(tmp$PIP)),] 
}
} 

snps1 = do.call("rbind", sentinel)

## Extract SNPs to query

list = unique(snps$Marker)

## Determine overlap between SNPs  

length(which(list %in% eqtls$SNP))
list = list[which(list %in% eqtls$SNP)]

## Subset SNP list to those with entries in eQTL dataset - will help resolve issues with gene names 

snps = snps[which(snps$Marker %in% list), ]


## Step 1 - Query the SNP 

query <- list()

for(i in list){ 
## pick out the gene related to the queried SNP 
tmp_snps <- snps[snps$Marker %in% i, ]
gene <- tmp_snps$Gene_of_Somamer
## find snps in eqtl dataset and store  
tmp_eqtl = eqtls[eqtls$SNP %in% i,]
## enter gene name - this will help find the ones with same name in both datasets 
tmp_eqtl$Gene_GWAS <- unique(gene)
query[[i]] <- tmp_eqtl 
print(i)
}

query <- do.call("rbind", query)

## Find Genes that overlap exactly 

genes_matched = query[query$Gene_GWAS == query$GeneSymbol,"Gene_GWAS"] 
unique(genes_matched)

[1] "C1QC"     "ANXA2"    "F7"       "EPHA1"    "CCL23"    "PCSK7"    "KYNU"
[8] "PDCD1LG2"


## Find Genes that do not overlap (exactly) 

genes_not_matched = query[-which(query$Gene_GWAS == query$GeneSymbol),"Gene_GWAS"] 
genes_not_matched = unique(genes_not_matched[-which(genes_not_matched %in% genes_matched)])

	 [1] "MATN3" "OLFM2" "CST2"  "APOA5" "FRZB"  "CNDP1" "FCRL4" "TIMP4" "GDF15"
	[10] "ASAH2"

## Find genes that might be similar in name but not exact match 
eqtls_not_matched = query[-which(query$Gene_GWAS == query$GeneSymbol),]
eqtls_not_matched = eqtls_not_matched[which(eqtls_not_matched$Gene_GWAS %in% genes_not_matched), ] 
eqtls_not_matched$GeneSymbol[order(unique(eqtls_not_matched$GeneSymbol))]


## Subset original SNP list to those which have matching gene names in 

snps1 = snps[which(snps$Gene_of_Somamer %in% c(genes_matched)),] 

## Now that we have eQTLs which definitely overlap in protein/mRNA transcript, we can set up colocalisation analysis 

## load data.table library - helps read in large GWAS files 
library(data.table)

## install/load coloc library for colocalisation analyses 
if(!require("remotes"))
   install.packages("remotes") # if necessary
library(remotes)
install_github("chr1swallace/coloc")

library(coloc)


## Genes that match up directly to eqtlGen consortium 

## Create output to accept coloc results 

output <- matrix(nrow = nrow(snps1), ncol = 9)
output <- as.data.frame(output)
names(output) <- c("SNP", "Somamer", "Gene", "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")

## Set up loop for colocalisation analyses 

for(i in 1:length(snps1$Marker)){ 

## Extract SNP 
tmp = snps1[i, ]
## Extract the SNPs coordinates 
position = as.numeric(tmp$hg19)
chr = as.numeric(tmp$Chr_of_Hit)

## Subset eqtl file to the chr of the queried SNP 
tmp_eqtl <- eqtls[which(eqtls$SNPChr %in% chr), ]
## Subset eqtl file to region +/- 200 kb from queried SNP 
tmp_eqtl <- tmp_eqtl[which(tmp_eqtl$SNPPos >= (position - 2e5)),]
tmp_eqtl <- tmp_eqtl[which(tmp_eqtl$SNPPos <= (position + 2e5)),]

tmp_eqtl = tmp_eqtl[tmp_eqtl$GeneSymbol %in% tmp$Gene_of_Somamer,]  
if(nrow(tmp_eqtl) == 0){ NULL }

else{ 

## Only those entries which relate to gene in question 

## Go to GWAS outputs and extract Somamer 

GWAS = fread(paste0("/Outputs/", tmp$Somamer, ".txt"))
GWAS = as.data.frame(GWAS)

print(tmp$Somamer)

## Subset GWAS file to chr + position 
GWAS <- GWAS[which(GWAS$CHR %in% chr), ]
## Subset eqtl file to region +/- 200 kb from queried SNP 
GWAS <- GWAS[which(GWAS$POS >= (position - 2e5)),]
GWAS <- GWAS[which(GWAS$POS <= (position + 2e5)),]
## Get MAF - coloc needs this 
GWAS$AF1 <- as.numeric(GWAS$AF1)
GWAS$MAF <- ifelse(GWAS$AF1 > 0.5, 1-GWAS$AF1, GWAS$AF1)



## Clean up datasets 
## Missing MAF 
if(length(which(is.na(tmp_eqtl$MAF))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(is.na(tmp_eqtl$MAF)),]
} 

if(length(which(is.na(GWAS$MAF))) >= 1) { 
  GWAS <- GWAS[-which(is.na(GWAS$MAF)),]
} 

## Duplicated SNPs
if(length(which(duplicated(tmp_eqtl$SNP))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(duplicated(tmp_eqtl$SNP)),]
} 

if(length(which(duplicated(GWAS$rsid))) >= 1) { 
  GWAS <- GWAS[-which(duplicated(GWAS$rsid)),]
} 

## Missing SNPs
if(length(which(is.na(tmp_eqtl$SNP))) >= 1) { 
  tmp_eqtl <- tmp_eqtl[-which(is.na(tmp_eqtl$SNP)),]
} 

if(length(which(is.na(GWAS$rsid))) >= 1) { 
  GWAS <- GWAS[-which(is.na(GWAS$rsid)),]
} 

if(all(is.na(as.numeric(tmp_eqtl$se)))){

dataset1 = list(snp = as.character(tmp_eqtl$SNP), N = as.numeric(tmp_eqtl$NrSamples),  pvalues = as.numeric(tmp_eqtl$Pvalue), MAF = as.numeric(tmp_eqtl$MAF), type = "quant", sdY =1 )
dataset2 = list(snp = as.character(GWAS$rsid), N = 1064, pvalues = as.numeric(GWAS$P), MAF = as.numeric(GWAS$MAF), type = "quant", sdY= 1)
coloc = coloc.abf(dataset2, dataset1)

## Create output file 

ind = i
print(ind)
output[ind, 1] <- tmp$Marker
output[ind, 2] <- tmp$Somamer
output[ind, 3] <- tmp$Gene_of_Somamer
output[ind, 4] <- coloc[[1]][[1]]
output[ind, 5] <- coloc[[1]][[2]]
output[ind, 6] <- coloc[[1]][[3]]
output[ind, 7] <- coloc[[1]][[4]]
output[ind, 8] <- coloc[[1]][[5]]
output[ind, 9] <- coloc[[1]][[6]]
} else { 
  
  ## Missing MAF or beta 
  if(length(which(is.na(as.numeric(tmp_eqtl$SE)))) >= 1) { 
    tmp_eqtl <- tmp_eqtl[-which(is.na(as.numeric(tmp_eqtl$SE))),]
  }
  
  }
} 
}


## Tidy up output 
output = output[rev(order(output$PP.H4)),] 

output[,c(5:9)] <- signif(output[,c(5:9)],2)

