## Colocalisation - methylation QTLs 

library(phenoscanner)
library(data.table)
library(coloc)

## Load in SNP information 

Sys.setlocale("LC_ALL", "C")
snps = read.csv("GWAS.csv")

## Load in mQTLs for lookups 

mqtls <- read.csv("mQTL_lookup.csv")



## Set up loop 

list = mqtls$SNP

## Create output 
output <- list()

for(i in 1:length(list)){ 
  
  tmp = mqtls[i,]
  
  
  
  print(i)
  
  ## Extract the SNP 
  snp = tmp$SNP
  ## Extract the Somamer 
  somamer = tmp$Somamer
  
  ## Extract the SNPs coordinates 
  position = snps[snps$Marker %in% snp, "hg19"]
  position = unique(position)
  chr = snps[snps$Marker %in% snp, "Chr_of_Hit"]
  
  
  
  ## Subset eqtl file to the chr of the queried SNP 
  tmp_mqtl <- phenoscanner(regionquery = paste0("chr", chr, ":", (position-2e5),"-", (position+2e5)), catalogue = "methQTL", proxies= "EUR")
  tmp_mqtl <- tmp_mqtl$results 
  tmp_mqtl <- tmp_mqtl[which(tmp_mqtl$tissue %in% c("Whole blood")), ]
  
  
  tmp_mqtl = tmp_mqtl[tmp_mqtl$marker %in% tmp$CpG,]  
  if(nrow(tmp_mqtl) == 0){ NULL } else{ 
    
    
    GWAS = fread(paste0("/GWAS_Outputs/", somamer, ".txt"))
    GWAS = as.data.frame(GWAS)
    
    
    
    ## Subset GWAS file to chr + position 
    GWAS <- GWAS[which(GWAS$CHR %in% chr), ]
    ## Subset eqtl file to region +/- 200 kb from queried SNP 
    GWAS <- GWAS[which(GWAS$POS >= (position - 2e5)),]
    GWAS <- GWAS[which(GWAS$POS <= (position + 2e5)),]
    
    ## Get MAF - coloc needs this 
    tmp_mqtl$eur <- as.numeric(tmp_mqtl$eur) 
    tmp_mqtl$MAF <- ifelse(tmp_mqtl$eur > 0.5, 1-tmp_mqtl$eur, tmp_mqtl$eur)
    GWAS$AF1 <- as.numeric(GWAS$AF1)
    GWAS$MAF <- ifelse(GWAS$AF1 > 0.5, 1-GWAS$AF1, GWAS$AF1)
    
    
    
    
    
    ## Missing SNPs
    
    if(length(which(is.na(tmp_mqtl$rsid))) >= 1) { 
      tmp_mqtl <- tmp_mqtl[-which(is.na(tmp_mqtl$rsid)),]
    } 
    
    
    if(length(which(is.na(GWAS$rsid))) >= 1) { 
      GWAS <- GWAS[-which(is.na(GWAS$rsid)),]
    } 
    
    ## Duplicated SNPs
    if(length(which(duplicated(tmp_mqtl$rsid))) >= 1) { 
      tmp_mqtl <- tmp_mqtl[-which(duplicated(tmp_mqtl$rsid)),]
    } 
    
    if(length(which(duplicated(GWAS$rsid))) >= 1) { 
      GWAS <- GWAS[-which(duplicated(GWAS$rsid)),]
    } 
    
    ## Clean up datasets 
    ## Missing MAF or beta 
    if(length(which(is.na(as.numeric(tmp_mqtl$MAF)))) >= 1) { 
      tmp_mqtl <- tmp_mqtl[-which(is.na(as.numeric(tmp_mqtl$MAF))),]
    } 
    
    if(all(is.na(as.numeric(tmp_mqtl$se)))){ 
      
      
      dataset1 = list(snp = as.character(tmp_mqtl$rsid), N = as.numeric(tmp_mqtl$n), pvalues = as.numeric(tmp_mqtl$p), MAF = tmp_mqtl$MAF, type = "quant")
      dataset2 = list(snp = as.character(GWAS$rsid), N = 1064, pvalues = as.numeric(GWAS$P), MAF = as.numeric(GWAS$MAF), type = "quant")
      
      coloc = coloc.abf(dataset1, dataset2)
      
      ## Create output file 
      
      out = coloc[[1]]
      out$SNP = tmp$SNP
      out$CpG = tmp$CpG 
      out$Somamer = tmp$Somamer
      out$Gene = tmp$Gene
      output[[i]] <- as.data.frame(out)
      
      
      
      
    } else { 
      
      if(length(which(is.na(as.numeric(tmp_mqtl$se)))) >= 1) { 
        tmp_mqtl <- tmp_mqtl[-which(is.na(as.numeric(tmp_mqtl$se))),]
      } 
      
      if(length(which(is.na(as.numeric(GWAS$SE)))) >= 1) { 
        GWAS <- GWAS[-which(is.na(as.numeric(GWAS$SE))),]
      } 
      
      
      
      dataset1 = list(snp = as.character(tmp_mqtl$rsid), N = as.numeric(tmp_mqtl$n), beta = as.numeric(tmp_mqtl$beta), varbeta = as.numeric(tmp_mqtl$se)^2, pvalues = as.numeric(tmp_mqtl$p), MAF = tmp_mqtl$MAF, type = "quant")
      dataset2 = list(snp = as.character(GWAS$rsid), N = 1064, pvalues = as.numeric(GWAS$P),  beta = as.numeric(GWAS$BETA), varbeta = as.numeric(GWAS$SE)^2, MAF = as.numeric(GWAS$MAF), type = "quant")
      
      coloc = coloc.abf(dataset1, dataset2)
      
      ## Create output file 
      out = coloc[[1]]
      out$SNP = tmp$SNP
      out$CpG = tmp$CpG 
      out$Somamer = tmp$Somamer
      out$Gene = tmp$Gene
      output[[i]] <- as.data.frame(out)
      
    }
  } 
}

output1 <- as.data.frame(do.call("rbind", output))

## Tidy up output 
output1 = output1[rev(order(output1$PP.H4)),] 

output1 = output1[,c(1,7,8,9,2:6)]

output1[,c(5:9)] <- signif(output1[,c(5:9)],2)
