library(kinship2)
library(coxme)
library(data.table)


ped = read.csv("clinical/pedigree.csv")

ped$father <- as.numeric(ped$father)
ped$mother <- as.numeric(ped$mother)
ped$father[ped$father==0] <- NA
ped$mother[ped$mother==0] <- NA
table(ped$sex)
ped$sex <- as.numeric(as.factor(ped$sex))
ped$sex[ped$sex==2] <- 0
ped$sex <- ped$sex+1
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 


# Function to Extract Lmekin Results	
extract_coxme_table <- function (mod){
beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}


## Read in Methylation File 

meth = fread("meth_file.txt")
meth = as.data.frame(meth)
## Read phenotypes file 

pheno = read.csv("pheno.csv")

names(pheno) <- gsub("\\.", "_", names(pheno))


ids = meth$IID 
pheno = pheno[match(ids, pheno$id), ]

## Read in Covariates 

covariates <- read.table("Covariates.cov", header = T)
covariates$FID <- NULL 
covariates$IID <- gsub(".*GS_", "", covariates$IID)
names(covariates) <- c("IID", "CT1", "CT2", "CT3", "CT4", "CT5")
covariates <- covariates[match(ids, covariates$IID), ]

## Read in EWAS hits for lmekin loop 

ewas <- read.csv("EWAS.csv")

list1 <- ewas$Marker
list2 <- ewas$Somamer


output <- matrix(nrow=nrow(ewas),ncol=6)
output <- as.data.frame(output)
for(i in 1:nrow(ewas)){ 
  mod = lmekin(scale(meth[,list1[i]]) ~ scale(pheno[,list2[i]]) + covariates$CT1 + covariates$CT2 + covariates$CT3 + covariates$CT4 + covariates$CT5 + (1|meth$IID), varlist = kin_model*2)
 mod1 <- as.data.frame(extract_coxme_table(mod))

 output[i,1] <- as.character(list1[i])
    output[i,2] <- as.character(list2[i])
    output[i,3:6] <- round(mod1[2,1:4],2)

print(i)
}
output$V7 <- ewas$Beta
 output$V6[output$V6 == 0] <- "<2e-16"
names(output) <- c("CpG", "SOMAmer", "Beta", "SE", "Z", "P", "BayesR_Beta")
output <- output[,c(1,2,3,7,4,5,6)]

output$BayesR_Beta <- signif(output$BayesR_Beta,2)

write.csv(output, "Lmekin_sensitivity_EWAS.csv", row.names = F)
