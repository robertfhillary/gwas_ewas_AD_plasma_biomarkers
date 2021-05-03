

library(data.table)

## Open ewascatalog 

ewascatalog <- fread("EWAS_Catalog_03-07-2019.txt") 
ewascatalog <- as.data.frame(ewascatalog)

## Read in cpgs 

cpgs <- read.csv("EWAS.csv")

## Trim Gene info for CpGs 

cpgs$Gene_of_Hit <- gsub(";.*", "", cpgs$Gene_of_Hit)

## Set up loop 

list = cpgs$Marker 

output <- list()

## Query loop 

for(i in list){ 

print(i)

output[[i]] = ewascatalog[which(ewascatalog$CpG %in% i), c("CpG", "Trait", "Gene", "P")]

}

output <- do.call("rbind", output)

## Subset to Genome-wide significant 

output1 <- output[output$P < 3.6e-8, ]


## Tidy up CpG and Trait Information 

list = unique(output1$CpG)

output <- matrix(nrow = length(list), ncol = 6)
output <- as.data.frame(output)
names(output) <- c("CpG", "CpG Gene", "No. of Somamers associated with CpG", "Somamers", "Somamer Gene(s)", "Trait in EWAS Catalog at P < 3.6e-8")


for(i in list){ 

ind = which(list %in% i) 

output[ind,1] <- i
output[ind,2] <- unique(cpgs[which(cpgs$Marker %in% i), "Gene_of_Hit"])
output[ind, 3] <- length(unique(cpgs[which(cpgs$Marker %in% i),"Somamer"]))
output[ind, 4] <- paste(cpgs[which(cpgs$Marker %in% i),"Somamer"], collapse = ", ")
output[ind,5] <- paste(unique(cpgs[which(cpgs$Marker %in% i),"Gene_of_Somamer"]), collapse = ", ")
output[ind, 6] <- paste(unique(output1[which(output1$CpG %in% i), "Trait"])[order(unique(output1[which(output1$CpG %in% i), "Trait"]))], collapse = ", ")

}

## Save file 

write.csv(output, "EWAS_Catalog.csv", row.names = F)
