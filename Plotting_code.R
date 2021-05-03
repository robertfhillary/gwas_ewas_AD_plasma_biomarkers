## Install Prerequisite Packages 

if(!require(ggplot2)){
  install.packages("ggplot2")
}

library(ggplot2)


## GWAS Atlas 

gen = read.csv("GWAS.csv", check.names = F)
chr = read.csv("CHR_Length.csv")
chr_hit <- chr
names(chr)[2] <- "Max_Chromosome_Biomarker"
names(chr_hit)[2] <- "Max_Chromosome_Hit"

gen = merge(gen,chr_hit,by.x="Chromsome of Hit", by.y="CHR")
gen = merge(gen,chr,by.x="Chromsome of Hit", by.y="CHR")

gen$`Chromosome of SOMAmer` <- as.character(gen$`Chromosome of SOMAmer`)
gen[gen$`Chromosome of SOMAmer` %in% "X", "Chromosome of SOMAmer"] <- "23"
gen$`Chromosome of SOMAmer` = as.numeric(gen$`Chromosome of SOMAmer`)

gen$Relative_Hit_Position <- gen$`Position of Hit`/gen$Max_Chromosome_Hit
gen$Relative_Gene_Position <- gen$`Position of SOMAmer`/gen$Max_Chromosome_Biomarker

gen$Relative_Hit_Position <- gen$Relative_Hit_Position + gen$`Chromsome of Hit`
gen$Relative_Gene_Position <- gen$Relative_Gene_Position + gen$`Chromosome of SOMAmer`

pdf("GWAS_Atlas.pdf", width = 7.6, height = 6.3)
p = ggplot(gen, aes(gen$Relative_Hit_Position, gen$Relative_Gene_Position))
q = p + geom_jitter(aes(colour = as.factor(gen$Type)), size =2) + xlab("pQTL Position") + ylab("Protein Position")
q = q + scale_x_continuous(limits = c(1,24), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24), labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X", "")) + scale_y_continuous(limits = c(1,24), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24), labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X", "")) + theme_bw()
q = q + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic", size = 12)) + theme(panel.background = element_rect(colour = "black")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.1, color = "darkslategray4")
q + theme(axis.text=element_text(size=12))
dev.off()





## EWAS Atlas 

gen = read.csv("EWAS.csv", check.names = F)
chr = read.csv("CHR_Length.csv")
chr_hit <- chr
names(chr)[2] <- "Max_Chromosome_Biomarker"
names(chr_hit)[2] <- "Max_Chromosome_Hit"

gen = merge(gen,chr_hit,by.x="Chromosome_of_Hit", by.y="CHR")
gen = merge(gen,chr,by.x="Chromosome_of_Hit", by.y="CHR")

gen$Chromosome_of_Hit <- as.character(gen$Chromosome_of_Hit)
gen[gen$Chromosome_of_Hit %in% "X", "Chromosome_of_Hit"] <- "23"
gen$Chromosome_of_Somamer = as.character(gen$Chromosome_of_Somamer)
gen[gen$Chromosome_of_Somamer %in% "X", "Chromosome_of_Somamer"] <- "23"
gen$Chromosome_of_Hit = as.numeric(gen$Chromosome_of_Hit)
gen$Chromosome_of_Somamer = as.numeric(gen$Chromosome_of_Somamer)

gen$Relative_Hit_Position <- gen$Position_of_Hit/gen$Max_Chromosome_Hit
gen$Relative_Gene_Position <- gen$Position_of_Somamer/gen$Max_Chromosome_Biomarker

gen$Relative_Hit_Position <- gen$Relative_Hit_Position + gen$Chromosome_of_Hit
gen$Relative_Gene_Position <- gen$Relative_Gene_Position + gen$Chromosome_of_Somamer

pdf("EWAS_Atlas.pdf", width = 7.6, height = 6.3)
p = ggplot(gen, aes(gen$Relative_Hit_Position, gen$Relative_Gene_Position))
q = p + geom_jitter(aes(colour = as.factor(gen$Type)), size = 3) + xlab("CpG Position") + ylab("Protein Position")
q = q + scale_x_continuous(limits = c(1,24), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24), labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X", "")) + scale_y_continuous(limits = c(1,24), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24), labels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X", "")) + theme_bw()
q = q + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic", size = 12)) + theme(panel.background = element_rect(colour = "black")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.1, color = "darkslategray4")
q + theme(axis.text=element_text(size=12))
dev.off()


## Proportion of Variants and their Function 
gen$Function <- gsub(",.*", "", gen$Function)
gen$Function[gen$Function %in% "missense_variant"] <- "exon variant"
gen$Function[gen$Function %in% "synonymous_variant"] <- "exon variant"
gen$Function[gen$Function %in% "non_coding_transcript_exon_variant"] <- "exon variant"
gen$Function[gen$Function %in% "downstream_gene_variant"] <- "downstream variant"
gen$Function[gen$Function %in% "upstream_gene_variant"] <- "upstream variant"
gen$Function[gen$Function %in% "5_prime_UTR_variant"] <- "5'UTR variant"
gen$Function[gen$Function %in% "3_prime_UTR_variant"] <- "3'UTR variant"
gen$Function[gen$Function %in% "intergenic_variant"] <- "intergenic variant"
gen$Function[gen$Function %in% "intron_variant"] <- "intron variant"

a = as.data.frame(table(gen$Function))
a <- a[order(rev(a$Freq)),] 
a$Proportion_of_Variants <- a$Freq/56


names(a)[1] <- "Function"
a$Function <- factor(a$Function, levels = a$Function[order(a$Proportion_of_Variants)])

pdf("Proportion_of_variants_by_function.pdf", width = 7, height = 7)
bar = ggplot(data = a, aes(a$Function, a$Proportion_of_Variants,fill = a$Function)) +
  geom_bar(stat="identity") + ylab("Proportion of Variants")
bar1 = bar + theme(axis.title.x=element_blank()) + theme(axis.text.x = element_text(size = 12, colour = "grey20")) + theme(axis.title.y = element_text(size = 14)) + theme(axis.text.y = element_text(size = 12))
bar1 = bar1 + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + scale_y_continuous(limits = c(0,0.5), breaks = c(0,0.1,0.2,0.3,0.4,0.5))                                                                                 
bar1
dev.off()


## Replication Beta Correlation

betas <- read.csv("pQTLs_from_Literature_cleaned.csv")
pdf("Correlation_of_pQTLs_with_previous_pQTLs.pdf", width = 8, height = 6.5)
beta = ggplot(betas, aes(betas$beta, betas$Hillary_Beta))
beta1 = beta + geom_point(aes(color = factor(betas$Study)))

beta2 = beta1 + 
  geom_smooth(method = "lm", se = FALSE) + 
  xlab("Betas from Literature") + 
  ylab(expression(paste("Betas from Hillary ", italic("et al.")))) +
  ggtitle("Correlation of Betas for pQTLs") +
  theme(plot.title= element_text(hjust =0.5)) + 
  labs(colour="Study")
beta3 = beta2 + geom_abline(intercept = 0, slope = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
beta3 + annotate(x = 1.65, y = 0.55, label = expression(paste(italic("r"), " = 0.82")), size = 4.6, geom = "text") + theme(legend.title = element_text(hjust = 0.5))
dev.off()
