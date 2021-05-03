## Combined BayesR Variance Output

var = read.csv("Combined_Variance_bayes.csv")
var = var[order(var$Mean_Combined_Variance_Explained), ]
protein <- rep(var$SOMAmer, each=3)
condition <- rep(c("Epigenetics" , "Genetics" , "Combined"), 276)
data <- data.frame(protein,condition)
data$value <- 0
var[,2] <- var[,2]*100
var[,5] <- var[,5]*100
var[,8] <- var[,8]*100
data$value[data$condition %in% "Epigenetics"] <- var$Mean_Epigenetic_Variance_Explained
data$value[data$condition %in% "Genetics"] <- var$Mean_Genetic_Variance_Explained
data$value[data$condition %in% "Combined"] <- var$Mean_Combined_Variance_Explained
data$condition <- factor(data$condition, levels = c("Epigenetics", "Genetics", "Combined"))


# Grouped

pdf("Additional_File4.pdf", width = 9, height = 6.5)

for(i in 0:45){ 
data1 <- data[(1+(18*i)):(18+(18*i)),]
x = ggplot(data1, aes(fill=condition, y=value, x=protein)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  +
  scale_fill_manual(values = wes_palette("Cavalcanti1", n = 3))+ xlab("SOMAmer") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "Data Type") 
print(x)
} 
dev.off()
