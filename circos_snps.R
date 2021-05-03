library(circlize)


snps = read.csv("GWAS.csv")

snps = snps[snps$Type %in% "Trans", ]

## Make first bed file - positions of pQTLs

bed1 <- matrix(nrow = nrow(snps), ncol = 4)
bed1 <- as.data.frame(bed1)
names(bed1) <- c("chr", "start", "end", "value1")

## Set up loop to populate bed1 file 

list = snps$Marker

for(i in 1:length(list)){ 

## get chromosome 
chr = snps[i, "Chr_of_Hit"]
chr = paste0("chr", chr)

## get start/end 
start.end <- snps[i, "Position_of_Hit"]

## value = beta
value <- snps[i, "PIP"]

## populate the table

bed1[i, "chr"] <- chr
bed1[i, "start"] <- start.end
bed1[i, "end"] <- start.end+5e5
bed1[i, "value1"] <- value


}


## Create bed2 file - target of interest (genes of interest)

bed2 <- matrix(nrow = nrow(snps), ncol = 3)
bed2 <- as.data.frame(bed2)
names(bed2) <- c("chr", "start", "end")

for(i in 1:length(list)){ 

## get chromosome 
chr = snps[i, "Chr_of_Somamer"]
chr = paste0("chr", chr)

## get start/end 
start.end <- snps[i, "Position_of_Somamer"]


## populate the table

bed2[i, "chr"] <- chr
bed2[i, "start"] <- start.end
bed2[i, "end"] <- start.end+5e5

}


## Create annotation file 

names1 <- bed1
names2 <- bed2

## give value1 column in bed1 (now names1) gene names of hit + denote its a SNP 

names1$value1 <- snps$Gene_of_Hit
names1$value1 <- paste(names1$value1, "(pQTL)")

## give value1 column in bed2 (now names2) gene names of Somamer + denote its a Somamer

names2$value1 <- 0
names2$value1 <- snps$Gene_of_Somamer
names2$value1 <- paste(names2$value1, "(SOMAmer)")

## create bed file - annotation 

bed = rbind(names1, names2)



pdf("file_here.pdf", width = 7.5, height = 7.5)

## initalize plot 

circos.initializeWithIdeogram()

circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), track.height = 0.1, bg.border = NA)
i_track = get.cell.meta.data("track.index") # remember this empty track, we'll come back



circos.genomicTrackPlotRegion(bed, ylim = c(0, 1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicText(region, value, y = 1, labels.column = 1,
                                                   facing = "clockwise", adj = c(1, 0.5),
                                                   posTransform = posTransform.text, cex = 0.5, padding =2)
                              }, track.height = 0.1, bg.border = NA)
tr_track = get.cell.meta.data("track.index") # position transformation track
# because `circos.genomicPosTransformLines` is implemented by
# `circos.trackPlotRegion`, it accepts `track.index` argument.
circos.genomicPosTransformLines(bed,
                                posTransform = function(region, value)
                                posTransform.text(region, y = 1, labels = value[[1]],
                                                  cex = 0.5, padding = 1, track.index = tr_track),
                                direction = "inside", track.index = i_track
)


## link the trans pQTls to their target Somamer/gene 

circos.genomicLink(bed1,bed2, col = rand_color(nrow(bed1)), border= NA)


dev.off()

