#Johnathan Shih
# 7/29/2020

#Based off Decontam tutorial from https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

#Installing packages
# install.packages("jsonlite", repos="http://cran.r-project.org")
# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("rhdf5")
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq); packageVersion("phyloseq")

#install.packages("ggplot2")
library(ggplot2); packageVersion("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("decontam")
library(decontam); packageVersion("decontam")






#Creating phyloseq object
library(readxl)

ta <- read.csv("L:/GILAB/Johnathan/Decontam_script/table.from_biom.csv", header =T)
tb <- read.csv("L:/GILAB/Johnathan/Decontam_script/032420JSillcus515F-mapping.csv", header = T)


row.names(ta) <- ta$sample
pa <- as.matrix(ta)
class(pa) <- "numeric"

pt <- subset(pa, select = -c(sample))



row.names(tb) <- tb$sample
tb <- subset(tb, select = -c(sample))

OTU = otu_table(pt, taxa_are_rows = FALSE)
sampledata = sample_data(tb)
pss = phyloseq(OTU, sampledata)


# sample_names(sampledata)
# sample_names(OTU)
# sample_data(pss)


df <- as.matrix(OTU) # Put sample_data into a ggplot-friendly data.frame


# Testing for contaminents

#Frequency method
conc <- tb$Concentration
conc

contamdf.freq <- isContaminant(pss, method="frequency", conc=conc, threshold = 0.5)

head(contamdf.freq)

table(contamdf.freq$contaminant)

table(which(contamdf.freq$contaminant))


#Prevalence method

sample_data(pss)$is.neg <- sample_data(pss)$Sample_or_Control == "CONTROL SAMPLE"
contamdf.prev <- isContaminant(pss, method="prevalence", neg="is.neg", threshold = 0.5)
table(contamdf.prev$contaminant)
table(which(contamdf.prev$contaminant))
head(contamdf.prev)



#Exporting data
contam <- as.data.frame(contamdf.prev)
library(xlsx)
write.xlsx2(contam, "L:/GILAB/Johnathan/Decontam_script/contam_prevalence.xlsx")



#Plotting Frequency vs DNA concentration

plot_frequency(pss, taxa_names(pss)[c(21,29,31,54,63,73,87,127,141,186,259,301,328,343,385,397,412,421,428,450,552,558,560,608,615,623,671,726,744,760,807,810,824,827,830,849,869,885,890,946,980,1016,1028,1107,1127,1143)], conc=conc) + 
  xlab("DNA Concentration ")
# plot_frequency(pss, taxa_names(pss)[c(1,2,3,4,5,6,25,80)], conc=conc) + 
#   xlab("DNA Concentration ")



#Plotting Prevalence

pss.pa <- transform_sample_counts(pss, function(abund) 1*(abund>0))
pss.pa.neg <- prune_samples(sample_data(pss.pa)$Sample_or_Control == "CONTROL SAMPLE", pss.pa)

pss.pa.pos <- prune_samples(sample_data(pss.pa)$Sample_or_Control == "TRUE SAMPLE", pss.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(pss.pa.pos), pa.neg=taxa_sums(pss.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
