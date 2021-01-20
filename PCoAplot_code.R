#Johnathan Shih
#7/12/2020


#Obtaining file
library(readxl)
distance_matrix <- read_excel("L:/GILAB/Johnathan/Plots/distance-matrix.xlsx")

#Performing PCoA
PC <- cmdscale(dist(distance_matrix))
x<- PC[,1]
y<- PC[,2]

#Defining factors
groups <- factor(distance_matrix$Dysbiosis)
cols = c('orange','red','green','darkgreen')[groups]
col1<-c("orange", "red", "green", "darkgreen")

groups
#Plotting PCoA
par(fig=c(0,0.8,0,0.8), new=FALSE)
plot(PC[,1],PC[,2],xlab= 'PC1', ylab= 'PC2',  pch=17, col = cols)


#Creating legend
legend("bottomleft",
       c("N Day 0","N Day 14","Y Day 0", "Y Day 14"),
       fill=c("orange","red","green","darkgreen"), cex = 0.8
)

#Plotting Boxplots
par(fig=c(0,0.8,0.55,1), new=TRUE)
boxplot(x ~ factor(distance_matrix$Dysbiosis), horizontal=TRUE, axes=FALSE, col = col1, outcolor = col1,outcex=0.5)
axis(side=1, labels=F, tcl=0, at=c(-19.8,31), tick=T, cex=0.7, lwd=1)
axis(side=2, labels=F, tcl=0, at=c(0.43,2.5), tick=T, cex=0.7, lwd=1)


par(fig=c(0.7,0.9,0,0.8), new=T)
boxplot(y ~ factor(distance_matrix$Dysbiosis), horizontal=FALSE, axes=FALSE, col = col1, outcolor = col1,outcex=0.5)
axis(side=1, labels=F, tcl=0, at=c(0.42,2.8), tick=T, cex=0.7, lwd=1)
axis(side=2, labels=F, tcl=0, at=c(-23.1,15), tick=T, cex=0.7, lwd=1)

mtext("PCoA on Distance matrix data (Dysbiosis plus time)", cex=1.5, side=3, outer=TRUE, line=-3)


