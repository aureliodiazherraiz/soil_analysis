## Standardised Major Axis Estimation 
## Stoichiometry Analysis
## Pablo Salazar Z
## University of Piura
## 17/09/2017

## Libraries required

library(plyr)
library(readxl)
library(readr)
library(xlsx)
library(ggthemes)
library(devtools)
#install_github("vqv/ggbiplot") ## Only for the first time
library(ggbiplot)
library(RColorBrewer)
library(Hmisc)
library(psych) ## to run corr.test
library(factoextra)
library(agricolae)
library(gridExtra)
library(ggplot2)
library(grid)
library(scales)
library(corrplot)
library(ggpubr)
library(lavaan) ## To make SEM 
library(AICcmodavg) ## To calculate AICC differences between models
library(semPlot)
library(tidyverse)


soil_data <- read_excel("E:/Dropbox/Fincyt/Estequiometria/Datos/Datos suelo.xlsx") ## Para laptop
PCAandLeaf_log <- read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/PCAandLeaf_log.csv")
PCAandSoil_log <- read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/PCAandSoil_log.csv")
soil_leaf_PCA12 <- read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/soil_leaf_PCA12.csv")

colnames(soil_data)[c(2)] <- c("plot")
colnames(soil_data)[c(7)] <- c("EC")
setwd("E:/Dropbox/Fincyt/Estequiometria/Analisis") ## Para laptop

names(soil_data)
## Open soil_and_traits.csv if you need functional traits too


## Variance explained ------------

VE <- read_excel("E:/Dropbox/Fincyt/Estequiometria/Datos/Variance explained.xlsx")
summary(VE)
names(VE)
l<-rep(1:10, times=3) ## order vector
l
VE<-cbind(VE,l)
## a is Tree, C is Population

a<- ggplot(data=VE[1:30,], aes(x=reorder(variable,l), y=variance, fill=scale)) + 
  geom_bar(stat="identity") + 
  xlab("\n Soil concentration variability") +
  ylab("Variance explained (%)")+
  labs(title="Soil concentration variability")+
  labs(fill=" ")

a<-a+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
           panel.background = element_blank(), axis.line = element_line(colour = "black"),
           axis.text=element_text(size=16),
           axis.title=element_text(size=20,face="bold"),
           plot.title=element_text(size=20,face="bold", hjust = 0.5),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           legend.position="top",
           axis.title.x = element_blank(),
           axis.text.x=element_text(angle = 45, hjust = 1)
)
#a <- a + scale_fill_brewer(palette = "Greens")
cols <- c("a" = "darkolivegreen1", "b" = "lightgreen", "c"="green4")
a<-a+scale_fill_manual(labels=c("Tree", "Plot", "Population"),values=cols)

a

pdf("Fig1aa.pdf", width=6, height=6)
par(mar=c(4.1,4.4,0.3,0.3))
a
dev.off()


b<- ggplot(data=VE[31:60,], aes(x=reorder(variable,l), y=variance, fill=scale)) + 
  geom_bar(stat="identity") + 
  xlab("\n leaf concentration variability") +
  ylab("Variance explained (%) \n") +
  labs(fill=" ")

b<-b+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
           panel.background = element_blank(), axis.line = element_line(colour = "black"),
           axis.text=element_text(size=16),
           axis.title=element_text(size=20,face="bold"),
           plot.title=element_text(size=20,face="bold", hjust = 0.5),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           legend.position="top",
           axis.title.x = element_blank(),
           axis.text.x=element_text(angle = 45, hjust = 1)
)

cols <- c("a" = "darkolivegreen1", "b" = "lightgreen", "c"="green4")
b<-b+scale_fill_manual(labels=c("Tree", "Plot", "Population"),values=cols)
b

pdf("Fig1bb.pdf", width=6, height=6)
par(mar=c(4.1,4.4,0.3,0.3))
b
dev.off()


### EC mean and statistics------------
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

CECSum <- summarySE(soil_data, measurevar="C.E.", groupvars=c("Zona"),na.rm=TRUE)
CECSum
pHSum <- summarySE(soil_data, measurevar="pH", groupvars=c("Zona"),na.rm=TRUE)
pHSum
write.csv(CECSum,file="cecsum.csv")
write.csv(pHSum,file="phsum.csv")
claySum <- summarySE(soil_data, measurevar="Clay", groupvars=c("Zona"),na.rm=TRUE)
sarSum <- summarySE(soil_data, measurevar="SAR", groupvars=c("Zona"),na.rm=TRUE)

CV<-sd(soil_data$C.E., na.rm=TRUE)/mean(soil_data$C.E., na.rm=TRUE)*100
cvpH<-sd(soil_data$pH, na.rm=TRUE)/mean(soil_data$pH, na.rm=TRUE)*100
cvN<-sd(soil_data$N, na.rm=TRUE)/mean(soil_data$N, na.rm=TRUE)*100
cvP<-sd(soil_data$P, na.rm=TRUE)/mean(soil_data$P, na.rm=TRUE)*100
cvK<-sd(soil_data$K, na.rm=TRUE)/mean(soil_data$K, na.rm=TRUE)*100
cvCa<-sd(soil_data$Ca, na.rm=TRUE)/mean(soil_data$Ca, na.rm=TRUE)*100
cvMg<-sd(soil_data$Mg, na.rm=TRUE)/mean(soil_data$Mg, na.rm=TRUE)*100
cvFe<-sd(soil_data$Fe, na.rm=TRUE)/mean(soil_data$Fe, na.rm=TRUE)*100
cvMn<-sd(soil_data$Mn, na.rm=TRUE)/mean(soil_data$Mn, na.rm=TRUE)*100
cvCu<-sd(soil_data$Cu, na.rm=TRUE)/mean(soil_data$Cu, na.rm=TRUE)*100
cvZn<-sd(soil_data$Zn, na.rm=TRUE)/mean(soil_data$Zn, na.rm=TRUE)*100

cvleaf<-c(cvN,cvP,cvK,cvCa,cvMg,cvFe,cvMn,cvCu,cvZn)
m<-mean(cvleaf)
cvleaf<-c(cvN,cvP,cvK,cvCa,cvMg,cvFe,cvMn,cvCu,cvZn,m)
write.csv(cvleaf,file="cv_leaf.csv")
names(soil_data)

cvN<-sd(soil_data$N_suelo, na.rm=TRUE)/mean(soil_data$N_suelo, na.rm=TRUE)*100
cvP<-sd(soil_data$P_suelo, na.rm=TRUE)/mean(soil_data$P_suelo, na.rm=TRUE)*100
cvK<-sd(soil_data$K_suelo, na.rm=TRUE)/mean(soil_data$K_suelo, na.rm=TRUE)*100
cvCa<-sd(soil_data$Ca_suelo, na.rm=TRUE)/mean(soil_data$Ca_suelo, na.rm=TRUE)*100
cvMg<-sd(soil_data$Mg_suelo, na.rm=TRUE)/mean(soil_data$Mg_suelo, na.rm=TRUE)*100
cvFe<-sd(soil_data$Fe_suelo, na.rm=TRUE)/mean(soil_data$Fe_suelo, na.rm=TRUE)*100
cvMn<-sd(soil_data$Mn_suelo, na.rm=TRUE)/mean(soil_data$Mn_suelo, na.rm=TRUE)*100
cvCu<-sd(soil_data$Cu_suelo, na.rm=TRUE)/mean(soil_data$Cu_suelo, na.rm=TRUE)*100
cvZn<-sd(soil_data$Zn_suelo, na.rm=TRUE)/mean(soil_data$Zn_suelo, na.rm=TRUE)*100

cvsoil<-c(cvN,cvP,cvK,cvCa,cvMg,cvFe,cvMn,cvCu,cvZn)
m<-mean(cvsoil)
cvsoil<-c(cvN,cvP,cvK,cvCa,cvMg,cvFe,cvMn,cvCu,cvZn,m)
write.csv(cvsoil,file="cv_soil.csv")
cv<-cbind(cvleaf,cvsoil)
Nut1<-c("N Soil", "P Soil","K Soil", "Ca Soil", "Mg Soil", "Fe Soil", "Mn Soil","Cu Soil", "Zn Soil", "Soil mean")
Nut2<-c("N Leaf", "P Leaf","K Leaf", "Ca Leaf", "Mg Leaf", "Fe Leaf", "Mn Leaf","Cu Leaf", "Zn Leaf", "Leaf mean")
cv<-cbind(cvleaf,cvsoil,Nut1,Nut2)
write.csv(cv,file="cv.csv")

cv <- read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/cv.csv")
summary(cv)

c<- ggplot(data=cv, aes(x=reorder(Nut1,X1), y=cvsoil)) + 
  geom_bar(stat="identity") + 
  xlab("\n Soil concentration variability") +
  ylab("Coefficient of variance (%)")
c<-c+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
           panel.background = element_blank(), axis.line = element_line(colour = "black"),
           axis.text=element_text(size=16),
           axis.title=element_text(size=20, face="bold"),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           ##legend.position="top",
           axis.title.x = element_blank(),
           axis.text.x=element_text(angle = 45, hjust = 1)
)
c<- c + scale_y_continuous(limits = c(0, 130))
c

d<- ggplot(data=cv, aes(x=reorder(Nut2,X1), y=cvleaf)) + 
  geom_bar(stat="identity") + 
  xlab("\n Soil concentration variability") +
  ylab("Coefficient of variance (%)")
d<-d+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
           panel.background = element_blank(), axis.line = element_line(colour = "black"),
           axis.text=element_text(size=16),
           axis.title=element_text(size=20),
           legend.text=element_text(size=16),
           legend.title=element_text(size=16),
           ##legend.position="top",
           axis.title.x = element_blank(),
           axis.text.x=element_text(angle = 45, hjust = 1)
)
d<- d + scale_y_continuous(limits = c(0, 130))
d

z<- ggarrange(b,d,labels = c("A)", "B)"),
              ncol=2,nrow=1,align="h",
              font.label = list(size = 22, face = "bold"),
              label.x=0.12,
              common.legend = TRUE,
              hjust=c(1,0.7)
                            )
z


tiff("Plot1_d.tiff", width = 12, height = 4.5, units = 'in', res = 300, compression = 'lzw')
z
dev.off()


################## Soil PCA Without C ------------------

log.ir <- log(soil_data[, 9:17]) ## Log Transformation of the variables included
colnames(log.ir)[1] <- "N"
colnames(log.ir)[2] <- "P"
colnames(log.ir)[3] <- "K"
colnames(log.ir)[4] <- "Ca"
colnames(log.ir)[5] <- "Mg"
colnames(log.ir)[6] <- "Fe"
colnames(log.ir)[7] <- "Mn"
colnames(log.ir)[8] <- "Cu"
colnames(log.ir)[9] <- "Zn"
is.na(log.ir) <- sapply(log.ir, is.infinite) ## This is to omit negative values after the log transformation
ir.zone <- soil_data[, 1] ## just renamed the categorical variable
ir.zone<- ir.zone[-c(58, 113),,drop=F] ## to omit the NAs in the categorical variable
ir.zone$Zona <- as.factor(ir.zone$Zona)
ir.pca <- prcomp(na.omit(log.ir),
                 center = TRUE,
                 scale = TRUE
) ## PCA operation with center ans scale active to reduce skewness effect
print(ir.pca)
plot(ir.pca,type="barplot")
summary(ir.pca)


ir.pca$x[,1] <- -ir.pca$x[,1] ## To flip the PC1
ir.pca$rotation[,1] <- -ir.pca$rotation[,1] ## Also hae to flip the eigenvectors
print(ir.pca)

# 
# ir.pca$x[,2] <- -ir.pca$x[,2] ## To flip the PC2
# ir.pca$rotation[,2] <- -ir.pca$rotation[,2] ## Also hae to flip the eigenvectors
# print(ir.pca)


ss<-get_eigenvalue(ir.pca)
test<-get_pca_var(ir.pca) ## variable contributioin to the PCA 
test$contrib
test$coord
test$cos2
test$cor

test2 <- get_pca_ind(ir.pca) ## Individual contribution to the PCA
test2$contrib
test2$coord
test2$cos2

Indi<-cbind(test2$coord,ir.zone)

log.ir2<-log.ir[-c(58, 113),,drop=F]
SP<-cbind(Indi,log.ir2)  
write.csv(SP,file="PCAandSoil_log_C.csv")


Anova<-aov(Dim.1~Zona,data=Indi)
Anova
summary(Anova)
TukeyHSD(Anova)
a<-HSD.test(Anova, "Zona", group=TRUE)
a
PCA1<-a$groups
IndiS <- Indi

Anova<-aov(Dim.2~Zona,data=Indi)
Anova
summary(Anova)
TukeyHSD(Anova)
b<-HSD.test(Anova, "Zona", group=TRUE)
b
PCA2<-b$groups

g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.zone$Zona, ellipse = TRUE, 
              circle = FALSE, repel =TRUE)
g<- g + labs(x=expression(PCA[1]~~group("(",63.3~"%",")")),
             y= expression(PCA[2]~~group("(",12.6~"%",")")))
g<-g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
           panel.background = element_blank(), axis.line = element_line(colour = "black"),
           axis.text=element_text(size=16),
           legend.text=element_text(size=16),
           legend.direction = 'horizontal',
           legend.position="top",
           axis.title.x=element_text(size=16),
           axis.title.y=element_text(size=16)
)

g <- g + scale_x_continuous(limits = c(-7,7))
g<- g + scale_color_brewer(palette="Paired")
g<-g + annotate("text", label="ND", y=-1.05, x= -6.2, colour= "red", size =2) 
g<-g+ annotate("text", label="PI", y=-0.87, x=-6.2, colour= "red", size =2)
g<-g + annotate("text", label="NW", y=-0.52, x=-6.2, colour= "red", size =2)
g<-g + annotate("text", label="QS", y=-0.02, x=-6.2, colour= "red", size =2)
g<-g + annotate("text", label="RS", y=0.19, x=-6.2, colour= "red", size =2)
g<-g + annotate("text", label="RI", y=0.34, x=-6.2, colour= "red", size =2)
g<-g+ annotate("text", label="LO", y=0.86, x=-6.2, colour= "red", size =2)
g<-g+ annotate("text", label="IT", y=1.12, x=-6.2, colour= "red", size =2)

g<-g + annotate("text", label="ND", y=-3, x=-2.2, colour= "red", size =2) 
g<-g+ annotate("text", label="PI", y=-3, x=0.74, colour= "red", size =2)
g<-g + annotate("text", label="NW", y=-3, x=-1.6, colour= "red", size =2)
g<-g + annotate("text", label="QS", y=-3, x=-2.6, colour= "red", size =2)
g<-g + annotate("text", label="RS", y=-3, x=3.1, colour= "red", size =2)
g<-g + annotate("text", label="RI", y=-3, x=1.7, colour= "red", size =2)
g<-g+ annotate("text", label="LO", y=-3, x=-0.25, colour= "red", size =2)
g<-g+ annotate("text", label="IT", y=-3, x=0.59, colour= "red", size =2)

g<- g + geom_segment(aes(x = -6.5, y = -1.2, xend = -6.5, yend = -0.3), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = -6.6, y = -0.9, xend = -6.6, yend = 0), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = -6.7, y = -0.6, xend = -6.7, yend = 0.45), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = -6.8, y = -0.1, xend = -6.8, yend = 0.98), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = -6.9, y = 0.23, xend = -6.9, yend = 1.2), colour= "black", size =0.8)

g<- g + geom_segment(aes(x = -2.7, y = -3.3, xend = -1.5, yend = -3.3), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = -2.3, y = -3.4, xend = -0.25, yend = -3.4), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = -0.4, y = -3.5, xend = 0.91, yend = -3.5), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = 0.58, y = -3.6, xend = 1.93, yend = -3.6), colour= "black", size =0.8)
g<- g + geom_segment(aes(x = 1.7, y = -3.7, xend = 3.32, yend = -3.7), colour= "black", size =0.8)

g


pdf("FigPCASoil.pdf", width=6, height=4.5)
par(mar=c(4.1,4.4,0.3,0.3))
g
dev.off()

colnames(test$cos2)[1] <- "PCA 1 \n (63 %)"
colnames(test$cos2)[2] <- "PCA 2 \n (13 %)"
colnames(test$cos2)[3] <- "PCA 3 \n (9 %)"

test_cos<-as.matrix(test$cos2, ncol=10)
test_cos<-test_cos[,(1:3)]


pdf("FigAnexSoilC.pdf", width=6, height=4.5)
par(mar=c(4.1,4.4,0.3,0.3))
corrplot(test_cos, is.corr=FALSE, win.asp =0.4, cl.cex = 1, cl.align.text = "l")
dev.off()

### Leaf PCA without C ----------

log.ir <- log(soil_data[, 29:37]) ## Log Transformation of the variables included
is.na(log.ir) <- sapply(log.ir, is.infinite) ## This is to omit negative values after the log transformation
ir.zone <- soil_data[, 1] ## just renamed the categorical variable
summary(log.ir)
which(is.na(log.ir))
log.ir<-log.ir[-c(36),,drop=F] ## to omit the crazy ND data
ir.zone<- ir.zone[-c(36),,drop=F] ## to omit the crazy ND data
ir.zone$Zona <- as.factor(ir.zone$Zona)
ir.pca <- prcomp(na.omit(log.ir),
                 center = TRUE,
                 scale = TRUE
) ## PCA operation with center ans scale active to reduce skewness effect
print(ir.pca)
plot(ir.pca,type="barplot")
summary(ir.pca)


ir.pca$x[,1] <- -ir.pca$x[,1] ## To flip the PC1
ir.pca$rotation[,1] <- -ir.pca$rotation[,1] ## Also hae to flip the eigenvectors
print(ir.pca)

ss<-get_eigenvalue(ir.pca)
test<-get_pca_var(ir.pca) ## variable contributioin to the PCA 
test$contrib
test$coord
test$cos2
test$cor
write.csv(test$coord,file= "PCAloadingLeaf.csv")

test2 <- get_pca_ind(ir.pca) ## Individual contribution to the PCA
test2$contrib
test2$coord
test2$cos2
summary(test2$coord)
##which(test2$coord<=-6.4000)
IndiLeaf<-cbind(test2$coord,ir.zone)
write.csv(IndiLeaf,file="PCAloadingIndividualsLeaf.csv")

SP<-cbind(IndiLeaf,log.ir)  

write.csv(SP,file="PCAandLeaf_log_C.csv")

Anova<-aov(Dim.1~Zona,data=IndiLeaf)
Anova
summary(Anova)
TukeyHSD(Anova)
a<-HSD.test(Anova, "Zona", group=TRUE)
a
PCA1<-a$groups


Anova<-aov(Dim.2~Zona,data=IndiLeaf)
Anova
summary(Anova)
TukeyHSD(Anova)
b<-HSD.test(Anova, "Zona", group=TRUE)
b
PCA2<-b$groups

h <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.zone$Zona, ellipse = TRUE, 
              circle = FALSE, repel =TRUE)
h<-h + labs(x=expression(PCA[1]~~group("(",46.7~"%",")")),
             y= expression(PCA[2]~~group("(",14.7~"%",")")))
h<-h+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
           panel.background = element_blank(), axis.line = element_line(colour = "black"),
           axis.text=element_text(size=16),
           legend.text=element_text(size=16),
           legend.direction = 'horizontal',
           legend.position="top",
           axis.title.x=element_text(size=16),
           axis.title.y=element_text(size=16)
)
h<-h + scale_color_brewer(palette="Paired")
h<-h + scale_x_continuous(limits = c(-7,7))
h<-h + annotate("text", label="ND", y=-0.1, x= -6.2, colour= "red", size =2) 
h<-h+ annotate("text", label="PI", y=-0.3, x=-6.2, colour= "red", size =2)
h<-h + annotate("text", label="NW", y=-1.3, x=-6.2, colour= "red", size =2)
h<-h + annotate("text", label="QS", y=0.07, x=-6.2, colour= "red", size =2)
h<-h + annotate("text", label="RS", y=0.2, x=-6.2, colour= "red", size =2)
h<-h + annotate("text", label="RI", y=0.6, x=-6.2, colour= "red", size =2)
h<-h+ annotate("text", label="LO", y=0.8, x=-6.2, colour= "red", size =2)
h<-h+ annotate("text", label="IT", y=0.4, x=-6.2, colour= "red", size =2)

##h<-h + annotate("text", label="ND", y=-2.9, x= 0.43, colour= "red", size =2) 
h<-h+ annotate("text", label="PI", y=-2.9, x=0.9, colour= "red", size =2)
h<-h + annotate("text", label="NW", y=-2.9, x=-0.03, colour= "red", size =2)
h<-h + annotate("text", label="QS", y=-2.9, x=-1.15, colour= "red", size =2)
h<-h + annotate("text", label="ND \n RS RI ", y=-2.9, x=0.3, colour= "red", size =2)
##h<-h + annotate("text", label="RI", y=-2.9, x=0.6, colour= "red", size =2)
h<-h+ annotate("text", label="LO", y=-2.9, x=1.9, colour= "red", size =2)
h<-h+ annotate("text", label="IT", y=-2.9, x=-2.8, colour= "red", size =2)

h<-h + geom_segment(aes(x = -6.5, y = 0.9, xend = -6.5, yend = -0.4), colour= "black", size =0.8)
h<-h + geom_segment(aes(x = -6.6, y = -0.05, xend = -6.6, yend = -1.4), colour= "black", size =0.8)

h<-h + geom_segment(aes(x = 2, y = -3.2, xend = 0.2, yend = -3.2), colour= "black", size =0.8)
h<-h + geom_segment(aes(x = 1.1, y = -3.3, xend = -0.3, yend = -3.3), colour= "black", size =0.8)
h<-h + geom_segment(aes(x = 0.8, y = -3.4, xend = -1.4, yend = -3.4), colour= "black", size =0.8)
h<-h + geom_segment(aes(x = -0.9, y = -3.5, xend = -2.9, yend = -3.5), colour= "black", size =0.8)

h


 
 pdf("FigPCAleaf_.pdf", width=6, height=4.5)
 par(mar=c(4.1,4.4,0.3,0.3))
 h
 dev.off()


colnames(test$cos2)[1] <- "PCA 1 \n (47 %)"
colnames(test$cos2)[2] <- "PCA 2 \n (15 %)"
colnames(test$cos2)[3] <- "PCA 3 \n (14 %)"

test_cos<-as.matrix(test$cos2, ncol=10)
test_cos<-test_cos[,(1:3)]


pdf("FigAnexLeafC.pdf", width=6, height=4.5)
par(mar=c(4.1,4.4,0.3,0.3))
corrplot(test_cos, is.corr=FALSE, win.asp =0.4, cl.cex = 1, cl.align.text = "l")
dev.off()

z<- ggarrange(g,h, nrow=2,align="h", labels = c("A) Soil", "B) Leaf"), vjust = c(20,20),
              font.label = list(size = 22, face = "bold"), hjust = -0.2,
              common.legend=TRUE, legend="top")
z

tiff("PlotPCA.tiff", width = 7.5, height = 10, units = 'in', res = 300, compression = 'lzw')
z
dev.off()

#### Restarting the soil fertility data to exclude C --------------

PCASoil<-read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/PCAandSoil_log_C.csv") ## Laptop
PCALeaf<-read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/PCAandLeaf_log_C.csv") ## laptop

index<-(1:120)
index<-matrix(index, ncol=1)
index_S<-index[-c(58, 113),,drop=F]
index_L<-index[-c(36),,drop=F]

colnames(PCASoil) <- paste("Soil", colnames(PCASoil), sep = " ")
colnames(PCASoil)[2]<-"S1"
colnames(PCASoil)[3]<-"S2"

colnames(PCALeaf) <- paste("Leaf", colnames(PCALeaf), sep = " ")
colnames(PCALeaf)[2]<-"L1"
colnames(PCALeaf)[3]<-"L2"

names(PCASoil)

soil_data2<-cbind(PCASoil[11:20],PCASoil$S1,PCASoil$S2,index_S)
colnames(soil_data2)[11]<-"S1 (63%)"
colnames(soil_data2)[12]<-"S2 (13%)"
colnames(soil_data2)[13]<-"index"

leaf_data<-cbind(PCALeaf[11:20],PCALeaf$L1,PCALeaf$L2,index_L)

colnames(leaf_data)[11]<-"L1 (47%)"
colnames(leaf_data)[12]<-"L2 (15%)"
colnames(leaf_data)[13]<-"index"

soil_leaf<-merge(soil_data2,leaf_data,by="index", all=T)
soil_leaf<-cbind(soil_leaf[2:13],soil_data$pH,soil_data$EC_P,soil_data$SAR,soil_leaf[24:25],soil_leaf[15:23])
soil_leaf<-cbind(soil_leaf[2:13],soil_data$pH,soil_data$EC_P,soil_data$SAR,soil_leaf[24:25],soil_leaf[15:23])
colnames(soil_leaf)[13]<-"pH"
colnames(soil_leaf)[14]<-"EC"
colnames(soil_leaf)[15]<-"SAR"
write.csv(soil_leaf,"soil_leaf_PCA_C.csv") ## This has log transformation

## correlation between soil, leaf and both PCA --------------

SL <- read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/soil_leaf_PCA_C.csv")
SL1<-SL[,-1,drop=F]

res <- rcorr(as.matrix(SL1))

res

r<-res$r
p<-res$P

write.csv(res$r,"res_C.csv")
write.csv(res$P,"res2_C.csv")

## define notions for significance levels; spacing is important.
mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ",ifelse(p<.1,"a ", " "))))

## trunctuate the matrix that holds the correlations to two decimal
R <- format(round(cbind(rep(-1.11, ncol(SL1)), r), 2))[,-1] 

## build a new matrix that includes the correlations with their apropriate stars 
Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(SL1)) 
diag(Rnew) <- paste(diag(R), " ", sep="") 
rownames(Rnew) <- colnames(1) 
colnames(Rnew) <- paste(colnames(SL1), "", sep="")

Rnew

write.csv(Rnew,"Cor_C.csv")

## Climatic correlations -------

SL2<-SL[,c(11,12,20,21)]
SL2<- cbind(soil_data$Zona,SL2)
colnames(SL2) <- c("Zone", "S1", "S2", "L1", "L2")
SL3 <- SL2 %>%  group_by(Zone) %>%  summarise_all(funs(mean))
SL3 <- SL2 %>%  group_by(Zone) %>%  summarise_all(list(~mean(., na.rm=TRUE)))

SL <- read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/soil_leaf_PCA_C.csv")

## Corr plot about the correlations --------------
summary(res)

tiff("Figcoors_C.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
corrplot(res$r, method = "ellipse",p.mat=res$P, insig ="blank", type = "upper"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)

dev.off()



## PCA vs plant nutrition correlations -----------------
fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}


## Leaf table of mean values ----------
names(soil_data)
colnames(soil_data)[29] <- "Ni"
MeanC <- summarySE(soil_data, measurevar="C", groupvars=c("Zona"),na.rm=TRUE)
MeanN <- summarySE(soil_data, measurevar="Ni", groupvars=c("Zona"),na.rm=TRUE)
MeanP <- summarySE(soil_data, measurevar="P", groupvars=c("Zona"),na.rm=TRUE)
MeanK <- summarySE(soil_data, measurevar="K", groupvars=c("Zona"),na.rm=TRUE)
MeanCa <- summarySE(soil_data, measurevar="Ca", groupvars=c("Zona"),na.rm=TRUE)
MeanMg <- summarySE(soil_data, measurevar="Mg", groupvars=c("Zona"),na.rm=TRUE)
MeanFe <- summarySE(soil_data, measurevar="Fe", groupvars=c("Zona"),na.rm=TRUE)
MeanMn <- summarySE(soil_data, measurevar="Mn", groupvars=c("Zona"),na.rm=TRUE)
MeanCu <- summarySE(soil_data, measurevar="Cu", groupvars=c("Zona"),na.rm=TRUE)
MeanZn <- summarySE(soil_data, measurevar="Zn", groupvars=c("Zona"),na.rm=TRUE)

MeanC$X<-paste(round(MeanC$C,digits=2),"\u00B1",round(MeanC$se,digits=2))
MeanCu$X<-paste(round(MeanCu$Cu,digits=2),"\u00B1",round(MeanCu$se,digits=2))
MeanFe$X<-paste(round(MeanFe$Fe,digits=2),"\u00B1",round(MeanFe$se,digits=2))
MeanMn$X<-paste(round(MeanMn$Mn,digits=2),"\u00B1",round(MeanMn$se,digits=2))
MeanN$X<-paste(round(MeanN$Ni,digits=2),"\u00B1",round(MeanN$se,digits=2))
MeanK$X<-paste(round(MeanK$K,digits=2),"\u00B1",round(MeanN$se,digits=2))
MeanCa$X<-paste(round(MeanCa$Ca,digits=2),"\u00B1",round(MeanN$se,digits=2))
MeanMg$X<-paste(round(MeanMg$Mg,digits=2),"\u00B1",round(MeanN$se,digits=2))
MeanP$X<-paste(round(MeanP$P,digits=2),"\u00B1",round(MeanP$se,digits=2))
MeanZn$X<-paste(round(MeanZn$Zn,digits=2),"\u00B1",round(MeanZn$se,digits=2))

Means<-data.frame(MeanC$Zona,MeanC$X,MeanN$X,MeanP$X, MeanK$X,MeanCa$X,
                  MeanMg$X,MeanFe$X,MeanMn$X,MeanCu$X,MeanZn$X)

Means<-rbind(Means,cv)

write.csv(Means,"Means_leaf.csv")


## CE and pH graphs ------

names(soil_data)
# soil_leaf_PCA12 <- cbind(soil_data,soil_leaf[11:12],soil_leaf[16:17])
# write.csv(soil_leaf_PCA12,"soil_leaf_PCA12.csv")
soil_leaf_PCA12 <- read_csv("E:/Dropbox/Fincyt/Estequiometria/Analisis/soil_leaf_PCA12.csv")
names(soil_leaf_PCA12)
colnames(soil_leaf_PCA12)[53] <- "S1"
colnames(soil_leaf_PCA12)[54] <- "S2"
colnames(soil_leaf_PCA12)[55] <- "L1"
colnames(soil_leaf_PCA12)[56] <- "L2"
colnames(soil_leaf_PCA12)[52] <- "EC"

temas <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               panel.border =element_rect(size = 0.5, linetype="solid",fill=NA, colour ="black"),
               axis.text=element_text(size=16),
               axis.title=element_text(size=20),
               legend.text=element_text(size=16),
               legend.title=element_text(size=16),
               legend.position="none"
)

p <- ggplot(soil_leaf_PCA12, aes(x = EC))
p <- p + geom_point(aes(y = S1, colour = "PC1 S"), size = 2)
p <- p + geom_point(aes(y = L1, colour = "PC1 L"), size = 2)
p<- p + geom_smooth(aes(EC,S1), color= "brown", method=lm, se=FALSE)
p<- p + geom_smooth(aes(EC,L1), color= "darkgreen", method=lm, se=FALSE)
cols <- c("PC1 S" = "brown", "PC1 L" = "darkgreen")
p <- p + scale_colour_manual(values = cols, labels=c("PCA1 L", "PCA1 S"))
p <- p + labs(y = expression(bold(PCA~loading)),
              x = expression(bold(Electrical~Conductivity~group("(",mu*S~cm^{-1},")"))),
              colour = "")
p
p<-p+ temas
p<- p + theme(axis.text=element_blank(),axis.title = element_blank())
p<-p+geom_text(aes(x=1300, y=-3), size=8, 
               label="S~r==0.34*'***'", color="brown", parse=TRUE)
p<-p+geom_text(aes(x=1300, y=-4), size=8, 
               label="L~r==0.19*'***'", color="darkgreen", parse=TRUE)
p


q <- ggplot(soil_leaf_PCA12, aes(x = pH))
q <- q + geom_point(aes(y = S2, colour = "PC2 S"), size = 2)
q <- q + geom_point(aes(y = L2, colour = "PC2 L"), size = 2)
q <- q + scale_colour_manual(values = c("darkgreen", "brown"))
q <- q + labs(y = expression(bold(PCA[2]~loading)),
              x = expression(bold(Soil~pH)),
              colour = " ")
q<- q + geom_smooth(aes(pH,S2), color= "brown", method=lm, se=FALSE)
q<- q + geom_smooth(aes(pH,L2), color= "darkgreen", method=lm, se=FALSE)
q<-q+temas + theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))
q<-q+geom_text(aes(x=8.4, y=2.7), size=8, 
               label="S~r==-0.56*'***'", color="brown",parse=TRUE)
q<-q+geom_text(aes(x=8.2, y=3.3), size=8, 
               label="L~r==-0.19*'*'", color="darkgreen",parse=TRUE)
q<-q+scale_x_continuous(limits=c(6,9))
q

r <- ggplot(soil_leaf_PCA12, aes(x = EC))
r <- r + geom_point(aes(y = S2, colour = "PC2 S"), size = 2)
r <- r + geom_point(aes(y = L2, colour = "PC2 L"), size = 2)
r<- r + geom_smooth(aes(EC,S2), color= "brown", method=lm, se=FALSE)
r<- r + geom_smooth(aes(EC,L2), color= "darkgreen", method=lm, se=FALSE)
r <- r + scale_colour_manual(values = c("darkgreen", "brown"))
r <- r + labs(y = expression(bold(PCA~loading)),
              x = expression(bold(Electrical~Conductivity~group("(",mu*S~cm^{-1},")"))),
              colour = " ")
r<-r+temas
r<- r + theme(axis.text.y=element_blank(),axis.title.y = element_blank(),
              axis.title.x = element_text(margin = margin(t = 5, r = 0, b =0, l = 0)))
r<-r+geom_text(aes(x=1300, y=2), size=8, 
               label="S~r==-0.50*'***'", color="brown",parse=TRUE)
r<-r+geom_text(aes(x=1300, y=3), size=8, 
               label="L~r==-0.22*'*'", color="darkgreen",parse=TRUE)
r

w <- ggplot(soil_leaf_PCA12, aes(x = pH))
w <- w + geom_point(aes(y = S1, colour = "PC1 S"), size = 2)
w <- w + geom_point(aes(y = L1, colour = "PC1 L"), size = 2)
w <- w + scale_colour_manual(values = c("darkgreen", "brown"))
w <- w + labs(y = expression(bold(PCA[1]~loading)),
              x = expression(bold(Soil~pH)),
              colour = " ")
w<- w + geom_smooth(aes(pH,S1), color= "brown", method=lm, se=FALSE)
w<- w + geom_smooth(aes(pH,L1), color= "darkgreen", method=lm, se=FALSE)
w<-w+temas
w<- w + theme(axis.text.x=element_blank(),axis.title.x = element_blank())
w<-w+geom_text(aes(x=8.3, y=4), size=8, 
               label="S~r==-0.49*'***'", color="brown",parse=TRUE)
w<-w+geom_text(aes(x=8.2, y=5), size=8, 
               label="L~r==-0.03", color="darkgreen",parse=TRUE)
w<-w+scale_x_continuous(limits=c(6,9))

w


s <- ggplot(soil_leaf_PCA12, aes(x = SAR)) + 
  geom_point(aes(y = S1, colour = "PC1 S"), size = 2) + 
  geom_point(aes(y = L1, colour = "PC1 L"), size = 2) + 
  scale_colour_manual(values = c("darkgreen", "brown")) + 
  labs(y = expression(bold(PCA~loading)),
              x = expression(bold(Soil~SAR)),
              colour = " ") + 
  geom_smooth(aes(SAR,S1), color= "brown", method=lm, se=FALSE) + 
  geom_smooth(aes(SAR,L1), color= "darkgreen", method=lm, se=FALSE) + 
  temas
s<- s + theme(axis.text=element_blank(),axis.title = element_blank())
s<-s+  geom_text(aes(x=3, y=4), size=8, 
               label="S~r==-0.07", color="brown",parse=TRUE) + 
  geom_text(aes(x=3, y=5), size=8, 
               label="L~r==-0.16", color="darkgreen",parse=TRUE)
s

v <- ggplot(soil_leaf_PCA12, aes(x = SAR)) + 
  geom_point(aes(y = S2, colour = "PC2 S"), size = 2) + 
  geom_point(aes(y = L2, colour = "PC2 L"), size = 2) + 
  scale_colour_manual(values = c("darkgreen", "brown")) + 
  labs(y = expression(bold(PCA~loading)),
       x = expression(bold(Soil~SAR)),
       colour = " ") + 
  geom_smooth(aes(SAR,S2), color= "brown", method=lm, se=FALSE) + 
  geom_smooth(aes(SAR,L2), color= "darkgreen", method=lm, se=FALSE) + 
  temas
v<- v + theme(axis.text.y=element_blank(),axis.title.y = element_blank(),
              axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))
v<-v+geom_text(aes(x=3, y=2.5), size=8, 
               label="S~r==-0.54*'***'", color="brown",parse=TRUE)
v<- v+geom_text(aes(x=3, y=3.5), size=8, 
            label="L~r==-0.42*'***'", color="darkgreen",parse=TRUE)
#v<-v+scale_x_continuous(limits=c(6,9))
v

z<- ggarrange(w,p,s,q,r,v, ncol=3, nrow=2, 
              labels = c("A)", "B)", "C)","D)","E)","F)"), hjust=c(-2.5,-0.5,-0.5,-2.5,-0.5,-0.5),
              font.label = list(size = 22, face = "bold"),
              common.legend = FALSE)

tiff("Plot5c.tiff", width = 18, height = 9, units = 'in', res = 300, compression = 'lzw')
z
dev.off()

## SEM analysis --------

## USE soil_leaf_PCA12

colnames(soil_leaf_PCA12)[53] <- "S1"
colnames(soil_leaf_PCA12)[54] <- "S2"
colnames(soil_leaf_PCA12)[55] <- "L1"
colnames(soil_leaf_PCA12)[56] <- "L2"
colnames(soil_leaf_PCA12)[52] <- "EC"

SL1<-soil_leaf_PCA12

names(SL1)

lvmod.1<-"
# Regressions
S1 ~ pH
S2 ~ pH + EC
SAR ~ EC
L1 ~ EC
L2 ~ pH + SAR
#Error covariance
"

lvmod.1.fit<-sem(lvmod.1,data=SL1, estimator ="ML", std.lv = TRUE)
varTable(lvmod.1.fit) 
EC_T<-SL1$EC/100
SL3 <- cbind(SL1,EC_T)

lvmod.2<-"
# Regressions
S1 ~ pH + EC_T
S2 ~ pH + EC_T
SAR ~ EC_T + pH
L1 ~ EC_T
L2 ~ SAR
#Error covariance
S1 ~~ SAR
"

lvmod.2.fit<-sem(lvmod.2,data=SL3, estimator ="ML", std.lv = TRUE)
varTable(lvmod.2.fit) 

summary(lvmod.2.fit, rsq=T,standardized=T) 
mi<-modindices(lvmod.2.fit) 
print(mi[mi$mi>3.0,]) 
mi

lbls<-c("Soil[N]","Soil[Ca]","Soil[P]","Soil[Mg]","pH","Leaf[N]","Leaf[P]",
        "Leaf[K]","Leaf[Ca]", "Leaf[Mg]", "Soil[Ca]", "EC", "Leaf[Mn]", "Soil\n Fertility", "Plant\n Nutrition")
grps<-list(Plant=c("N","P", "K", "Ca", "Mg","PlantNutrition"),
           Soil=c("N_suelo","P_suelo_Trans","K_suelo_Trans","Mg_suelo_Trans","Ca_suelo_Trans", "SoilFertility"),
           Prop=c("EC_Trans","pH"),
           Mn="Mn_Trans")
tiff("FigSem2_.tiff", width=6, height=4.5, res= 300, units="in", compression ="lzw")
par(mar=c(4.1,4.4,0.3,0.3))
semPaths(lvmod.2.fit,layout="circle3",what="std",residuals=FALSE, exoCov=FALSE, fade=FALSE,
         edge.label.position=0.20, nCharNodes=0, bifactor = "EC", intercepts = FALSE,
  #       nodeLabels=lbls, groups=grps,color=c("green4", "brown4","yellow","gray"), 
         sizeMan=10,sizeMan2=6,edge.label.cex=0.8, legend=FALSE)
text(0.9,0.9,labels="Chi-square = 2 \n df = 6 \n p = 0.9")
dev.off()



lvmod.3<-"
# Regressions

S1 ~ pH + EC_T + Clay
S2 ~ pH + EC_T
SAR ~ EC_T + pH
L1 ~ EC_T
L2 ~ SAR
#Error covariance
S1 ~~ SAR
SAR ~~ Clay
S2 ~~ Clay
"

lvmod.3.fit<-sem(lvmod.3,data=SL3, estimator ="ML", std.lv = TRUE, test = "bootstrap")
varTable(lvmod.3.fit) 

summary(lvmod.3.fit, rsq=T,standardized=T) 
mi<-modindices(lvmod.3.fit) 
print(mi[mi$mi>3.0,]) 
mi
fitmeasures(lvmod.3.fit)

lbls<-c("Soil[N]","Soil[Ca]","Soil[P]","Soil[Mg]","pH","Leaf[N]","Leaf[P]",
        "Leaf[K]","Leaf[Ca]", "Leaf[Mg]", "Soil[Ca]", "EC", "Leaf[Mn]", "Soil\n Fertility", "Plant\n Nutrition")
grps<-list(Plant=c("N","P", "K", "Ca", "Mg","PlantNutrition"),
           Soil=c("N_suelo","P_suelo_Trans","K_suelo_Trans","Mg_suelo_Trans","Ca_suelo_Trans", "SoilFertility"),
           Prop=c("EC_Trans","pH"),
           Mn="Mn_Trans")
tiff("FigSem3_.tiff", width=6, height=4.5, res= 300, units="in", compression ="lzw")
par(mar=c(4.1,4.4,0.3,0.3))
semPaths(lvmod.3.fit,layout="circle",what="std",residuals=FALSE, exoCov=FALSE, fade=FALSE,
         edge.label.position=0.20, nCharNodes=0, #bifactor = "S1",
         #       nodeLabels=lbls, groups=grps,color=c("green4", "brown4","yellow","gray"), 
         sizeMan=10,sizeMan2=6,edge.label.cex=0.8, legend=FALSE)
text(0.9,-0.9,labels="Chi-square = 7.2 \n df = 9 \n p = 0.6")
dev.off()


## New SEM model with latent variables -----------


colnames(soil_leaf_PCA12)[53] <- "S1"
colnames(soil_leaf_PCA12)[54] <- "S2"
colnames(soil_leaf_PCA12)[55] <- "L1"
colnames(soil_leaf_PCA12)[56] <- "L2"
colnames(soil_leaf_PCA12)[52] <- "EC"

SL1<-soil_leaf_PCA12

names(SL1)


lvmod.1<-"
#latent variable

Soil =~ N_suelo_T + P_suelo + K_suelo_T + Ca_suelo_T + Mg_suelo_T + Fe_suelo + Mn_suelo + Cu_suelo_T + Zn_suelo
leaf =~ N + P + K + Mg + Fe_T + Mn_T + Cu + Zn

# Regressions

Soil ~ pH + EC_T + Clay + SAR
SAR ~ EC_T + pH
leaf ~ EC_T



# residual covariance

Mn_suelo ~~ Cu_suelo_T
N ~~  P
leaf ~~ SAR + Soil
Cu ~~ Zn
P_suelo ~~ K_suelo_T
Ca_suelo_T ~~   Mn_suelo
N_suelo_T ~~    P_suelo
"

varTable(fit) 
N_suelo_T <- SL1$N_suelo*10
Ca_suelo_T <- SL1$Ca_suelo/1000
K_suelo_T <- SL1$K_suelo/100
Mg_suelo_T <- SL1$Mg_suelo/100
Cu_suelo_T <- SL1$Cu_suelo*10
Fe_T <- SL1$Fe/100
Mn_T <- SL1$Mn/10
EC_T<-SL1$EC/100

SL2 <- cbind(SL1, N_suelo_T, Ca_suelo_T, K_suelo_T, Mg_suelo_T, Cu_suelo_T, Mn_T, Fe_T,EC_T)

fit1<- sem(lvmod.1, data=SL2)

summary(fit1, rsq=T,standardized=T) 
mi<-modindices(fit1) 
print(mi[mi$mi>10.0,]) 
mi



lvmod.2<-"
#latent variable

SoilO =~ N_suelo_T + P_suelo_T + K_suelo_T
leafO =~ N + P + K 

SoilI =~ Ca_suelo_T + Mg_suelo_T + Fe_suelo_T + Mn_suelo + Cu_suelo_T + Zn_suelo
leafI =~ Mg + Fe_T + Mn_T + Cu + Zn

# Regressions

EC_T ~ SoilO
SoilI ~ Clay + pH
SAR ~ EC_T + pH +SoilO
pH ~ SoilO
leafO ~ EC_T + leafI
leafI ~ SAR 

# residual covariance

Mn_suelo ~~ Cu_suelo_T
N ~~  P
Cu ~~ Zn
P_suelo_T ~~ K_suelo_T
Ca_suelo_T ~~ Mn_suelo
N_suelo_T ~~  P_suelo_T
SoilI ~~ leafI 
SoilI ~~ SoilO
Mg_suelo_T ~~ EC_T
Fe_suelo_T ~~ Cu_suelo_T
SoilO ~~ SAR
leafI ~~ EC_T
"

varTable(fit) 
N_suelo_T <- SL1$N_suelo*10
Ca_suelo_T <- SL1$Ca_suelo/1000
K_suelo_T <- SL1$K_suelo/100
Mg_suelo_T <- SL1$Mg_suelo/100
Cu_suelo_T <- SL1$Cu_suelo*10
P_suelo_T <- SL1$P_suelo/10
Fe_suelo_T <- SL1$Fe_suelo/10
Fe_T <- SL1$Fe/100
Mn_T <- SL1$Mn/10
EC_T<-SL1$EC/100

SL2 <- cbind(SL1, N_suelo_T, Ca_suelo_T, K_suelo_T, Mg_suelo_T, P_suelo_T, Fe_suelo_T, Cu_suelo_T, Mn_T, Fe_T,EC_T)

fit2<- sem(lvmod.2, data=SL2)

summary(fit2, rsq=T,standardized=T) 
mi<-modindices(fit2) 
print(mi[mi$mi>10.0,]) 
mi


lvmod.3<-"
#Composite variables

Soil <~ N_suelo_T + P_suelo_T + K_suelo_T + Ca_suelo_T + Mg_suelo_T + Fe_suelo_T + Mn_suelo + Cu_suelo_T + Zn_suelo
leaf =~ N + P + K + Mg + Fe_T + Mn_T + Cu + Zn

# Regressions

leaf ~ EC_T + pH + Clay + SAR + Soil

"


varTable(fit) 
N_suelo_T <- SL1$N_suelo*10
Ca_suelo_T <- SL1$Ca_suelo/1000
K_suelo_T <- SL1$K_suelo/100
Mg_suelo_T <- SL1$Mg_suelo/100
Cu_suelo_T <- SL1$Cu_suelo*10
P_suelo_T <- SL1$P_suelo/10
Fe_suelo_T <- SL1$Fe_suelo/10
Fe_T <- SL1$Fe/100
Mn_T <- SL1$Mn/10
EC_T<-SL1$EC/100

SL2 <- cbind(SL1, N_suelo_T, Ca_suelo_T, K_suelo_T, Mg_suelo_T, P_suelo_T, Fe_suelo_T, Cu_suelo_T, Mn_T, Fe_T,EC_T)

fit3<- sem(lvmod.3, data=SL2)

summary(fit3, rsq=T,standardized=T) 
mi<-modindices(fit2) 
print(mi[mi$mi>10.0,]) 
mi


lvmod.4<-"
#Composite variables

Soil <~ N_suelo_T + P_suelo_T + K_suelo_T + Ca_suelo_T + Mg_suelo_T + Fe_suelo_T + Mn_suelo + Cu_suelo_T + Zn_suelo
leaf =~ N + P + K + Mg + Fe_T + Mn_T + Cu + Zn

# Regressions

leaf ~ Soil

"

fit4<- sem(lvmod.4, data=SL2)

summary(fit4, rsq=T,standardized=T) 
