
setwd(dir="C:/Users/Aurelio Diaz/OneDrive - Universidad de Córdoba/Doctorado/Iberian_forest/Soil_analisys/inputs_database")

datsoil <- read.table("290520_Resumen_suelos forestales_HR.csv", stringsAsFactors = F, sep = ";", header = T, dec = ",")
str(datsoil)

rownames(datsoil) = datsoil$ï..Local
head(datsoil)

#remove the first column because now have the local names in the index
datsoil$ï..Local<-NULL
datsoil<-datsoil[, -c(1,2)]
head(datsoil)
chart.Correlation(datsoil)

#we can study the possible correlations between variables
cor(datsoil)

#applied in the columns(2) of the dataset the statistical variance
apply(datsoil, 2, var)

#like have a hight variance in some variables, apply log transformation anda reapply variance
datsoil1<-log(datsoil)
apply(datsoil1,2, var)
head(datsoil1)

#another option maybe centrering anda scaling the variables with prcomp
pca.soil<-prcomp(datsoil, center = T, scale. = T)
print(pca.soil)

fviz_eig(pca.soil)
summary(pca.soil)#80% the variance proportion is explained by the principles four components 

#we try apply prcomp anda comparise with the logistic transformation (in datsoil1)
# but we dont got "Error in svd(x, nu = 0, nv = k) : infinite or missing values in 'x'"
#pca.soil1<-prcomp(datsoil1, center = T, scale. = T)
#print(pca.soil1)
cor(datsoil1)
chart.Correlation(datsoil1)

#for graphics
biplot(pca.soil, scale = 0)

fviz_pca_var(pca.soil,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


eig.val <- get_eigenvalue(pca.soil)
eig.val





''' 
#cambiamos el nombre de las filas por la columa local previamente montada, 
#no puede haber valores repetidos lo cual me ha complicado la vida.

datdig1<-datdig1[,-c(1,5)]
head(datdig1)

#realizamos una transformacion en los par?metros de mayores valores absolutos
# asi al m?ximo del K, Ca y Mn se le retira el valor de cada muestra
datdig1$k<-max(datdig1$k) -datdig1$k
datdig1$Na<-max(datdig1$Na) -datdig1$Na
datdig1$Ca<-max(datdig1$Ca) -datdig1$Ca

head(datdig1)


#para ver las diferentes variables
library(TeachingDemos) 

faces2(datdig1, nrows=7)


panel.hist <- function(x, ...) 
{ 
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(usr[1:2], 0, 1.5) ) 
  h <- hist(x, plot = FALSE) 
  breaks <- h$breaks; nB <- length(breaks) 
  y <- h$counts; y <- y/max(y) 
  rect(breaks[-nB], 0, breaks[-1], y, col="blue", ...) 
} 
pairs(datdig1,diag.panel=panel.hist)

library(PerformanceAnalytics)
chart.Correlation(datdig1)
rs.cor <- (cor(x=datdig1, method="pearson",use="complete.obs"))
rs.cor

cor(datdig1)
dig.pc<-princomp(datdig1,cor=TRUE) 
summary(dig.pc,loadings=TRUE)

# Es lo mismo que calcular los autovalores y autovectores de S
S = cor(datdig1) 
eigen(S)

dig.pc$scores[,1:3]

par(pty="s") 
plot(dig.pc$scores[,1],dig.pc$scores[,2], 
     ylim=range(dig.pc$scores[,1]), 
     xlab="PC1",ylab="PC2",type="n",lwd=2) 
text(dig.pc$scores[,1],dig.pc$scores[,2], 
     labels=abbreviate(row.names(datdig1)),cex=0.5,lwd=2)


plot(dig.pc$scores[,1],dig.pc$scores[,3], 
     ylim=range(dig.pc$scores[,1]), 
     xlab="PC1",ylab="PC3",type="n",lwd=4) 
text(dig.pc$scores[,1],dig.pc$scores[,3], 
     labels=abbreviate(row.names(datdig1)),cex=0.5,lwd=4)


plot(dig.pc$scores[,2],dig.pc$scores[,3], 
     ylim=range(dig.pc$scores[,2]), 
     xlab="PC2",ylab="PC3",type="n",lwd=2) 
text(dig.pc$scores[,2],dig.pc$scores[,3], 
     labels=abbreviate(row.names(datdig1)),cex=0.5,lwd=2)


#para relacionarlo con una de las variables
par(mfrow=c(2,3)) 
plot(dig.pc$scores[,1], datdig$P, xlab="PC1") 
plot(dig.pc$scores[,2], datdig$P, xlab="PC2") 
plot(dig.pc$scores[,3], datdig$P, xlab="PC3")

plot(dig.pc$scores[,1], datdig$Ca, xlab="PC1") 
plot(dig.pc$scores[,2], datdig$Ca, xlab="PC2") 
plot(dig.pc$scores[,3], datdig$Ca, xlab="PC3")

summary(lm(datdig$P~dig.pc$scores[,1]+dig.pc$scores[,2]+dig.pc$scores[,3]))
summary(lm(datdig$Ca~dig.pc$scores[,1]+dig.pc$scores[,2]+dig.pc$scores[,3]))

plot(dig.pc$scores[,2],datdig$P,xlab="PC1",ylab="P") 


pca <- prcomp(datdig1, scale = TRUE)
names(pca)
pca$center
pca$rotation
head(pca$x)
biplot(x = pca, scale = 0, cex = 0.4, col = c("blue4", "brown3"))
'''


