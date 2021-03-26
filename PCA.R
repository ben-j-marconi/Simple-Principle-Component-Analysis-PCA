setwd("~/Utah PhD/Projects/Range Creek Canyon/Data/XRF/RCC_CM_2018")


require(ggcorrplot)
require(GGally)
require(reshape2)
require(ggplot2)
require(ggfortify)

## Load this package 
require(factoextra)

##PCA
photon <- read.csv("RCC_CM_18_22.9_0-100_H_Y.csv", header=TRUE)


## Get rid of depth
photon <-(photon[,-1])

## Do this to make sure no missing (Na) is read and confuses the PCA
photon <-photon[which(complete.cases(photon)),]

head(photon,4)


photon.pca<-prcomp(photon, center = TRUE, scale = TRUE)

summary(photon.pca)



##autoplot-method ############################

autoplot(photon.pca, data = photon.pca,
         x = 3,
         y = 4,
         #xlim = c(-.025, .075),
         #ylim = c(-.1, .1),
        # linetype = TRUE,
        size = .001,
                colour = 'Black',
                loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)


########## This is the plots you want ######################

####### Factoextra method ##############################################################
## http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

pdf("test.pdf")


setwd("D:/XRF/PCA_plots")


#### This plots all 900 some rows
fviz_pca_ind(photon.pca,
             geom = c("point"),
             axes = c(1,2),
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("gray80", "gray60", "gray40", "gray10"),
             repel = FALSE, # Avoid text overlapping
             #geom = c("point")
)


########## This is the plots you want ######################
### This figures out dimensions by column/element  ######################

pdf("CoreB_PCA1_2pdf")
pdf("PCA1_3.pdf")

fviz_pca_var(photon.pca,
             axes = c(3,4),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("gray80", "gray50", "gray20", "black"),
             repel = TRUE,
             fill.var= "black",
             alpha.var = .75,
             col.circle = "black",
             title = "PCA1 & 2",
             #geom = c("arrow"),
             #geom = c("text"),
            
)
             #main = "Core B - PCA")

dev.off()

#### Biplot 

pdf("CoreA_PCA3_4_biplot.pdf")


fviz_pca_biplot(photon.pca, 
                repel = TRUE,
                axes = c(1,2),
                geom = c("point"),
                col.var = "#2E9FDF", # Variables color
                col.circle = "black",
                title = "Core A PCA 3 & 4",
                #col.ind = "#696969",  # Individuals color
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

dev.off()

### Access scores

# Eigenvalues
eig.val <- get_eigenvalue(photon.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(photon.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(photon.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation





# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <- photon.pca$rotation
sdev <- photon.pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:4])


# Compute Cos2
#::::::::::::::::::::::::::::::::::::::::
var.cos2 <- var.coord^2
head(var.cos2[, 1:4])



# Compute contributions
#::::::::::::::::::::::::::::::::::::::::
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])







##ggpairs method - See which variables are highly correlated 
plot1 = ggpairs(photon,
        columns=1:13)
        upper = list(continuous = "density"),
        lower = list(combo = "facetdensity"),
        title="B1He data")
print(plot1)

##Pairs method 
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

pairs(photon, upper.panel = panel.cor)

pairs(photon)

##actual PCA
photon.pca<-prcomp(photon, center = TRUE)
plot(photon.pca)
photon.pca
## photon.pca2<-princomp(photon, scores = TRUE, cor=TRUE)

## coefficients (eigenvectors- principle directions) of the linear combinations of the continuous variables.
print(photon.pca)

##summary stats of each PCA - how much variability in data is explained by each PCA 
summary(photon.pca)  ## PC1 or component one explains xxx% proportion of the variance in the dataset

##summary stats graph - shows the variances of each component (stdev)
plot(photon.pca, type="b")


print(photon.pca$x[,1:6])

##effect of each variable on each PCA
photon.pca$rotation[,1:5]

plot(photon.pca$rotation[,2:3])


photon.pca$x[,1:3]
sort(photon.pca$x[,1:4])

##PCA 1v2
pdf("test1.pdf")

biplot(photon.pca, xlim=c(-.05,.22)) ## Shows which elements are aligning on component 
       main="PCA1 & PCA2 of coreB He filter (Na-Fe)")

dev.off()
##PCA 2v3
biplot(photon.pca$x[,1:2])

##biplots w/o invidiual numbers, just category titles
pdf("biplot.pdf")
biplot(photon.pca, col=c(0,1),
       main="PCA1 & PCA2 of coreB He filter (Na-Fe)")
dev.off()

##individual PCAs vs depth of core

pdf("test4.pdf")

plot(photon$depth, photon.pca$x[,2], 
     xlab="Depth (mm)", 
     ylab="PCA Axis 1",
     main="PCA2 coreB He Filter (Na-Fe)")
     
dev.off()
##loess smoother for PCA vs depth of core graph
line1<-loess(photon.pca$x[,1]~photon$depth, data=photon, span=.2)



# THis is the code that makes the nice PCA with a box 
require(FactoMineR) 
require(missMDA)
# PCA with function PCA
setwd("C:/Users/Joshua/Desktop/Research_/Red Lab/XRF/CorAnalysis/B")


datMy <- read.csv("PCA_YF_B1.combined.csv", header = TRUE)
datMy <- read.csv("PCA_B.csv", header = TRUE)
p14 <-datMy[-1]

#read the tab file using the read table function.
pdf("test.PCAHe.pdf")

pca <- PCA(photo, scale.unit=TRUE,
                 ncp=6,
                 graph=T, axes=c(1,2))

plot(pca, habillage = "none", col.hab=c("red"))

dev.off()


#scale all the features,  ncp: number of dimensions kept in the results (by default 5)

pdf("test2.PCAHe.pdf")
pca <- imputePCA(p14, scale=TRUE, ncp=2)

plot(pca, habillage = "none", col.hab=c("green","blue"))

#This line of code will sort the variables the most linked to each PC. It is very useful when you have many variables.
dimdesc(pca)

## ggbiplot
library(devtools)
install.packages("digest")
install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggfortify)
library(cluster)
library(vqv)

print(ggbiplot(photon.pca, choices=1:2, obs.scale = 1, alpha = 0, var.scale = 1, ellipse = TRUE, circle = TRUE))

g <- ggbiplot(photon.pca, choices=3:4, scale = 1, obs.scale = .5, alpha=0, 
              var.scale = TRUE, groups = NULL, 
              varname.size = 7, labels.size=6,
              ellipse = TRUE, ellipse.prob = 0.99,
              circle = TRUE, var.axes=TRUE, varname.adjust=6)

dev.off()
              

g <- g + scale_color_discrete(name='')
g <- g + theme(text = element_text(size=20),
               legend.direction = 'horizontal',
               legend.position = 'top')


print(g)



#predict(pca,
        #newdata=tail(log.pca,2))

##Simple linear regressions
setwd("C:/Users/Joshua/Desktop/Research_/Red Lab/XRF/")


B1 <- read.csv("PCA_HE_B1.csv")


## Columns from dataset
Mg <-a[,3]
Al <-a[,4]
Si <-a[,5]
P <-a[,6]
S <-a[,7]
K <-a[,8]
Ca<-a[,9]
Ba <-a[,10]
Ti<-a[,11]
Cr <-a[,12]
Mn<-a[,13]
Fe<-a[,14]



##
cor.test(Si,Fe)



cor.test(c,e)


plot(c, e)


##Multiple linear regression - Lab 3 

summary(B1)

B1 <- B1[,-1]

mod.1 <- lm(Si ~ ., data=B1)
anova(mod.1)
summary(mod.1)


##Add row names somehow 
plot(mod.1)


##See if variables have colinearity
  #(i.e. correlations between the X variables)

cor(B1)


##Automatic variable selection 
  ##For datasets with multiple explanatory variables 
  ##Stepwise variable selection - adding or removing variables in an iterative way
  ##until the best model is reached.

mod.1 <- lm(Si ~ ., data=B1)
## Null model
mod.0 <- lm(Si ~ 1, data=B1)
step(mod.0, scope=formula(mod.1))

##Models with interaction varialbeS
Si.lm5 = lm(Si ~ K * Ti, data=B1)
summary(Si.lm5)

##Plots 
plot(Si ~ K, data=B1, col=1, pch=16)
coef(Si.lm5)

##Figure out how add abline - keep looking over old labs from Simons class
abline(coef(Si.lm5)[1], coef(Si.lm5)[2])
##Perform a Partial least square (PLS) to see which elements are related? 
##Perform factor analysis to see how variable are related 

require(ggplot)
require(ggfortify)

score <- which(photon(heptathlon) == "score")

