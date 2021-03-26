#Principle Component Analysis for DImensional Decomposition 
setwd("~/Utah PhD/Projects/Range Creek Canyon/Data/XRF/RCC_CM_2018/Stacked_22.9_75_95")

require(ggcorrplot)
require(GGally)
require(reshape2)
require(ggplot2)
require(ggfortify)
require(factoextra)

#Data for PCA
photon <- read.csv("xrf_data.csv", header=TRUE)

#Remove Depth
photon <-(photon[,-1])

#Read no missing XRF Data
photon <-photon[which(complete.cases(photon)),]
head(photon,4)
photon.pca<-prcomp(photon, center = TRUE, scale = TRUE)
summary(photon.pca)

#Plots
#AutoPlot
autoplot(photon.pca, data = photon.pca,
         x = 1,
         y = 2,
         #xlim = c(-.025, .075),
         #ylim = c(-.1, .1),
         # linetype = TRUE,
         size = .001,
         colour = 'Black',
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5)

#Circular Greyscale Plot
fviz_pca_var(photon.pca,
             axes = c(2,3),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("gray80", "gray50", "gray20", "black"),
             repel = TRUE,
             fill.var= "black",
             alpha.var = .75,
             col.circle = "black",
             title = "Stacked Cores: PCA 2 & 3",
             #geom = c("arrow"),
             #geom = c("text"),
             
)
#main = "Core B - PCA")

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

# Eigenvalues
eig.val <- get_eigenvalue(photon.pca)
eig.val
