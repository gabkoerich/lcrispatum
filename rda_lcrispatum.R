# RDA Analysis 
# Plot based in https://rpubs.com/seblun/rda

setwd("C:/Users/gabri/Google Drive/Mestrado/Experimento/Resultados")

val <- read.csv("RDA_val_pos_test.csv", row.names = 1)
trat <- read.csv("RDA_trat.csv", row.names = 1)

library(vegan)
# Check if matrices have the same number of observations, or lines 
dim(val) #response variables
dim(trat) #explanatory variables
nrow(val) - nrow(trat)

val.chi <- decostand(val, "chi.square")

val.rda <- rda(val.chi ~ ., data=trat) # RDA constrained
val.rda

coef(val.rda)
summary(val.rda)
(R2 <- RsquareAdj(val.rda)$r.squared) 
(R2adj <- RsquareAdj(val.rda)$adj.r.squared) # R2-adjusted

# global test
anova.cca(val.rda, permutations=how(nperm=999)) 
# axis test (by="axis")
anova.cca(val.rda, by="axis", permutations=how(nperm=999)) 

# check if correlations between variables as below 3
vif.cca(val.rda)
#ok 

# Importance of components:
#   RDA1    RDA2     RDA3      PC1       PC2       PC3       PC4
# Eigenvalue            0.08339 0.01483 0.001176 0.008892 0.0001271 7.759e-05 4.865e-05
# Proportion Explained  0.76828 0.13663 0.010832 0.081922 0.0011711 7.149e-04 4.482e-04
# Cumulative Proportion 0.76828 0.90491 0.915744 0.997666 0.9988369 9.996e-01 1.000e+00


tiff("rda.tiff", width=250, height=250, bg="white", res=300, unit="mm")

par(mar=c(4,4,2,2))
plot(val.rda, scaling=2, display=c("cn", "lc", "sp"), type="n", 
     xlab="RDA1 (76.83 %)", ylab="RDA2 (13.66 %)", 
     xlim=c(-1,1.6), ylim=c(-1,1), cex.axis=1.5, cex.lab=1.5)

# Points for values
sites <- scores(val.rda, choices=1:2, scaling=2, display="lc")
points(sites, pch=1, cex=2)

# Points for variables
sp <- scores(val.rda, choices=1:2, scaling=2, display="sp")
points(sp, pch=17, cex=2, col="gray25")

# Arrows
env <- scores(val.rda, choices=1:2, scaling=2, display="bp")
arrows(0,0, env[1:3,1], env[1:3,2], lty=1, lwd=2.5, length=0.1, col="gray25")
env.names <- c("Temperature", "Time", "Nutrients + pCO2")
text(env[1:2,], env.names[1:2], cex=2, font=2, pos=4, col="gray25")
text(x = 1.72, y=0.09, env.names[3], cex=2, font=2, pos=2, col="gray25")

# Bottom text
text(0.9, -1.25, labels="F:94.195  R²:0.92  p:0.001", cex=2, col="gray25")

# Other labels added in Inkscape

dev.off()
