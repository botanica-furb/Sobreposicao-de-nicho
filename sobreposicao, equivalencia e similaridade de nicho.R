install.packages("ecospat")
library(ecospat)
# Importe os valores ambientais de cada espécie/população
sp <- exemplo.sobreposicao

pca.env <- dudi.pca(rbind(sp)[,4:10],scannf=F,nf=2)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

# PCA scores for the whole study area
scores.globclim <- pca.env$li
# PCA scores for the sp1 distribution
scores.sp.sp <- suprow(pca.env,sp[which(sp[,17]==1),4:10])$li
# PCA scores for the species sp2 distribution
scores.sp.sp2 <- suprow(pca.env,sp[which(sp[,18]==1),4:10])$li
# PCA scores for the whole sp1 study area
scores.clim.sp <- suprow(pca.env,sp[,4:10])$li
# PCA scores for the whole sp2 study area
scores.clim.sp2 <- suprow(pca.env,sp[,4:10])$li


grid.clim.sp <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                      glob1=scores.clim.sp,
                                      sp=scores.sp.sp, R=500,
                                      th.sp=0)

grid.clim.sp2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.sp2,
                                       sp=scores.sp.sp2, R=500,
                                       th.sp=0)

D.overlap <- ecospat.niche.overlap (grid.clim.sp, grid.clim.sp2, cor=T)$D


D.overlap

eq.test <- ecospat.niche.equivalency.test(grid.clim.sp, grid.clim.sp2,
                                          rep=500, alternative = "greater")

eq.test$p.D        # p.D < 0.05 = Equivalência, p.D > 0.05 = Não equivalente

sim.test <- ecospat.niche.similarity.test(grid.clim.sp, grid.clim.sp2,
                                                  rep=500, alternative = "greater",
                                                  rand.type=1)

sim.test$p.D       #p.D < 0.05 = Similaridade, p.D > 0.05 não similar


ecospat.plot.overlap.test(eq.test, "D", "Equivalency")

ecospat.plot.overlap.test(sim.test, "D", "Similarity")

niche.dyn <- ecospat.niche.dyn.index (grid.clim.sp, grid.clim.sp2, intersection = 0.1)

ecospat.plot.niche.dyn(grid.clim.sp, grid.clim.sp2, quant=0.25, interest=2,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")


ecospat.plot.niche(grid.clim.sp)
ecospat.plot.niche(grid.clim.sp2)

ecospat.shift.centroids(scores.sp.sp, scores.sp.sp2, scores.clim.sp, scores.clim.sp2)


