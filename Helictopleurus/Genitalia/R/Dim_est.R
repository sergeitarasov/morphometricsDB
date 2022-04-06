
#   ____________________________________________________________________________
#   Dim Estimation and Phylo Inference                                      ####

# back to main
rstudioapi::navigateToFile("main.R")

##  ............................................................................
##  Get Fourier XY data                                                     ####

#-- calliper alighment
s4 <- coo_aligncalliper(semi.out)
s4 %>% stack
f4 <- efourier(s4, 10, norm=F)
PCA(f4) %>% plot(morpho=T)

xyz <- as_df(f4)
spp <- names(s4)
str(xyz)


##  ............................................................................
##  Estimate Dimension                                                      ####

dim.mle <- mle_twonn(X = xyz)
dim.lst <- linfit_twonn(X = xyz, trimmed = F, alpha_trimmed = 0)
dim.mle
dim.lst

dim.est <- round(dim.mle[2],0)


#   ____________________________________________________________________________
#   UMAP Projection                                                         ####

#ump <- umap(xyz, n_components = 3, n_neighbors =2, learning_rate = 0.5, bandwidth = 1, init = "spectral", n_epochs = 20)
ump <- umap(xyz, n_components = 2, scale = T, n_neighbors =2, learning_rate = 0.5, bandwidth = 1, init = "random", n_epochs = 20)
rownames(ump) <- names(s4)
plot(ump)
ump



#   ____________________________________________________________________________
#   Infer tree                                                             ####
# BM tree inference using Phylip as library(Rphylip)

# set path to Phylip executables
setPath(path='/Users/taravser/Documents/Soft/phylip-3.695/exe')
#clearPath()

rm(phy)
phy <- Rcontml(ump, quiet = T)
plot(phy)
str(phy)



#   ____________________________________________________________________________
#   PCA                                                                     ####


f1 <- efourier(semi.out, 10, norm=T)
pca <- PCA(f1)
pca %>% plot(morpho=T)
str(pca)
pca.mt <- pca$x[,1:2]
plot(pca.mt)

# infer tree
rm(phy)
phylp <- Rcontml(pca.mt, path='/Users/taravser/Documents/Soft/phylip-3.695/exe')
plot(phy)

