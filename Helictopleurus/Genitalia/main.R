# Morphometric workflow steps: --------------------------------------
#
# 1. Get body part images in JPG
# 2. Get body part outlines in JPG (most likely filled)
# 3. Use StereoMoprph to set key landmarks
# 4. Use key landmarks on outlines to sample equally distant semilandmarks on outlines
# 4. Use semilandmarks for Procrustes
# 5. Use outlines for ellyptic Fourier


#   ____________________________________________________________________________
#   Morphometric Workflow                                                   ####

# Federica edits this file

project.dir <- getwd()
source('R/dependencies.R')


#   ____________________________________________________________________________
#   Adobe Illustrator for processing images                                 ####

rstudioapi::navigateToFile("Ai/README_Ai.txt")


#   ____________________________________________________________________________
#   StereoMorph: Digitize Landmarks                                         ####

# Go to the respective script
rstudioapi::navigateToFile("R/stereoMorph.R")


#   ____________________________________________________________________________
#   Process ldks and other data                                             ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Parameters to adjust                                                    ####


stereo.ldk <- readRDS('data/stereo_ldk/stereo.ldk.RDS') # read stereo landmarks

file.list <- strsplit(names(stereo.ldk), '_im.txt') %>% unlist  # files to go through obtained from stereo obj
reorder.ldk <- 1                                                # reorder outlines given first ldk
imh <- load.image('data/im/viridans_im.JPG')                    # image height, use the first image
hei=height(imh)                                                 # image height, use the first image

dirs <- list(
  file.im= '_im.JPG',    # image files extension
  file.out='_out.JPG',   # outline files extension
  file.stereo='_im.txt', # stereo files extension
  dir.im='im',           # image folder in 'path'
  dir.out='out',         # outline folder in 'path'
  path = 'data'          # path to data
)
####


##  ............................................................................
##  Read outlines and landmarks from a bunch of files                       ####

dt <- readOutlines2(file.list, dirs, stereo.ldk, hei, reorder.ldk=1)
str(dt)

# plot 1st of them
plot(dt$outl[[1]], type='l')
points(dt$ldk[[1]] , col='red')

# save dt
saveRDS(dt, 'data/saved/dt.RDS')

##  ............................................................................
##  Sample equdistant semilandmarks along outlines                          ####

# the first ldk in each outline matrix is used as a starting point
semi.ldk0 <- digit.curves_List(dt$outl, nPoints=200)
str(semi.ldk0)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Read file HERE                                                          ####
# save
saveRDS(semi.ldk0, 'data/saved/semi.ldk0.RDS')
# read
semi.ldk0 <- readRDS('data/saved/semi.ldk0.RDS')


# translate to Momocs object, outline + 1st ldk !!!! ERORS
# semi.out <- Momocs::Out(semi.ldk0, ldk=as.list(rep(1,dim(semi.ldk0)[3])))
# names(semi.out) <- dimnames(semi.ldk0)[[3]]
# semi.out

# Error in coo_check.Coo(Coo) : 13, 28 do not pass coo_check

semi.ldk01 <- semi.ldk0[,,-c(13,28)]
semi.out <- Momocs::Out(semi.ldk01, ldk=as.list(rep(1,dim(semi.ldk01)[3])))
names(semi.out) <- dimnames(semi.ldk01)[[3]]
semi.out

#--------
# we will need to repair!!!!!!!!!!
# Error in coo_check.Coo(Coo) : 13, 28 do not pass coo_check
dim(semi.ldk0)
semi.ldk0[,,13]

#------------
saveRDS(semi.out, 'data/saved/semi.out.RDS')

# Plot outlines and ldk
semi.out %>% panel(names=TRUE, points=T, points.cex = 10, points.pch = 9)
inspect(semi.out)

semi.out %>%
  coo_center %>% coo_scale %>%
  coo_alignxax() %>% coo_slidedirection("up") %T>%
  print() %>% stack()


##  ............................................................................
##  Procrustes superimposition                                              ####

# semi.ldk0 # original sampled ldk array
semi.out  # Momocs outline obj
semi.out.ldk <- semi.out %$% l2a(coo) # Momocs outline obj -> back to ldk array

#semi.ldk <- Ldk(coo=semi.ldk0 %>% a2l) # Momocs ldk
semi.ldk <- Ldk(coo=semi.ldk01 %>% a2l)

# Procrustes superimposition using different methods
# using package geomorph
gpag <- gpagen(semi.out.ldk, ProcD = F, max.iter = 50)
# using package  Momocs
fgProc <- fgProcrustes(semi.ldk)
# using package Morpho
procSy <- procSym(semi.out.ldk, reflect = FALSE)


# Plot

plot(gpag)
Ldk(coo=gpag$coords %>% a2l) %>% panel()
Ldk(coo=gpag$coords %>% a2l) %>% stack()

fgProc %>% Out %>% stack()
fgProc %>% panel()

Ldk(coo=procSy$rotated %>% a2l) %>% Out %>% stack()
Ldk(coo=procSy$rotated %>% a2l) %>% panel()


##  ............................................................................
##  Procrustes using sliders                                                ####

sliders = define.sliders(c(1:100,1), write.file = FALSE)
gpag.s <- gpagen(semi.out.ldk, ProcD = F, curves = sliders)

plot(gpag.s)
Ldk(coo=gpag.s$coords %>% a2l) %>% panel()
Ldk(coo=gpag.s$coords %>% a2l) %>% stack()
Ldk(coo=gpag.s$coords %>% a2l) %>% Out %>% stack()

##  ............................................................................
##  Elliptic Fourier                                                        ####

semi.out %>% stack()

#-- normalization
f1 <- efourier(semi.out, 10, norm=T)
PCA(f1) %>% plot(morpho=T)
as_df(f1)

#-- manual alognment
semi.out.al <- semi.out %>%
  coo_center %>% coo_scale %>%
  coo_alignxax() %>% coo_slidedirection("up")

semi.out.al %>% stack()
f2 <- efourier(semi.out.al, 10, norm=F)
PCA(f2) %>% plot(morpho=T)

#-- Procrustes alighment
fgProc %>% Out %>% stack
f3 <- efourier(fgProc %>% Out, 10, norm=F)
PCA(f3) %>% plot(morpho=T)

#-- calliper alighment
s4 <- coo_aligncalliper(semi.out)
s4 %>% stack
f4 <- efourier(s4, 10, norm=F)
PCA(f4) %>% plot(morpho=T)



#   ____________________________________________________________________________
#   Dim Estimation and Phylo Inference                                      ####

rstudioapi::navigateToFile("R/Dim_est.R")

#   ____________________________________________________________________________
#   Image Alignment                                                         ####

rstudioapi::navigateToFile("R/Image2Proc.R")



# Helpful Functions -----
# methods(class=Out)
# shp <- shapes[4]
# apropos("coo_")
# coo_plot(shp)
# coo_plot(shp, col="grey80", border=NA, centroid=FALSE, main="Meow")
# coo_plot(coo_center(shp), main="centered Meow")
# coo_plot(coo_sample(shp, 64), points=TRUE, pch=20, main="64-pts Meow")
# shapes[4] %>% coo_smooth(5) %>% coo_sample(64) %>% coo_scale() %>% coo_plot()
#
###
