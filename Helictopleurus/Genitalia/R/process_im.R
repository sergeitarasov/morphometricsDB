# To do
# 1. scale images to W=2K


setwd("~/Documents/My_papers/Dim_reduction/Stereomorph")
library(imager, shapes)
library(StereoMorph)
library(ggplot2)
library(matlib)
library(Momocs)
library(geomorph)
library(Morpho)
library(abind)
source('ImportSM.R')
source('Parachorius/Functions_morpho.R')

setwd("~/Documents/My_papers/Dim_reduction/Stereomorph/Helictopleurus")

file <- 'data/im/viridans_im.JPG'
im1 <- load.image(file)
plot(im1)


file <- 'data/out/viridans_out.JPG'
jpg <- import_jpg1(file, auto.notcentered =T)
jpg
plot(jpg, type='l')


# To imager coords
plot(im1)
hei=height(im1)
jpg.n <- im2car(jpg, hei) # to imager coord
#jpg.n <-reoderLM(jpg.n, poiA[1,,1]) # reoder given firsr Lm
lines(jpg.n, col='red')



#_----------------------

# read stereo landmarks
stereo.gen <- readRDS('data/stereo_gen.RDS')

stereo.ldk <- stereo.gen
file.list <- strsplit(names(stereo.ldk), '_im.txt') %>% unlist # files to go through obtained from stereo obj
reorder.ldk <- 1 # reoder outline given first ldk

file.base <- 'unifasciatus'
hei=height(im1) # image height


dirs <- list(
  file.im= '_im.JPG',
  file.out='_out.JPG',
  file.stereo='_im.txt',
  dir.im='im',
  dir.out='out',
  path = 'data'
  )



# read outlines with ldk
readOutlines1(file.base, dirs, stereo.ldk, hei, reorder.ldk=1)
readOutlines1 <- function(file.base, dirs, stereo.ldk, hei, reorder.ldk=1){

  # get ldk from stereo obj
  file <- paste0(file.base, dirs[['file.stereo']])
  ldk <- stereo.ldk[[file]] %>% as.matrix()
  ldk <-im2car(ldk, hei) # to the same coors as outline

  file <- file.path(dirs[['path']], dirs[['dir.out']], paste0(file.base, dirs[['file.out']]) )
  outl <- import_jpg1(file, auto.notcentered =T)
  outl <- reoderLM(outl , ldk[reorder.ldk,]) # reoder given firsr Ldk

  # change direction of the curve points to be clockwise
  del <- outl[1,1]-outl[2,1] # if del is negative: anticlockwise, if pos: clockwise
  # change direction if del anticlockwise
  if (del<0){
    tt <- outl[2:nrow(outl),]
    tt <-cbind(rev(tt[,1]), rev(tt[,2]))
    outl <- rbind(outl[1,], tt)
  }

  list(file.base=file.base, ldk=ldk, outl=outl)
  #plot(outl, type='l')
  #points(ldk , col='red')
}



readOutlines2 <- function(file.list, dirs, stereo.ldk, hei, reorder.ldk=1){

  file.base.Out <- c()
  ldk.Out <- list()
  outl.Out<- list()

  #i=1
  for (i in 1:length(file.list)){
    flb=file.list[i]
    cat('Reading outlines:', flb, '\n')
    readOut1 <- readOutlines1(file.base=flb, dirs, stereo.ldk, hei, reorder.ldk=1)
    #plot(readOut1$outl, type='l')
    #points(readOut1$ldk , col='red')

    file.base.Out <-c(file.base.Out, readOut1$file.base)
    #ldk.Out <- abind(ldk.Out, readOut1$ldk, along=3)
    #outl.Out<-abind(outl.Out, readOut1$outl, along=3)
    ldk.Out[[i]] <- readOut1$ldk
    outl.Out[[i]]<-readOut1$outl
  }

  names(ldk.Out) <- file.base.Out
  names(outl.Out) <- file.base.Out

  ddd <- list(files=file.base.Out, ldk=ldk.Out, outl=outl.Out)
  #str(ddd)
  return(ddd)

}

dt <- readOutlines2(file.list, dirs, stereo.ldk, hei, reorder.ldk=1)
# plot read outlines
str(dt)
plot(dt$outl[[2]], type='l')
points(dt$ldk[[2]] , col='red')

# Sample Equidistant landmarks
# over one matrix
outl <- dt$outl[[1]]
ldk.sample <- digit.curves(outl[1,], outl, 100, closed = T)
plot(ldk.sample, type='l')
plot(ldk.sample)
points(dt$ldk[[1]] , col='red')

#over list
outl.list <- dt$outl
semi.ldk <- digit.curves_List(dt$outl, nPoints=100)
str(semi.ldk)

digit.curves_List <- function(outl.list, nPoints=100){

  ldk.Out <- c()
  i=1
  for (i in 1:length(outl.list)){
    curve <- outl.list[[i]]
    ldk.sample <- digit.curves(curve[1,], curve, nPoints, closed = T)
    ldk.Out <- abind(ldk.Out, ldk.sample, along=3)
  }
  dimnames(ldk.Out)[[3]] <- names(outl.list)
  ldk.Out
}

#----- Procrustes alighment

semi.ldk <- digit.curves_List(dt$outl, nPoints=150)
semi.out <- Out(semi.ldk)
semi.out %>%
  coo_center %>% coo_scale %>%
  coo_alignxax() %>% coo_slidedirection("up") %T>%
  print() %>% stack()
semi.out %>% panel

coo_plot(semi.out[3])
plot(semi.out[3])
points(semi.out[3][1,1], semi.out[3][1,2], col='blue')
points(semi.out[3][2,1], semi.out[3][2,2], col='red')
semi.out[3][1,]
semi.out[3][2,]

coo_plot(semi.out[2])
Ldk(coo=semi.ldk %>% a2l) %>% panel()


Ldk(coo=semi.ldk %>% a2l) %>%
  coo_center %>% coo_scale %>%
  coo_alignxax() %>% coo_slidedirection("up") %T>%
  print() %>% stack()

ll <- Ldk(coo=semi.ldk %>% a2l) %>%
  coo_center %>% coo_scale %>%
  coo_alignxax() #%>% stack()
ll %>% stack()

sliders = define.sliders(c(1:150,1))
symproc <- gpagen(semi.ldk[,,c(2:3)], ProcD = F, curves = sliders)
symproc <- gpagen(semi.ldk[,,c(1:4)], ProcD = F, max.iter = 50)

# PLot
plot(symproc)
Ldk(coo=symproc$coords %>% a2l) %>% panel()
Ldk(coo=symproc$coords %>% a2l) %>% stack()
coo_plot(semi.ldk[,,3])

symproc2 <- fgProcrustes(semi.ldk[,,c(1:4)])
symproc2 <- fgProcrustes(ll)
symproc2 %>%  str()
Ldk(coo=symproc2$rotated %>% a2l) %>% stack()

symproc3 <- procSym(semi.ldk, reflect = FALSE)
Ldk(coo=symproc3$rotated %>% a2l) %>% stack()

Ldk(coo=semi.ldk[,,c(3,3)] %>% a2l) %>%
  coo_center %>% coo_scale %>%
  coo_alignxax() %>% stack()

#--- Outlines in Momocs
ol <- Out(jpg)
plot(ol[1])

coo_plot(ol[1])
ol.f <- efourier(semi.out, 10)
boxplot(ol.f, drop=1)

bot.p <- PCA(ol.f)
summary(bot.p)
  class(bot.p)        # a PCA object, let's plot it
plot(bot.p)

ol.i <- efourier_i(ol.f, nb.h=25)
coo_draw(ol.i, border='red', col=NA)
#----

#   ____________________________________________________________________________
#   img                                                                     ####

# read image
file <- file.path(path, dirs[['dir.im']], paste0(file.base, dirs[['file.im']]) )
im1 <- load.image(file)

#plot
plot(im1)
#jpg.n <-reoderLM(jpg.n, poiA[1,,1]) # reoder given firsr Lm
lines(im2car(outl, hei), col='red')
points(im2car(ldk, hei), col='blue')
