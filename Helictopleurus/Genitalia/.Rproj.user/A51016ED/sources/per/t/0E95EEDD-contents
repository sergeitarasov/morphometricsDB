
#---- Transformation Image into Procrutes coordinates
# Given:
# Im0 orignal Image  as imager object
# L0 original Landmraks for Im0 (cartesian normalized, y=y_org*-1)
# LP Procructes Landmark of L0
# beta0 the scaling parameter
# output.center the center of the output (cartesian normalized)
#
# Produces:
# Im1 aligned imager object
# LP1 Procructes Landmarks for Im1 (imager normalized)
#
# Algorithm
# T1: L0->LP->LP1,
# Im0*T1->Im1

# # cartesian normalized center
# output.center <- c(round(width(im1)/2, 0), -round(height(im1)/2, 0))
# # cartesian normalized L0
# L0 <- poiA.t[,,1]
# Im0 <- im1
# LP <- symproc$coords[,,1]
# beta0 <- 386
# plot(LP)
#
# out.ctr <- c(round(width(im1)/2, 0), -round(height(im1)/2, 0))
#
# im.aligned <- image2proc(Im0=im1, L0=poiA.t[,,1], LP=symproc$coords[,,1],
#                          beta0=386, output.center=out.ctr)
# plot(im.aligned$imageProc)
# points(im.aligned$coordLandm, col='red')

image2proc <- function(Im0, L0, LP, beta0, output.center){

  #----- 1. LP->LP1: take LP, scale and translate to center
  # scale
  LP1 <- LP*beta0
  # rotate, LP1 is cartesian normalized
  LP1 <- t(t(LP1)+output.center)

  # plot(Im0)
  # points(LP1 %*% M.coor, col='red')
  # text(LP1 %*% M.coor,labels=rownames(poi$gen1.txt), cex=0.6, pos=4, col="red")
  # plot(Im0)
  # points(L0 %*% M.coor, col='red')
  # text(L0 %*% M.coor,labels=rownames(poi$gen1.txt), cex=0.6, pos=4, col="red")

  #----- 2. find alignment transform from L0 to LP1
  rot <- rotonto(LP1, L0, scale = TRUE)

  # Given the alignment transfrom, transform the image

  # General way to deal with output of rotonto()
  # #t(t(L0) - rot$transy)
  # yrot <- sweep(L0, 2, rot$transy)
  # yrot <-rot$bet*yrot%*%rot$gamm
  # t(t(yrot) + rot$trans)
  # LP1
  #
  # we rearrange it as (scale+rotation) + M.trans:translation
  # M.trans <- -1*rot$bet*rot$transy%*%rot$gamm + rot$trans
  # L0.1 <-rot$bet*L0%*%rot$gamm
  # t(t(L0.1)+c(M.trans))
  # LP1

  #----- 3. find global translation M.trans
  M.trans <- -1*rot$bet*rot$transy%*%rot$gamm + rot$trans

  #----- 4. Transform image
  #-- 4.1 scale
  Im1 <- imresize(Im0, rot$bet, interpolation = 3) #%>% plot
  #plot(Im1)
  #points((rot$bet*L0) %*% M.coor, col='red')

  #-- 4.2 rotate
  x <- c(1,0)
  y <- x%*%rot$gamm
  ang <- matlib::angle(x, as.vector(y), degree = TRUE)
  ang.direction <- rot$gamm[2,1]/abs(rot$gamm[2,1])
  Im1 <-imrotate(Im1, ang*ang.direction, cx=0, cy=0) #%>% plot()
  #plot(Im1)
  #points((rot$bet*L0%*%rot$gamm) %*% M.coor, col='red')

  #-- 4. translate with M.trans
  # remeber to change direction of y by -1
  Im1 <- imshift(Im1, M.trans[1], -1*M.trans[2], boundary=1) # %>% plot
  # this matrix is reverse redirects y for Cartesian -> Imager coors
  M.coor <- matrix(c(1,0, 0,-1), 2, 2)
  #plot(Im1)
  #points(LP1 %*% M.coor, col='red')
  #points(t(t(L0.1)+c(M.trans)) %*% M.coor, col='red')

  list(imageProc=Im1, coordLandm=LP1 %*% M.coor)

}

# Miltiple alignment
image2procMulti <- function(Im0.list, L0.array, LP.array, beta0, output.center){

  n.im <- length(Im0.list)

  cat(paste0('Processing N(images)=', n.im, '\n'))
  cat(paste0(names(Im0.list), collapse = ' '), '\n')
  # n1=names(Im0.list)
  n2=dimnames(L0.array)[[3]]
  # n3=dimnames(LP.array)[[3]]

  IM <- imlist()
  CR <- c()

  i=1
  for (i in 1:n.im){
    cat(paste0('Processing Image ', i, '\n'))
    aligned1 <- image2proc(Im0=Im0.list[[i]], L0=L0.array[,,i], LP=LP.array[,,i],
                           beta0=beta0, output.center=output.center)

    IM[[i]] <- aligned1$imageProc
    CR <-abind(CR, aligned1$coordLandm, along=3)
  }


  names(IM) <- n2
  dimnames(CR)[[3]]<- n2

  out <- list(imageProc=IM, coordLandm=CR)
  return(out)
}



#--- converting between coordinates
# im2car.matrix(poiA[,,1], hei)
# im2car(poiA[,,1], hei)
# im2car(poiA[1,,1], hei)

im2car <- function(M, max){
  UseMethod("im2car", M)
}

im2car.matrix <- function(M, max){
  cbind(M[,1], - (M[,2]-max))
}

# M <- array(matrix(1:12), c(3,2,2))
# max=10
im2car.array <- function(M, max){
  #cbind(M[,1], - (M[,2]-max))
  dm <- dim(M)
  ar <- array(rbind(M[,1,], -(M[,2,]-max)), dim=dm)
  dimnames(ar)[3] <- dimnames(M)[3]
  return(ar)
}

im2car.numeric <- function(M, max){
  c(M[1], - (M[2]-max))
}

# Given array make Y coord negative
negateY <- function(ldk){
  M.coor <- matrix(c(1,0, 0,-1), 2, 2)
  ldk.new <- apply(ldk, 3, function(x) x %*% M.coor) %>% array(., dim(ldk))
  dimnames(ldk.new)[3] <- dimnames(ldk)[3]
  return(ldk.new)
}

#---

# pp <- poiA[1,,1]
# mt <- jpg.n
# get closest point to landmark
get_closest<-function(mt, pp){
  dis <- t(t(mt)-pp)^2
  dis <- sqrt(dis[,1]+dis[,2])
  which(dis==min(dis))[1]
}

# reoder Lm to make matrix starting with the fisrs LM
reoderLM <- function(ar, lm){
  p.id <- get_closest(ar, lm )
  ar1 <- ar[p.id:nrow(ar),]
  ar2 <- ar[-(p.id:nrow(ar)),]
  rbind(ar1,ar2)
}

# read outlines with ldk over one ldk
#readOutlines1(file.base, dirs, stereo.ldk, hei, reorder.ldk=1)
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


# read outlines with ldk over many ldk
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



# Sample Equidistant landmarks
# over one matrix
#outl <- dt$outl[[1]]
#ldk.sample <- digit.curves(outl[1,], outl, 100, closed = T)
#over list
#outl.list <- dt$outl
#semi.ldk <- digit.curves_List(dt$outl, nPoints=100)

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
