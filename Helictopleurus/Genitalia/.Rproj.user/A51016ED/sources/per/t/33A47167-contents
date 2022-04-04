
# Go to the main script to get objects if necessary
rstudioapi::navigateToFile("main.R")


#   ____________________________________________________________________________
#   Pairwise image alignment                                               ####

#---------- parameters:

file.list  # files to go through obtained from stereo obj
dirs # directories see main.R
gpag # Procrustes ldk #gpag <- gpagen(semi.out.ldk, ProcD = F, max.iter = 50)
semi.out.ldk # original postprocessed outlines see main.R

file.base1 <- file.list[1]
file.base2 <- file.list[2]
###

file1 <- file.path(dirs[['path']], dirs[['dir.im']], paste0(file.base1, dirs[['file.im']]) )
file2 <- file.path(dirs[['path']], dirs[['dir.im']], paste0(file.base2, dirs[['file.im']]) )

im1 <- load.image(file1)
im2 <- load.image(file2)
hei=height(im1)

plot(im1, main = file1)
points(semi.out.ldk[,,file.base1])
plot(im2, main = file2)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Transform original outlines                                             ####

# for image2proc() work:
# 1. plotting Im0 and L0 with imager returns Ok results
# 2. Then transform L0 by negateY() and use it  image2proc()

# condition 1 holds by tranfoming semi.out.ldk
plot(im1)
im2car(semi.out.ldk[,,file.base1], hei) %>% points

# Transform original ldk to make L0 given cond. 2
L0.all <- ldk.new <- semi.out.ldk %>% im2car(., hei) %>% negateY
###


##  ............................................................................
##  Pairwise Alignment                                                       ####

# Transformation of Image into Procrutes coordinates
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

out.ctr <- c(round(width(im1)/2, 0), -round(height(im1)/2, 0)) # cartesian normalized center
beta0 <- 4000 # scaler

# Image 1
im.aligned <- image2proc(Im0=im1, L0=L0.all[,,file.base1], LP=gpag$coords[,,file.base1],
                         beta0=beta0, output.center=out.ctr)
plot(im.aligned$imageProc)
points(im.aligned$coordLandm, col='red')
str(im.aligned)

# Image 2
im.aligned2 <- image2proc(Im0=im2, L0=L0.all[,,file.base2], LP=gpag$coords[,,file.base2],
                         beta0=beta0, output.center=out.ctr)
plot(im.aligned2$imageProc)
points(im.aligned2$coordLandm, col='red')


# Plot together Im1, Im2
imdraw(im.aligned$imageProc, im.aligned2$imageProc, opacity = .3) %>% plot
lines(im.aligned$coordLandm, col='red')
lines(im.aligned2$coordLandm, col='blue')



#   ____________________________________________________________________________
#   Multiple Alignment                                                      ####

#---------- parameters:

file.list  # files to go through obtained from stereo obj
dirs # directories see main.R
gpag # Procrustes ldk #gpag <- gpagen(semi.out.ldk, ProcD = F, max.iter = 50)
semi.out.ldk # original postprocessed outlines see main.R

# Transform original ldk to make L0 given compatible for alignment
L0.all <- ldk.new <- semi.out.ldk %>% im2car(., hei) %>% negateY

# read all images in path
files <- file.path(dirs[['path']], dirs[['dir.im']])
im.list <- load.dir(path=files, pattern = NULL, quiet = FALSE)
str(im.list)

im1 <- im.list[[1]]
out.ctr <- c(round(width(im1)/2, 0), -round(height(im1)/2, 0)) # cartesian normalized center
beta0 <- 4000 # scaler

###




# multiple alignment
im.ali <- image2procMulti(Im0.list=im.list, L0.array=L0.all,
                          LP.array=gpag$coords, beta0=beta0, output.center=out.ctr)
#saveRDS(im.ali, 'data/saved/im.ali.RDS')

# Plot
i=4
plot(im.ali$imageProc[[i]])
lines(im.ali$coordLandm[,,i], col='red')

# Plot superimpose
imdraw(im.ali$imageProc[[1]], im.ali$imageProc[[2]], opacity = .3) %>% plot
lines(im.ali$coordLandm[,,1] %>% as.matrix(), col='red')
lines(im.ali$coordLandm[,,2] %>% as.matrix(), col='blue')

