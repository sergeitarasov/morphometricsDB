
##X is training data
##Xp is test data
##cl: labels for the rows of X
##k = number of neighbours
##Returns average value of k nearest neighbours
fknn <- function(X,Xp,cl,k=1)
{
  out <- nabor::knn(X,Xp,k=k)
  cl[as.vector(out$nn.idx)] %>% matrix(dim(out$nn.idx)) %>% rowMeans
}

im <- load.image("/Users/taravser/Documents/My_papers/Dim_reduction/morphometricsDB/Helictopleurus/Head/data/hel.jpg")
plot(im)


grabRect(im, output = "coord")

#Coordinates of a foreground rectangle (x0,y0,x1,y1)
fg <- c(510,521,670,671 )
#Background
bg <- c(791,   28, 1020,  194 )
#Corresponding pixel sets
px.fg <- ((Xc(im) %inr% fg[c(1,3)]) & (Yc(im) %inr% fg[c(2,4)]))
px.bg <- ((Xc(im) %inr% bg[c(1,3)]) & (Yc(im) %inr% bg[c(2,4)]))

plot(im)
highlight(px.fg)
highlight(px.bg,col="blue")
