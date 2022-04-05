
#-- calliper alighment
s4 <- coo_aligncalliper(semi.out)
s4 %>% stack
f4 <- efourier(s4, 40, norm=F)
PCA(f4) %>% plot(morpho=T)
xyz <- as_df(f4)
#-

ump <- umap(xyz, n_components = 2, n_neighbors = 2, learning_rate = 0.5, init = "spectral", n_epochs = 20)
plot(ump)
ump %>% as_tibble()

#----
iris30 <- iris[c(1:10, 51:60, 101:110), ]
iris_umap <- umap(iris30, n_neighbors = 5, learning_rate = 0.5, init = "random", n_epochs = 20)
plot(iris_umap)

umap.out <- umap(dt, n_neighbors = 5, bandwidth = 1, min_dist=.1, spread = 1, learning_rate = 0.5, init = "spectral", n_components = 2)
umap.out <-as_tibble(umap.out)
umap.out <-mutate(umap.out, label=nn,  f.id=pars$f.id, k.id=pars$k.id)

ggplot(data = umap.out, aes(x = V1, y = V2) )+
  geom_point(aes(colour = as.character(f.id), size=2 )) +
  geom_text(aes(label=label),hjust=.5, vjust=1)


#--------
mle_twonn(X = xyz)
linfit_twonn(X = xyz, trimmed = F, alpha_trimmed = 0)

#------
Swissroll <- as_tibble(Swissroll_maker(N = 1000))
mus_Swiss <- generate_mus(X = Swissroll)

lin11 <- linfit_twonn(X = Swissroll, trimmed = F, alpha_trimmed = 0)
lin12 <- linfit_twonn(X = Swissroll, trimmed = T, alpha_trimmed = 0.01)


dist_Eucl_Sw <- as.matrix(stats::dist(Swissroll))
mle_twonn(X = Swissroll)
mle_twonn(X = Swissroll, DM=dist_Eucl_Sw)
mle_twonn(X = Swissroll[c(1:3),])
