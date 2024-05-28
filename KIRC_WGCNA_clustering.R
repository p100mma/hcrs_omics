    tom<- function(A) {diag(A)=0; deg<-colSums(A);  (A + A%*%A)/( outer(deg,deg,pmin) +1 - A )  }

    dtom_hcl<- function(TOM,...) { dt= 1- TOM; message('done dt'); diag(dt)=0; rm(TOM); hclust(as.dist(dt), ...) }

    plot_dth<- function(dth)   plot(dth, labels = FALSE, check=FALSE)

source('KIRC_X_andV_X.R')

TM= tom(abs(cor(X))^5)
message('done TOM')
gc()
dtom_hcl(TM,'average')-> dhcl
message('done dhcl')
gc()
save.image('KIRC_WGCNA_clustering.RData')
message('saved')
