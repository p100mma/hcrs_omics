
X<- readRDS('~/ECAI_computations/rdsFiles/gene_expr_data.rds')
means_X<- colMeans(X)
svd(scale(X, center=TRUE, scale=FALSE) )->svd_X
PC= svd_X$u[,1:70] %*% diag(svd_X$d[1:70])
X_rec<- PC %*% t(svd_X$v[,1:70])
apply(X,2,var)-> V_X
apply(X_rec,2,var)-> V_rec
V_left = V_X - V_rec
X_rec_noise<- X_rec
set.seed(43124)
for (j in 1:ncol(X_rec_noise)) {
	X_rec_noise[,j]= X_rec_noise[,j] + rnorm(nrow(X_rec_noise), 0, sqrt( V_left[[j]] ) )
	X_rec_noise[,j]= X_rec_noise[,j] + means_X[[j]]
	X_rec[,j]= X_rec[,j] + means_X[[j]]
	}
saveRDS( list(PC=PC,
	      X_rec=X_rec,
	      X_rec_noise=X_rec_noise), 'plainSVDreconstruct70.rds')

PC= svd_X$u[,1:10] %*% diag(svd_X$d[1:10])
X_rec<- PC %*% t(svd_X$v[,1:10])
apply(X,2,var)-> V_X
apply(X_rec,2,var)-> V_rec
V_left = V_X - V_rec
X_rec_noise<- X_rec
set.seed(43184)
for (j in 1:ncol(X_rec_noise)) {
	X_rec_noise[,j]= X_rec_noise[,j] + rnorm(nrow(X_rec_noise), 0, sqrt( V_left[[j]] ) )
	X_rec_noise[,j]= X_rec_noise[,j] + means_X[[j]]
	X_rec[,j]= X_rec[,j] + means_X[[j]]
	}
saveRDS( list(PC=PC,
	      X_rec=X_rec,
	      X_rec_noise=X_rec_noise), 'plainSVDreconstruct10.rds')

PC= svd_X$u %*% diag(svd_X$d)
PCvars<- apply(PC,2,var)
varTot<- sum(PCvars)
varCumul<- cumsum(PCvars)
perc48= min(which( varCumul >= 0.48*varTot))
X_rec<- PC[,1:perc48] %*% t(svd_X$v[,1:perc48])
apply(X,2,var)-> V_X
apply(X_rec,2,var)-> V_rec
V_left = V_X - V_rec
X_rec_noise<- X_rec
set.seed(23423)
for (j in 1:ncol(X_rec_noise)) {
	X_rec_noise[,j]= X_rec_noise[,j] + rnorm(nrow(X_rec_noise), 0, sqrt( V_left[[j]] ) )
	X_rec_noise[,j]= X_rec_noise[,j] + means_X[[j]]
	X_rec[,j]= X_rec[,j] + means_X[[j]]
	}
saveRDS( list(PC=PC,
	      X_rec=X_rec,
	      X_rec_noise=X_rec_noise), 'plainSVDreconstruct48perc.rds')

