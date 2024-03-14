k = 5
source('getXandV_X.R')
TV<- sum(V_X)
source('blockwisePCA_R_engine.R')
readRDS('reconstruction_results.rds')-> reclist
X_tables_mcl<- list(
                cnorm=reclist$X_norm,
                cmeta=reclist$X_meta,
                cnorm_noise=reclist$X_norm_noise,
                cmeta_noise=reclist$X_meta_noise)
#remove variables not modelled from metrics calculation:
	noise_mask <- apply(X_tables_mcl[[1]],2,sd)==0
	X<- X[,!noise_mask]
	lapply(X_tables_mcl, function(Xtab) Xtab[,!noise_mask] )-> X_tables_mcl

C_o<- cor(X)
message('computed C_o')
lapply( seq_along(X_tables_mcl), function(j) {
decomp<- X_tables_mcl[[j]]
if (j<3)
		
message('next')	
btwCORsq<- cor(X, decomp)^2
maxBTWcor_sq<- max(btwCORsq)
meanBTWcor_sq<- mean(btwCORsq)
gc(); rm(btwCORsq)
res=c(btwMaxCS= maxBTWcor_sq,
      btwMeanCS=meanBTWcor_sq,
  var_explained= .48,
  n_cl= NA,
  k=k, 
  noiseAdded=ifelse(j<3,0,1),
  total_nPC= NA,
  sqrtCOR_mse= sqrt(MSE( C_o, cor( decomp) ))
	)
	gc(); return(res)
		} )-> summarySIM




saveRDS(summarySIM,
	sprintf('sim_MCL_nPC%d_summary.rds',k))


