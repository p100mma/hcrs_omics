args<- commandArgs(trailingOnly=TRUE)
k = as.integer(args[[1]])
source('getXandV_X.R')
source('blockwisePCA_R_engine.R')
gs<- readRDS(sprintf('vary_nPC/decomposition_MCL_nPC%d.rds',k))
#removing C_0
	noise_mask <- apply(gs[[1]]$rec,2,sd)==0
	X<- X[,!noise_mask]

TV<- sum(V_X)
C_o<- cor(X)

lapply( gs, function(decomp) {
message('next')	
btwCORsq<- cor(X, decomp$rec_noise[,!noise_mask])^2
maxBTWcor_sq<- max(btwCORsq)
meanBTWcor_sq<- mean(btwCORsq)
gc(); rm(btwCORsq)
res=c(btwMaxCS= maxBTWcor_sq,
      btwMeanCS=meanBTWcor_sq,
  var_explained= sum( apply(decomp$rec,2,var) ) / TV,
  n_cl= length(decomp$hclustering_list),
  k=k, 
  noiseAdded=1,
  total_nPC= k*length(decomp$hclustering_list),
  sqrtCOR_mse= sqrt(MSE( C_o, cor( decomp$rec_noise[,!noise_mask]) ))
	)
	gc(); return(res)
		} )-> summaryPerLVL_noise


lapply( gs, function(decomp) {
message('next')	
btwCORsq<- cor(X, decomp$rec[,!noise_mask])^2
maxBTWcor_sq<- max(btwCORsq)
meanBTWcor_sq<- mean(btwCORsq)
gc(); rm(btwCORsq)
res=c(btwMaxCS= maxBTWcor_sq,
      btwMeanCS=meanBTWcor_sq,
  var_explained= sum( apply(decomp$rec,2,var) ) / TV,
  n_cl= length(decomp$hclustering_list),
  k=k, 
  noiseAdded=0,
  total_nPC= k*length(decomp$hclustering_list),
  sqrtCOR_mse= sqrt(MSE( C_o, cor( decomp$rec[,!noise_mask]) ))
	)
	gc(); return(res)
		} )-> summaryPerLVL_clean



saveRDS(c(summaryPerLVL_noise, summaryPerLVL_clean),
	sprintf('vary_nPC/decomp_MCL_nPC%d_summary.rds',k))


