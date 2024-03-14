load('afterDecomposition.RData')
source('blockwisePCA_R_engine.R')

set.seed(1234)
X_norm<-get_cloned_reconstruction(X=X,hcluster_blockPCA=rec_blockPCA, hclustering_list=rec_hclust, 
			  cloned_PCmatrix=PCmat_cloned_norm, noise_variances=NULL, uncenter=TRUE) 
X_norm_noise<-get_cloned_reconstruction(X=X,hcluster_blockPCA=rec_blockPCA, hclustering_list=rec_hclust, 
			  cloned_PCmatrix=PCmat_cloned_norm, noise_variances=noise_variances, uncenter=TRUE)

library(rmetalog)
PC_metalogs<- list()
PCmat_hidden_norm<- PCmat_cloned_norm
for (j in 1:ncol(PCmat_orig)) {
	PC_metalogs[[j]] <- readRDS( sprintf('PCs/metalog_%d.rds',j))
	PCmat_hidden_norm[,j] <- readRDS( sprintf('PCs/PChidden_normal%d.rds',j))
				}
target_terms= rep(5, ncol(PCmat_orig))
PCmat_cloned_meta<- cloned_metalogPC_matrix( PChidden_normal=PCmat_hidden_norm, 
					     target_metalogs=PC_metalogs, target_terms=target_terms)

X_meta<-get_cloned_reconstruction(X=X,hcluster_blockPCA=rec_blockPCA, hclustering_list=rec_hclust, 
			  cloned_PCmatrix=PCmat_cloned_meta, noise_variances=NULL, uncenter=TRUE) 
X_meta_noise<-get_cloned_reconstruction(X=X,hcluster_blockPCA=rec_blockPCA, hclustering_list=rec_hclust, 
			  cloned_PCmatrix=PCmat_cloned_meta, noise_variances=noise_variances, uncenter=TRUE) 

saveRDS(list(rec_blockPCA=rec_blockPCA,
rec_hclust=rec_hclust,
X_rec=X_rec,
X_rec_noise=X_rec_noise,
PCmat_orig=PCmat_orig,
PCmat_cloned_norm=PCmat_cloned_norm,
X_norm=X_norm,
X_norm_noise=X_norm_noise,
PCmat_hidden_norm= PCmat_hidden_norm,
PC_metalogs=PC_metalogs,
PCmat_cloned_meta=PCmat_cloned_meta,
X_meta=X_meta,
X_meta_noise=X_meta_noise), 'reconstruction_results.rds')
