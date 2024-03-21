library(Hmisc)
library(Matrix) 
library(data.table)
source('blockwisePCA_R_engine.R')

args=commandArgs(trailingOnly=TRUE)
k = as.integer(args[[1]])

source('getXandV_X.R')

cldir=file.path(getwd(),'clusterings')
init_similarity_matrix_name='S_X'
 file.path(cldir, paste0(init_similarity_matrix_name,'.m0') )-> s_path
s_path_exists<- file.exists(s_path)
if ( !s_path_exists) {
print('no initial similarity, computing correlations')
rcorr_adjusted(X, rcorr.method='pearson', p.adjust.method="holm", 
			  diag.fill=0, r.diag.fill=NULL)-> rcorr.X
message('computed correlation')
cor_similarity( rcorr.adj=rcorr.X, addLoops=TRUE )-> S.X

 } else {S.X= read_mcl_abc2matrix(file_path=s_path, 
			          domain_size=ncol(X)); print('loaded initial similarity') }

initial_clustering_and_blockwisePCA(X=X,X_similarity_matrix=S.X, clustering_function=mcl_end2end,
					       min_cl_size=100, clustering_functionOtherArgs=list( inflation=2),
						f_E=0.7, X_variances=V_X,k=k)-> init_objects
init_objects$C_base -> C_base
init_objects$f -> f
init_objects$C_base_blockPCA -> C_base_blockPCA


set.seed(1234*k)

gs<- list( list(hclustering_list=C_base, hcluster_blockPCA=C_base_blockPCA) ) 
gs[[1]]$rec_noise<-get_full_reconstruction(X,gs[[1]]$hcluster_blockPCA, gs[[1]]$hclustering_list, 
add_noise=TRUE, uncenter=TRUE) 
gs[[1]]$rec<-get_full_reconstruction(X,gs[[1]]$hcluster_blockPCA, gs[[1]]$hclustering_list, 
add_noise=FALSE, uncenter=TRUE) 
for (a in 2:6) {
print(sprintf('building partition number %d',a))
subcluster_and_blockwisePCA3(X=X, X_similarity_matrix=S.X, 
				X_variances=V_X, hclustering_list=gs[[a-1]]$hclustering_list, 
				hcluster_PCA_blocks=gs[[a-1]]$hcluster_blockPCA,
				 k=k, f=0.7,
				clfun2=similarity_based_hclust,
				 clfun2OtherArgs_constant=list(method='complete'),
				  clfun2OtherArgs_ranges=list(n_group=c(7:2)))-> gs[[a]]
gs[[a]]$rec_noise<-get_full_reconstruction(X,gs[[a]]$hcluster_blockPCA, gs[[a]]$hclustering_list, 
add_noise=TRUE, uncenter=TRUE) 
gs[[a]]$rec<-get_full_reconstruction(X,gs[[a]]$hcluster_blockPCA, gs[[a]]$hclustering_list, 
add_noise=FALSE, uncenter=TRUE) 
}
saveRDS(gs, sprintf("./vary_nPC/decomposition_MCL_nPC%d.rds",k))
