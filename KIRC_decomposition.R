m_C=5
m_Csplit=6
library(Hmisc)
library(Matrix) 
library(data.table)
source('blockwisePCA_R_engine.R')
source('KIRC_X_andV_X.R')

cldir=file.path(getwd(),'clusterings')
init_similarity_matrix_name='KIRC_S_X'
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
					       min_cl_size=100, clustering_functionOtherArgs=list( inflation=2,
								similarity_matrix_name=init_similarity_matrix_name,
								output_name='KIRC_base'),
						f_E=0.7, X_variances=V_X,k=5)-> init_objects
init_objects$C_base -> C_base
init_objects$f -> f
init_objects$C_base_blockPCA -> C_base_blockPCA



gs<- list( list(hclustering_list=C_base, hcluster_blockPCA=C_base_blockPCA) ) 
for (a in 2:6) {
print(sprintf('building partition number %d',a))
subcluster_and_blockwisePCA3(X=X, X_similarity_matrix=S.X, 
				X_variances=V_X, hclustering_list=gs[[a-1]]$hclustering_list, 
				hcluster_PCA_blocks=gs[[a-1]]$hcluster_blockPCA,
				 k=5, f=0.7,
				clfun2=similarity_based_hclust,
				 clfun2OtherArgs_constant=list(method='complete'),
				  clfun2OtherArgs_ranges=list(n_group=c(7:2)))-> gs[[a]]
if (!is.null(gs[[a]]$hclustering_list)) gs[[a]]$rec<-get_reconstruction_at_g(X,gs[[a]]$hcluster_blockPCA, gs[[a]]$hclustering_list, g=a) 
}
set.seed(12412)
gs[[2]]$hcluster_blockPCA-> rec_blockPCA
gs[[2]]$hclustering_list -> rec_hclust
X_rec<- get_full_reconstruction(X=X,hcluster_blockPCA=rec_blockPCA,
	 hclustering_list=rec_hclust, add_noise=FALSE, uncenter=TRUE) 

noise_variances<- apply( X-X_rec,2,var) 
X_rec_noise<- get_full_reconstruction(X=X,hcluster_blockPCA=rec_blockPCA,
	 hclustering_list=rec_hclust, add_noise=TRUE, uncenter=TRUE) 

PCmat_orig<- PC_generator_matrix(hcluster_blockPCA= rec_blockPCA)
PCmat_cloned_norm<-cloned_normalPC_matrix(reference_PC_gen_matrix=PCmat_orig)
saveRDS(PCmat_cloned_norm, "PCs/KIRC_PCmat_cloned_norm.rds")
saveRDS(PCmat_orig, "PCs/KIRC_PCmat_orig.rds")
message('saving R session image')
save.image('KIRC_afterDecomposition.RData')

