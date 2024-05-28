
    tom<- function(A) {diag(A)=0; deg<-colSums(A);  (A + A%*%A)/( outer(deg,deg,pmin) +1 - A )  }

    dtom_hcl<- function(A,...) { TOM=tom(A); dt= 1- TOM; diag(dt)=0; rm(TOM); hclust(as.dist(dt), ...) }



    readRDS('KIRC_reconstruction_results.rds')-> reclist
source('KIRC_X_andV_X.R')
    noise_mask =  !( (1: ncol(X)) %in% c( reclist$rec_hclust[['1']],reclist$rec_hclust[['2']],reclist$rec_hclust[['3']]  ) )

readRDS( 'KIRC_plainSVDreconstruct93.rds')->SVD_res70
readRDS( 'KIRC_plainSVDreconstruct10.rds')->SVD_res10
# @ UP : list(PC=PC,
#	      X_rec=X_rec,
#	      X_rec_noise=X_rec_noise), 
colnames(SVD_res70$X_rec)<- colnames(X)
colnames(SVD_res70$X_rec_noise)<- colnames(X)


colnames(SVD_res10$X_rec)<- colnames(X)
colnames(SVD_res10$X_rec_noise)<- colnames(X)

    X_tables<- list(
                    rec=reclist$X_rec,
                    rec_noise=reclist$X_rec_noise,
                    cnorm=reclist$X_norm,
                    cnorm_noise=reclist$X_norm_noise,
                    cmeta=reclist$X_meta,
                    cmeta_noise=reclist$X_meta_noise)
    X_tables<-lapply(X_tables, function(x) x[,!noise_mask])
    X_tables$orig=X
    X_tables$SVD_93= SVD_res70$X_rec
    X_tables$SVD_10= SVD_res10$X_rec
    X_tables$SVD_93_noise= SVD_res70$X_rec_noise
    X_tables$SVD_10_noise= SVD_res10$X_rec_noise
    X_tables$WGCNA=readRDS( 'KIRC_WGCNA_sim_and_cl.rds')[[1]]
    X_tables_fancyNames<- c(
                    rec='hier. reconstruction',
                    rec_noise='hier. rec. + noise',
                    cnorm='prcedure I',
                    cnorm_noise='procedure I + noise',
                    cmeta='procedure II',
                    cmeta_noise='procedure II + noise',
		    orig='reference',
		    SVD_93='SVD: 93 PCs',
                    SVD_10='SVD: 10 PCs',
                    SVD_93_noise='SVD: 93 PCs + noise',
                    SVD_10_noise='SVD: 10 PCs + noise',
		    WGCNA='WGCNA'
			)


for (v in rev(names(X_tables))){
message(v)
COR<- cor(X_tables[[v]]);  gc()
dth_v<-  dtom_hcl(COR^2, method='average'); rm(COR);  gc()
saveRDS(dth_v, paste0('KIRC_',v,'_dtom_hclust.rds'))
rm(dth_v);
}
