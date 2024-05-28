readRDS('KIRC_reconstruction_results.rds')-> reclist
source('KIRC_X_andV_X.R')
WGCNA_res<- readRDS('KIRC_WGCNA_sim_and_cl.rds')
WGCNA_tab=WGCNA_res[[1]]
WGCNA_tab<- scale(WGCNA_tab, center=TRUE, scale=TRUE)
for (cname in colnames(WGCNA_tab))
    WGCNA_tab[, cname]=  (WGCNA_tab[, cname]*sd(X[,cname]) )+mean(X[,cname]) 

readRDS( 'KIRC_plainSVDreconstruct93.rds')->SVD_res70
# @ UP : list(PC=PC,
#	      X_rec=X_rec,
#	      X_rec_noise=X_rec_noise), 
colnames(SVD_res70$X_rec)<- colnames(X)
colnames(SVD_res70$X_rec_noise)<- colnames(X)


readRDS( 'KIRC_plainSVDreconstruct10.rds')->SVD_res10
# @ UP : list(PC=PC,
#	      X_rec=X_rec,
#	      X_rec_noise=X_rec_noise), 
colnames(SVD_res10$X_rec)<- colnames(X)
colnames(SVD_res10$X_rec_noise)<- colnames(X)



X_tables_mcl<- list(
                rec=reclist$X_rec,
                cnorm=reclist$X_norm,
                cmeta=reclist$X_meta,
                rec_noise=reclist$X_rec_noise,
                cnorm_noise=reclist$X_norm_noise,
                cmeta_noise=reclist$X_meta_noise)


noise_mask =  !( (1: ncol(X)) %in% c( reclist$rec_hclust[['1']],reclist$rec_hclust[['2']],reclist$rec_hclust[['3']]  ) )

X_tables_mcl<-lapply(X_tables_mcl, function(x) x[,!noise_mask])

X_tables<- c( list( orig=X,
		    WGCNA=WGCNA_tab,
		    SVD_93=SVD_res70$X_rec, 
		    SVD_10=SVD_res10$X_rec 
		  ),
	      X_tables_mcl[1:3],
	     list( SVD_93_noise=SVD_res70$X_rec_noise,
	           SVD_10_noise=SVD_res10$X_rec_noise
		),
	      X_tables_mcl[4:6]
	    )

gc()
message('loaded data')
source('blockwisePCA_R_engine.R')
C_orig<- cor(X_tables$orig) 
message('computed original correlations')
 
mse_deg_cc<- lapply(names(X_tables), function(variant)
#for (variant in names(X_tables))
	 { message(variant)
		if( variant=='orig') { meanSE=0; C=C_orig; diag(C)=0; clcoef= ccfast(C^2)
				       deg= colSums(C^2); rm(C); gc()
	         } else {
	   		 Xtab= X_tables[[variant]]
			 Xtab_subset<- colnames(X) %in% colnames(Xtab)
			 C=cor(Xtab)
			 meanSE= MSE( C, C_orig[Xtab_subset, Xtab_subset] )
			 diag(C)=0;
			 clcoef=ccfast(C^2)
			 deg= colSums(C^2); rm(C); gc()
			}	
		list( mse= meanSE,
		      cc= clcoef,
                      deg=deg)#-> partial_res
		#saveRDS(partial_res,paste0('cHeavyM',variant,'.rds'))
				}
		)

saveRDS(mse_deg_cc, 'KIRC_mse_deg_cc.rds')  


