    plot_dth<- function(dth,...)   plot(dth, labels = FALSE, check=FALSE,xlab="",ylab="",sub="",...)
   

X_tables_fNames<- c(
                    rec='HR-5-4-15',
                    rec_noise='HR-5-4-15+N',
                    cnorm='PI-5-4-15',
                    cnorm_noise='PI-5-4-15+N',
                    cmeta='PII-5-4-15',
                    cmeta_noise='PII-5-4-15+N',
		    orig='REF',
		    SVD_93='SVD93',
                    SVD_10='SVD10',
                    SVD_93_noise='SVD93+N',
                    SVD_10_noise='SVD10+N',
		    WGCNA='WGCNA'
			)


 
X_tables_fancyNames<- c(
                    rec='HCR-5-4-15',
                    rec_noise='HCR-5-4-15+N',
                    cnorm='HCS(n)-5-4-15',
                    cnorm_noise='HCS(n)-5-4-15+N',
                    cmeta='HCS(f)-5-4-15',
                    cmeta_noise='HCS(f)-5-4-15+N',
		    orig='REF',
		    SVD_93='SVD93',
                    SVD_10='SVD10',
                    SVD_93_noise='SVD93+N',
                    SVD_10_noise='SVD10+N',
		    WGCNA='WGCNA'
			)


for (v in seq_along(X_tables_fancyNames)){
plain_name<- names(X_tables_fancyNames)[[v]]
dth_v<-readRDS( paste0('KIRC_',plain_name,'_dtom_hclust.rds'))
jpeg_fname= X_tables_fNames[[v]]
jpeg(filename= paste0('KIRC_hclust_',jpeg_fname,'.jpg'), width=500, height=450)
par(cex=1.7)
par(mar = c(1, 2, 2, 1) + 0.1)
plot_dth(dth_v, main=NA)
title(X_tables_fancyNames[[v]], line=-0.5)
dev.off()
message(jpeg_fname)
}
