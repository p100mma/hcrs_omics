### results of BRCA_simulation.R and BRCA_decomposition.R
readRDS('reconstruction_results.rds')-> reclist


source('getXandV_X.R')

### results of WGCNA_simulation.R
WGCNA_res<- readRDS('WGCNA_sim_and_cl.rds')
WGCNA_tab=WGCNA_res[[1]]
### scale WGCNA simulation to adjust means and variances
### only useful for KS distance criteria (supplementary marerial)
WGCNA_tab<- scale(WGCNA_tab, center=TRUE, scale=TRUE)
for (cname in colnames(WGCNA_tab))
    WGCNA_tab[, cname]=  (WGCNA_tab[, cname]*sd(X[,cname]) )+mean(X[,cname])



### results of plain_SVD.R
readRDS( 'plainSVDreconstruct70.rds')->SVD_res70
# @ UP : list(PC=PC,
#             X_rec=X_rec,
#             X_rec_noise=X_rec_noise),
colnames(SVD_res70$X_rec)<- colnames(X)
colnames(SVD_res70$X_rec_noise)<- colnames(X)


readRDS( 'plainSVDreconstruct10.rds')->SVD_res10
# @ UP : list(PC=PC,
#             X_rec=X_rec,
#             X_rec_noise=X_rec_noise),
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
                    SVD_70=SVD_res70$X_rec,
                    SVD_10=SVD_res10$X_rec
                  ),
              X_tables_mcl[1:3],
             list( SVD_70_noise=SVD_res70$X_rec_noise,
                   SVD_10_noise=SVD_res10$X_rec_noise
                ),
              X_tables_mcl[4:6]
            )

gc()
message('loaded data')
source('blockwisePCA_R_engine.R')

## result of computation_heavy_metrics.R
mse_cc_deg<- readRDS('mse_deg_cc.rds')
names(mse_cc_deg)<-names(X_tables)

mse_cc_deg[[1]]-> reference_stats

######PLOT setup
########labels
plot_selection<- c(2:length(X_tables))
test_stats<- mse_cc_deg[ plot_selection]
names(test_stats)<- c('WGCNA', 'SVD70','SVD10','HR-5-3-11','PI-5-3-11','PII-5-3-11',
                       'SVD70+N','SVD10+N', 'HR-5-3-11+N','PI-5-3-11+N','PII-5-3-11+N')
test_titles<- c('WGCNA', 'SVD70','SVD10','HCR-5-3-11','HCS(n)-5-3-11','HCS(f)-5-3-11',
                       'SVD70+N','SVD10+N', 'HCR-5-3-11+N','HCS(n)-5-3-11+N','HCS(f)-5-3-11+N')

########colors
Tol_muted <- c('#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499', '#DDDDDD')
ref_color=Tol_muted[[2]]
t_color_RGB<-col2rgb(Tol_muted[[6]])/255
tcArgs<- as.list(t_color_RGB)
names(tcArgs)<- rownames(t_color_RGB)
tcArgs$alpha=0.6
test_color= do.call(rgb,tcArgs)

#########clustering coefficient (via WGCNA approach)
w=360
CEX=1.6
for (i in 1:length(test_stats)){
    jpeg(paste0('cc_',names(test_stats)[[i]],'.jpg'), width= w, height= 360 )
    par(mar = c(4, 5, 2, 1) + 0.1)
    gmin<-min(unlist(lapply(list(reference_stats, test_stats[[i]] ), function(x) min(x$cc) )))
    gmax<-max(unlist(lapply(list(reference_stats, test_stats[[i]] ), function(x) max(x$cc) )))
    gbreaks= seq(from=gmin, to=gmax, length.out = 20)
    table(cut(test_stats[[i]]$cc, gbreaks))-> heights_test
        heights_test= heights_test/sum(heights_test)
    table(cut(reference_stats$cc, gbreaks))-> heights_ref
        heights_ref= heights_ref/sum(heights_ref)
    YL= max(heights_ref, heights_test)+0.2
   barplot(height= heights_ref, xaxt='n' ,ylim=c(0,YL), main=test_titles[[i]], xlab='clustering coefficient',
          ylab='P(cc)',col=ref_color ,
           cex.lab=CEX*1.2, cex.axis=CEX, cex.main=CEX*1.3, cex.sub=CEX)
    barplot(height= heights_test, xaxt='n' , yaxt='n', ylim=c(0,YL), col=test_color, add=TRUE )
     axis(side = 1, at = seq_along(gbreaks) -1, tick = FALSE, labels = round(gbreaks,2), cex.axis=CEX)
        legend('topright', legend=c('reference', test_titles[[i]]), 
                                   y.intersp = 1, pch=16, 
                                   col= c(ref_color,test_color), cex=CEX*1.3,
                                      bg='transparent',bty='n')
    dev.off()
    }

#########weighted degree (aka Connectivity via WGCNA approach)
w=360
CEX=1.6
for (i in 1:length(test_stats)){
    jpeg(paste0('deg_',names(test_stats)[[i]],'.jpg'), width= w, height= 360 )
    par(mar = c(4, 5, 2, 1) + 0.1)
    gmin<-min(unlist(lapply(list(reference_stats, test_stats[[i]] ), function(x) min(x$deg) )))
    gmax<-max(unlist(lapply(list(reference_stats, test_stats[[i]] ), function(x) max(x$deg) )))
    gbreaks= seq(from=gmin, to=gmax, length.out = 20)
    table(cut(test_stats[[i]]$deg, gbreaks))-> heights_test
        heights_test= heights_test/sum(heights_test)
    table(cut(reference_stats$deg, gbreaks))-> heights_ref
        heights_ref= heights_ref/sum(heights_ref)
    YL= max(heights_ref, heights_test)+0.2
   barplot(height= heights_ref, xaxt='n' ,ylim=c(0,YL), main=test_titles[[i]], xlab='weighted degree',
          ylab='P(deg)',col=ref_color ,
           cex.lab=CEX*1.2, cex.axis=CEX, cex.main=CEX*1.3, cex.sub=CEX)
    barplot(height= heights_test, xaxt='n' , yaxt='n', ylim=c(0,YL), col=test_color, add=TRUE )
     axis(side = 1, at = seq_along(gbreaks) -1, tick = FALSE, labels = round(gbreaks,2), cex.axis=CEX)
        legend('topright', legend=c('reference', test_titles[[i]]), 
                                   y.intersp = 1, pch=16, 
                                   col= c(ref_color,test_color), cex=CEX*1.3,
                                      bg='transparent',bty='n')
    dev.off()
    }


##########################################################################################
################################ KS distances#############################################
##########################################################################################

lapply(names(X_tables[2:length(X_tables)]), function(tn)
      { Xtab<- X_tables[[tn]]
        message(tn)  
     unlist(lapply(colnames(Xtab), function(cn)
        {
             ref_ecdf<- ecdf(X_tables$orig[, cn  ])
         ecdf_distance(y= Xtab[,cn], ecdf_x= ref_ecdf, summarizeFun = max)}
        ))     
      }
    
      )->X_tables_KSdistances

######## HCS(f)
KS_distances_PC<- lapply(1:ncol(reclist$PCmat_orig), function(j){
    ecdf(reclist$PCmat_orig[,j])-> ref_ecdf
    ecdf_distance(y= reclist$PCmat_cloned_meta[,j], ecdf_x = ref_ecdf, summarizeFun = max)
    
    }
                         )

######## HCS(n)
KS_distances_PCI<- lapply(1:ncol(reclist$PCmat_orig), function(j){
    ecdf(reclist$PCmat_orig[,j])-> ref_ecdf
    ecdf_distance(y= reclist$PCmat_cloned_norm[,j], ecdf_x = ref_ecdf, summarizeFun = max)
    
    }
                         )


rbind(`PCs: procedure I`=summary(unlist(KS_distances_PCI) ),
    `PCs: procedure II`=summary(unlist(KS_distances_PC) ),
      do.call(rbind,lapply(X_tables_KSdistances,summary)))-> KS_distances_tab

write.table(x=round(KS_distances_tab,2), 'KS_distances.csv',quote = FALSE,sep = ',')
