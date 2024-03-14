
summaries<- list()

for (k in 2:5) summaries[[k-1]]= readRDS(sprintf('vary_nPC/decomp_MCL_nPC%d_summary.rds',k))
summaries[[length(summaries)+1]] = readRDS('./sim_MCL_nPC5_summary.rds')
summaries<- lapply(summaries, function(x) do.call(rbind,x))
do.call(rbind, summaries)-> summry
summry<- as.data.frame(summry)
summry[, 'g' ] = c( rep(1:6, 8), rep(2,4) )
#print(summry)
summry= summry[, c('k', 'g', 'n_cl', 'total_nPC', 'noiseAdded', 'var_explained',
		   'btwMaxCS', 'btwMeanCS','sqrtCOR_mse') ]
summry= summry[ summry$g %in% 1:4, ]
summry$var_explained= round(summry$var_explained,2)
summry$sqrtCOR_mse= round(summry$sqrtCOR_mse,2)
summry$btwMaxCS= round(summry$btwMaxCS,4)
summry$btwMeanCS= round(summry$btwMeanCS,4)
summry$btwMaxCS[ !is.na(summry$n_cl)]= round(summry$btwMaxCS[!is.na(summry$n_cl)],2)
summry$btwMeanCS[ !is.na(summry$n_cl)]= round(summry$btwMeanCS[!is.na(summry$n_cl)],2)
write.csv(summry, file='vary_nPC/comparison_summary.csv', col.names=FALSE, row.names=FALSE, sep=',', quote=FALSE)


