load('KIRC_WGCNA_clustering.RData')
library(dynamicTreeCut)
TM<- 1 - TM; diag(TM)=0 #to distance conversion
cl<- cutreeDynamic(dendro=dhcl, method="hybrid",distM=TM, minClusterSize=150, deepSplit=4)

simulate_cluster<- function( domain, cluster_location){
colMeans(domain[,cluster_location])->cl_means
with(svd( scale(domain[,cluster_location],center=TRUE,scale=FALSE) ), u[,1] * rep(d[[1]], nrow(u) ) )-> pc
pc_cor<- cor(pc, domain[, cluster_location])
n_neg<- sum(pc_cor < 0)
v_p<- var(pc)
simulated<- domain[,cluster_location]
for (j in 1:length(pc_cor)) # cor(sim_X,pc)^2== cor(X,pc)^2
	{
	numer= 1 - pc_cor[[j]]^2
	denom= pc_cor[[j]]^2
	 var_e<- v_p *(numer/denom)
	simulated[,j]= pc + rnorm( length(pc), 0, 
				   sd=  sqrt( 
					   var_e 
					    ) 
				  )
	}
neg_idx<- sample(1:ncol(simulated), n_neg,replace=FALSE)
for (j in neg_idx)
	simulated[,j]= simulated[,j]*(-1)
return(simulated)
}
set.seed(1234)
X.S<-X
X.S[,]=0
non0cl<- unique( cl[cl!=0] )
for (lbl in non0cl)
	X.S[,which(cl==lbl)]<-simulate_cluster( X, which(cl==lbl) )
rm(TM)
save.image('KIRC_WGCNA_simulation.RData')
X.S<- X.S[, cl!=0]
saveRDS(list(X.S, cl), 'KIRC_WGCNA_sim_and_cl.rds')
