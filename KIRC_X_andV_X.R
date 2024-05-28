X<-readRDS("KIRC_gene_expr_data.rds")
V_X=apply(X,2,var)
X=as.matrix(X)
hV= V_X >= quantile(V_X, 0.25)
X=X[,hV]
V_X=V_X[hV]
gc()

