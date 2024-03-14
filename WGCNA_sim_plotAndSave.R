load('WGCNA_simulation.RData')    

    plot_dth<- function(dth)   plot(dth, labels = FALSE, check=FALSE)


X.S<- X.S[, cl!=0]
saveRDS(list(X.S, cl), 'WGCNA_sim_and_cl.rds')
TM= tom(cor(X.S)^2)
message('done TOM')
gc()
dtom_hcl(TM,'average')-> dhcl
message('done dhcl')
gc()
rm(TM)

jpeg(filename= paste0('WGCNA','.jpg'), width=600, height=450)
plot_dth(dhcl)
dev.off()


