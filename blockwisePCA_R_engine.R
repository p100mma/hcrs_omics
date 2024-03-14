
###############################################################################################################
######CORR BASED SIMILARITY CALCULATION########################################################################
###############################################################################################################

#requires Hmisc
rcorr_adjusted<- function(X, rcorr.method='pearson', p.adjust.method="holm", 
			  diag.fill=NA, r.diag.fill=NULL){
cor.data<- rcorr(X)
cor.data$Pa<- cor.data$P
ltr.select<- lower.tri(cor.data$Pa)
padj<- p.adjust( cor.data$Pa [ ltr.select ] )
cor.data$Pa[ ltr.select ] <-  padj
cor.data$Pa<- t(cor.data$Pa)
cor.data$Pa[ ltr.select ] <- padj
cor.data$Pa<- t(cor.data$Pa)
if (!is.na(diag.fill))
for (k in 1:length(cor.data)){

if (!is.null( dim(cor.data[[k]]) )  )
	diag(cor.data[[k]]) = diag.fill
}
if (!is.null(r.diag.fill)) diag(cor.data$r)<- r.diag.fill
cor.data
}

cor_similarity<- function( rcorr.adj, addLoops=TRUE ) {
S<- rcorr.adj$r ^ 2
S[ rcorr.adj$Pa >= 0.05 ] = 0
if (addLoops) diag(S)<-1
S
}

###############################################################################################################
######MCL INTERACTION FUNCTIONS################################################################################
###############################################################################################################

#requires Matrix
mcl_abc_format<- function(S) {
S.s<- as(S, "sparseMatrix")
S.sl<-as.matrix( as.data.frame(summary(S.s)))
for (j in 1:2) S.sl[,j] = S.sl[,j] - 1 
not_loops<- S.sl[ S.sl[,1]!=S.sl[,2],  ]
nS.sl<- matrix(ncol=3, nrow= nrow(S.sl) + nrow(not_loops) )
nS.sl[1:nrow(S.sl),] = S.sl
nS.sl[c(nrow(S.sl)+1):nrow(nS.sl) ,] =  not_loops[, c(2,1,3) ]
as.data.frame(nS.sl)
}

prepare_MCL_qsub<- function(taskPath, taskName, headerTemplatePath, workdir,
			    inputPath, outputPath, inflation){
con= file(headerTemplatePath)
readLines(con)-> taskLines
close(con)
taskLines[[ length(taskLines) + 1 ]]<-  sprintf('# PBS -N %s',taskName) 
taskLines[[ length(taskLines) + 1 ]]<-  sprintf('cd %s',workdir)
taskLines[[ length(taskLines) + 1 ]]<- sprintf('mcl %s --abc  -I %.2f -scheme 7 -o %s', inputPath, inflation, outputPath)
con<- file(taskPath)
writeLines(taskLines, con)
close(con)
}

run_qsub<- function(taskPath) {

system( sprintf(' qsub %s', taskPath) )

}


#requires data.table
write_mcl_abc2file<- function( edge_df, file_path) {
fwrite(edge_df,file_path , sep=" ", row.names=FALSE, col.names=FALSE)
}

#requires data.table & matrix
mcl_abc2matrix<- function(edge_df, domain_size){
as.matrix(sparseMatrix( i=edge_df[[1]], j=edge_df[[2]], x= edge_df[[3]], index1=FALSE), dims=c(domain_size, domain_size) )
}

read_mcl_abc2matrix<- function(file_path, domain_size) {
edge_df<- fread( file_path, sep=" ")
mcl_abc2matrix(edge_df, domain_size)
}

run_mcl<- function(inputPath, outputPath, inflation){
system(sprintf('~/local/bin/mcl %s --abc  -I %.2f -scheme 7 -o %s', inputPath, inflation, outputPath))
}

read_mcl_clustering<- function( file_path, domain_size ){

clustering<- rep(NA,domain_size)

lines<-readLines(con<- file(file_path))
close(con)
cl_list<- lapply( strsplit( lines, "\t"), as.integer)

for (j in 1:length(cl_list))
{	
if( length( cl_list[[j]])>1)
clustering[ cl_list[[j]] + 1 ] = j else clustering [ cl_list[[j]] + 1 ]=0 
}
clustering[ is.na(clustering) ] = 0
clustering
}


mcl_end2end<- function( similarity_matrix, 
			inflation,
			files_dir=file.path(getwd(),'clusterings'),similarity_matrix_name='S_X',
			output_name='base')

{
sim_mat_fname<- paste0(similarity_matrix_name,'.m0')
in.file_path= file.path( files_dir, sim_mat_fname) 
o_fname<- paste0('C',output_name ,'p', inflation,'.mcl')
o.file_path= file.path( files_dir, o_fname  ) 
if (!file.exists(in.file_path)){
	mcl_abc_format(similarity_matrix)-> similarity_edges
	write_mcl_abc2file( edge_df=similarity_edges, 
		     file_path= in.file_path)
			    }
if (!file.exists(o.file_path) ) run_mcl(inputPath=in.file_path, 
				    outputPath= o.file_path, inflation=inflation)
clV<-read_mcl_clustering(o.file_path, domain_size= ncol(similarity_matrix))
attr(clV, "params")= list(inflation=inflation)
clV
}


###############################################################################################################
######BLOCKWISE PCA RECONSTRUCTION#############################################################################
###############################################################################################################

fastTable<- function(values)
{
    uqs<-unique(values)
    histo<-unlist(lapply(uqs, function(x) sum(values==x)))
    list(value=uqs,
         count=histo)
}



blockwise_PCA_reconstruction<- function(A,A_variances, clustering_list, f, k){

   lapply(clustering_list, function(C_i)
        {
            A.cl<- A[, C_i, drop=FALSE]
            A.cl<-scale(A.cl, center=TRUE, scale=FALSE)
              # K<- floor(ncol(A.cl)*perc)
	    sum(A_variances[ C_i])-> var.total
	    svd(A.cl)-> svd_A.cl
	    PC_A.cl<- with(svd_A.cl,  u %*% diag(d, ncol=length(d), nrow=length(d)) )
	    pc.vars<-apply(PC_A.cl,2,var)
            cumsum(pc.vars)-> pc.vars.cumul
	    if( length(C_i)==1 ) {k_i =1; k2use=1} else {
            k_i <-min(which(pc.vars.cumul/var.total >= f  ))
	    k2use<- min(k_i, k) }
	    lc_coefs_ALL<- t(svd_A.cl$v)
	    lc_coefs<- lc_coefs_ALL[1:k2use,, drop=FALSE]
	    PC_rec<- PC_A.cl[, 1: k2use, drop=FALSE ]
	    PC_rec %*%  lc_coefs -> A.cl_rc
#	    print( 'var reconstr - var PC[1:k_used]:')
#	    print( sum(apply(A.cl_rc,2,var)) - sum(pc.vars.cumul[k2use]))
            list(location=C_i,
                 total_var=var.total,
		 pc_var=pc.vars,
		 cumul_var=pc.vars.cumul,
                 k_used=k2use,
		 k_f=k_i,
		 f=f,
		 saturated= k_i <= k2use,
                 var_explained= pc.vars.cumul[[k2use]],
                 PC= PC_rec,
                 Vt=lc_coefs,
		 all_Vt=lc_coefs_ALL,
		 reconstruction=A.cl_rc
                 )
	})-> blocks_list; names(blocks_list)<- names(clustering_list)
		return(blocks_list)
								}


allocate_blocks<-function( where, blockwise_PCA){

    for (block in blockwise_PCA){
        where[, block$location] <- block$reconstruction
			    }
where
						}





###############################################################################################################
######INITIAL STEPS AND INITIAL SUBDIVISION####################################################################
###############################################################################################################


noise_separation<- function(clustering_vector, min_cl_size){
 cluster_counts<- fastTable(clustering_vector)
 for (i in seq_along(cluster_counts[[1]]))
	if (cluster_counts$count[[i]]< min_cl_size) 
		clustering_vector[ clustering_vector== cluster_counts$value[[i]] ] = 0
 return(clustering_vector)
}

# reindex nonzero labels to 1:N
tidyUpLabels<- function(labelV) {
uqn0<-unique(labelV[labelV!=0])
rmentL<- seq_along(uqn0)
rment<- rep(0, length(labelV))
for (l in seq_along(uqn0))
   rment[ labelV == uqn0[[l]] ]= rmentL[[l]]
attr(rment, "params")= attr(labelV, "params")
return(rment)
}

clustering_vector2list<- function(clustering_vector){

labels<- sort( unique(clustering_vector[clustering_vector!=0]))
clustering_list<-list()
for (i in labels)
	clustering_list[[i]]<- which(clustering_vector==i)
attr(clustering_list, 'domain_size')= length(clustering_vector)
attr(clustering_list, 'params')= attr(clustering_vector,'params')
names(clustering_list)<- as.character(seq_along(clustering_list))
clustering_list
}

clustering_list2vector<- function(clustering_list){
clustering_vector<- rep(0, attr(clustering_list, 'domain_size'))
for (i in seq_along(clustering_list))
	clustering_vector[ clustering_list[[i]] ] = i
clustering_vector
}


calculate_f<- function( clustering_vector, f_E, X_variances){

 v_total<- sum(X_variances)
 v_noiseless<- sum(X_variances[ clustering_vector!= 0 ])
 v_target<- f_E * v_total
 if (v_noiseless < v_target) {
	warning("initial target fraction f_E cannot be reached by using any fraction of noiseless clustering, returning NA")
	return(NA)
    }
 stopifnot(v_noiseless>0) 
 return( v_target/ v_noiseless)	
}

cluster_and_separate_noise<- function( similarity_matrix, clustering_function, min_cl_size, clustering_functionOtherArgs) {

clfARGS<- c(list( similarity_matrix), clustering_functionOtherArgs)
cl_vec<- do.call(clustering_function, clfARGS)
noise_separation(clustering_vector=cl_vec, min_cl_size=min_cl_size)-> cl_vec
tidyUpLabels(cl_vec)-> cl_vec
clustering_vector2list(cl_vec)
}


initial_clustering_and_blockwisePCA<- function(X,X_similarity_matrix, clustering_function,
					       min_cl_size, clustering_functionOtherArgs,
						f_E, X_variances,k){
C_base<-cluster_and_separate_noise( X_similarity_matrix, clustering_function, min_cl_size, clustering_functionOtherArgs) 
f=calculate_f( clustering_list2vector(C_base), f_E, X_variances)
if (is.na(f)) {warning(" stopped because provided clustering leaves too much variance in the noise part");
	       return( list(C_base=C_base,
			    f=f,
			    C_base_blockPCA=NA) ) }
print( lapply(C_base, length) )
print( f )
list(C_base=C_base, f=f,
C_base_blockPCA=blockwise_PCA_reconstruction(X,X_variances, clustering_list=C_base, f=f, k=k)
	)
} 

.p<- function(...) paste0(c(...),collapse=',')

subdivide_cluster<- function( hclustering_list,
			      cluster_tuple,
                                subdivision){

if (length(hclustering_list)==0) stop(" hclustering_list must be nonempty") 
cluster_tuple<- .p( cluster_tuple)
if (!(cluster_tuple %in% names(hclustering_list))) 
	stop("cluster to divide is not an element of hclustering_list!")
 stopifnot( all ( unlist(subdivision) %in% hclustering_list[[ cluster_tuple ]] ) )
 stopifnot( all ( unlist(lapply(subdivision, length)) > 0 ) )
 stopifnot(  length( unlist( unique(subdivision)) ) == length(unlist(subdivision)))
 for (i in seq_along(subdivision))
	hclustering_list [[ .p(cluster_tuple, i) ]]<- subdivision[[i]]
attr(subdivision,'params')-> s_params
names(s_params)<- paste0( cluster_tuple, names(s_params))
attr(hclustering_list, 'params') <- c(attr(hclustering_list, 'params'), s_params)
 hclustering_list
}
#see the last argument of the function initial_subcluster for details
mcl_OtherArgsDependence<- function(...) {
values2use<- list(...)
values2use[[1]] -> g # the number of the reconstruction right now examined (g=1,2 ... G_TOTAL)
values2use[[2]] -> cl_tuple # indentifier of the cluster to subdivide
list(
similarity_matrix_name= sprintf('S_%s.E%s', cl_tuple, g),
output_name= cl_tuple)-> dependentArgs
if (length(values2use) > 2) {
 finalArgs<- c( values2use[ 3: length(values2use) ] , dependentArgs ) 
} else { finalArgs<- dependentArgs }
return(finalArgs)
}

# if >0 then clusters are too small
assess_size_distribution<- function( clustering_list,varLeft2explain, residuals) {
	if (!length(clustering_list))   #no noiseless part
		return(2)


	#check the variance of noiseless residuals (those in C_i_subdivision) vs the variance left to explain
	noiseless_var<- sum(apply(  residuals[,unlist(clustering_list)], 2, var ))

	if ( noiseless_var < varLeft2explain)   #too small noiseless part
		return(1)

	if (length(clustering_list)==1)  #only 1 big cluster in noiseless part but with enough variance
		return(0)

	return(-1)
}

similarity_based_hclust<- function( similarity_matrix, n_group, max_sim=1, method="complete") {
hclust( as.dist( max_sim - similarity_matrix), method=method)-> hcl_object
clvec<-cutree(hcl_object, k=n_group)
attr(clvec, 'params')=list(n_group=n_group)
attr(clvec, 'domain_size')= ncol(similarity_matrix)
clvec
} 


normalized_cl_entropy<- function( clustering_list ) {
 unlist(lapply(clustering_list, length))-> n_c
 p_c<- n_c/sum(n_c)
 length(n_c)-> n_groups
 -sum( p_c * log( p_c ) )/(log(n_groups))
}

initial_subcluster_and_blockwisePCA<- function(X, X_similarity_matrix, 
				X_variances, C_base, C_base_PCA_blocks, k, f,
				clfun1, clfun2, min_cl_size, mc_beforeStop=1.5*min_cl_size,
				 clfun1OtherArgs_functionalDependence=NULL,
				 clfun1OtherArgs_constant=NULL,
				 clfun1OtherArgs_ranges=NULL, 
				 clfun2OtherArgs_constant=NULL,
				  clfun2OtherArgs_ranges=NULL){


X1<- X
X1[,]=0  #placeholder for the reconstruction of 1st order
allocate_blocks( where=X1, blockwise_PCA= C_base_PCA_blocks)-> X1
X1_variances<- apply(X1, 2, var)
E1 = X - X1 #residuals
E1_variances<- apply(E1, 2, var)
for (bl in 1:length(C_base_blockPCA)) #if 'exhausted' status was not set for a block, declare it to be FALSE
	if (is.null(C_base_blockPCA[[bl]]$exhausted)) C_base_blockPCA[[bl]]$exhausted = FALSE
#theese lists will be expanded by additional blocks/subdivisions
C_expanded<- C_base
C_expanded_blockPCA<- C_base_blockPCA


unsaturated_indexes<- which( unlist(lapply(C_base_blockPCA, function(x) !x$saturated) ))

exhausted_clusters<- which( unlist(lapply(C_base_blockPCA, function(x) x$exhausted) ))
unsaturated_indexes<- setdiff ( unsaturated_indexes, exhausted_clusters)    

print(sprintf(" number of unsaturated indexes = %d ", length(unsaturated_indexes) ) )
for (i in unsaturated_indexes) {
	C_i <- C_base[[i]]
        S_E1.C_i <- cor ( E1 [ ,C_i ] )^2 #building similarity of residuals E1 in cluster C_i
	S_E1.C_i[ X_similarity_matrix[C_i, C_i] == 0 ]= 0 #carry zeros from the similarity of X
	
        left2explain= f* sum(X_variances[ C_i ] ) - sum(X1_variances[ C_i] )


	###################################################	
	##Clustering choice################################
	###################################################
	if ( length(C_i) < mc_beforeStop) {  # case of exhausted clusters, only 1 one split remains for them,
	# blockPCA entries are later marked as exhausted=TRUE
		clf2Args<- c(list( S_E1.C_i), clfun2OtherArgs_constant, list(n_group=2) )
		clustering_vector2list(do.call(clfun2, clf2Args))-> C_i_subdivision
		cis_attrs<- attributes(C_i_subdivision)
		lapply(C_i_subdivision, function(C_i.j)  C_i[ C_i.j ] ) -> C_i_subdivision #translate local to global indexes
		attributes(C_i_subdivision)<- cis_attrs

	} else { sizeDistrState= 9999; n_tries= length(clfun1OtherArgs_ranges[[1]]); i_try=1
		 clfun1args= clfun1OtherArgs_functionalDependence(g=1, cl_tuple=i )
		 clfun1args<- c(clfun1args, clfun1OtherArgs_constant) 
		while ((i_try < n_tries + 1) && (sizeDistrState>0)) {
			for (pname in names(clfun1OtherArgs_ranges))
				{
				clfun1args[[pname]] = clfun1OtherArgs_ranges[[pname]][[i_try]]
				}
			print(sprintf("trying cluster nr %d with settings:", i))
			print(clfun1args[ names(clfun1args) %in% names(clfun1OtherArgs_ranges) ] )	
			cluster_and_separate_noise( similarity_matrix= S_E1.C_i, 
					   clustering_function=clfun1, 
					   min_cl_size=min_cl_size, 
					   clustering_functionOtherArgs= clfun1args) -> C_i_subdivision; cis_attrs<- attributes(C_i_subdivision)
			
			lapply(C_i_subdivision, function(C_i.j)  C_i[ C_i.j ] ) -> C_i_subdivision #translate local to global indexes
			attributes(C_i_subdivision)<- cis_attrs
			assess_size_distribution( clustering_list=C_i_subdivision,
						varLeft2explain=left2explain, residuals=E1)-> sizeDistrState 
			i_try=i_try+1
			print(sprintf("sizeDistrState= %d", sizeDistrState))
			print(sprintf("itry = %d", i_try))
			print(sprintf("n_tries = %d", n_tries))
			print((i_try < n_tries + 1) && (sizeDistrState>0) )
		}

		if (sizeDistrState==0) {   # if TRUE, then  C_i_subdivision has one element and 
					   # and also, this one cluster has enough variance to meet the left2explain target value
					   # therefore we have to split it using clfun2 instead of clfun1
					   # we switch the domain to subdivide from C_i to C_{i=/=0}
					   C_i<- C_i_subdivision[[1]]  
        				   S_E1.C_i <- cor ( E1 [ ,C_i ] )^2  #have to recalculate similarity matrix
					    S_E1.C_i[ X_similarity_matrix[C_i, C_i] == 0 ]= 0 #carry zeros from the similarity of X
					   # check: maybe C_{i=/=0} has less that mc_beforeStop elements...
					   if ( length(C_i) < mc_beforeStop) { # if yes, use clfun2 as a last splitting step and mark the subdivision of C_i as exhausted later 
						    clf2Args<- c(list( S_E1.C_i), clfun2OtherArgs_constant, list(n_group=2) )
						    clustering_vector2list(do.call(clfun2, clf2Args))-> C_i_subdivision
						    cis_attrs<- attributes(C_i_subdivision)
						    lapply(C_i_subdivision, function(C_i.j)  C_i[ C_i.j ] ) -> C_i_subdivision #translate local to global indexes
						    attributes(C_i_subdivision)<- cis_attrs
					    } else {  # otherwise C_{i=/=0} has enough elements, so we can apply clfun2 with normal parameters, so set sizeDistrState <0     
						     sizeDistrState= 1
					    }
				       }	
		if (sizeDistrState>0) {    #in that case, clfun1 did not work (too small clusters/ too big noise part)
					    #have to use clfun2 ("percollator") instead (without noise / non-noise distinction)
					    #either on the whole initial C_i or substituted for C_{i=/=0} in previous sizeDistrState==0 clause.
			clf2Args<- c(list( S_E1.C_i), clfun2OtherArgs_constant )
			settings_scores<-vector()
			candidate_cl_lists<-list()
			settings_m<- length(clfun2OtherArgs_ranges[[1]])
			for (i_settings in 1:settings_m) {
				for (pname in names(clfun2OtherArgs_ranges))
					clf2Args[[pname]] = clfun2OtherArgs_ranges[[pname]][[i_settings]]	
				clustering_vector2list(do.call(clfun2, clf2Args))-> candidate_cl
				can_attrs<- attributes(candidate_cl)
				lapply(candidate_cl, function(C_i.j)  candidate_cl[ C_i.j ] ) -> candidate_cl #translate local to global indexes
				attributes(candidate_cl)<- can_attrs
				candidate_cl_lists[[ i_settings ]] = candidate_cl
				settings_scores [[ i_settings ]] =  normalized_cl_entropy( candidate_cl )     # "best" precollating clustering is the one which proportions are closest
													      # ... to groups of even size
				}
			C_i_subdivision<- candidate_cl_lists [[ which.max( settings_scores) ]]	
				       }
	} # end else from if (size(C_i) < minimalsizeBeforeStop)

	#######################################################	
	##Clustering choice END################################
	#######################################################
	#at this point, whatever subdivision happened, the result is sensible

	noiseless_C_i_var<- sum(apply(  E1[,unlist(C_i_subdivision)], 2, var ))
	# rename each cluster 'j' in subdivision from 'j' to 'i,j' 
	names(C_i_subdivision) <- unlist(lapply( names(C_i_subdivision), function(j) .p( i, j ) ) )
		
	subdivide_cluster( hclustering_list=C_expanded,
			      cluster_tuple=i,
                                subdivision= C_i_subdivision) -> C_expanded #add subclusters to the hclust list

	f_i = left2explain/noiseless_C_i_var #target fraction of variance for the reconstruction of E1 in C_i
	blockwise_PCA_reconstruction(E1,E1_variances, clustering_list= C_i_subdivision,
				      f=f_i, 
					k=k) -> C_i.blockPCA	
	print( 'left2explain in cluster: ')
	print(left2explain)
	print(' new reconstruction var: ')	
	X_new <- X
	X_new[,]=0
	allocate_blocks( where=X_new, blockwise_PCA= C_i.blockPCA)  -> X_new
	print(sum(apply(X_new,2,var)))
	if ( length(C_i) < mc_beforeStop )
		for (bl in 1:length(C_i.blockPCA))
			C_i.blockPCA[[bl]]$exhausted=TRUE
	
	C_expanded_blockPCA<- c( C_expanded_blockPCA, C_i.blockPCA)
	}
return( list( hclustering_list= C_expanded,
	      hcluster_blockPCA= C_expanded_blockPCA) )
}








##TO DO
fill_in_noise_variables<- function ( final_reconstruction_container, base_clustering ) {
invisible()
}


###############################################################################################################
######SUBDIVISION IN GENERAL###################################################################################
###############################################################################################################


p.<-function(str_tuple) unlist(strsplit(str_tuple, split=','))

get_max_g<- function(hclustering_list){
if (length(hclustering_list)==0) stop(" hclustering_list must be nonempty") 
names(hclustering_list)-> all_tuples
 tuple_lengths<-unlist( lapply(all_tuples, function(tpl) length(p.(tpl)) ))
 max(tuple_lengths)
}


get_partition_at_g<- function( hclustering_list,
			      g){
 stopifnot(g>0)
if (length(hclustering_list)==0) stop(" hclustering_list must be nonempty") 
names(hclustering_list)-> all_tuples
 tuple_lengths<-unlist( lapply(all_tuples, function(tpl) length(p.(tpl)) ))
 stopifnot( g <= max(tuple_lengths))
 g_tuples<- all_tuples[ tuple_lengths == g ] 
 g_P<-hclustering_list [ g_tuples ] 
 attr(g_P,'params')<- attr(hclustering_list, 'params')
 g_P
}




subcluster_and_blockwisePCA<- function(X, X_similarity_matrix, 
				X_variances, hclustering_list, hcluster_PCA_blocks, k, f,
				clfun1, clfun2, min_cl_size, mc_beforeStop=1.5*min_cl_size,
				 clfun1OtherArgs_functionalDependence=NULL,
				 clfun1OtherArgs_constant=NULL,
				 clfun1OtherArgs_ranges=NULL, 
				 clfun2OtherArgs_constant=NULL,
				  clfun2OtherArgs_ranges=NULL){


Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order
g<- get_max_g(hclustering_list)

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_PCA_blocks [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
}

Xg_variances<- apply(Xg, 2, var)
Eg = X - Xg #residuals
Eg_variances<- apply(Eg, 2, var)
 g_partition<- get_partition_at_g( hclustering_list, g)
 g_blockPCA<- hcluster_PCA_blocks[ names(g_partition) ] 
for (bl in 1:length(g_blockPCA)) #if 'exhausted' status was not set for a block, declare it to be FALSE
	if (is.null(g_blockPCA[[bl]]$exhausted)) g_blockPCA[[bl]]$exhausted = FALSE
#theese lists will be expanded by additional blocks/subdivisions
hcl_expanded<- hclustering_list
hcl_expanded_blockPCA<- hcluster_PCA_blocks


unsaturated_tuples<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) !x$saturated) )) ]

exhausted_clusters<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) x$exhausted) )) ]
unsaturated_tuples<- setdiff ( unsaturated_tuples, exhausted_clusters)    

print(sprintf(" number of unsaturated tuples = %d ", length(unsaturated_tuples) ) )
if (!length(unsaturated_tuples)) stop(' all clusters are either exhausted, or unsaturated. stopping')
for (K in unsaturated_tuples) {
	C_K <-  g_partition[[K]]
        S_Eg.C_K <- cor ( Eg [ ,C_K ] )^2 #building similarity of residuals E1 in cluster C_i
	S_Eg.C_K[ X_similarity_matrix[C_K, C_K] == 0 ]= 0 #carry zeros from the similarity of X
	
        left2explain= f* sum(X_variances[ C_K ] ) - sum(Xg_variances[ C_K] )


	###################################################	
	##Clustering choice################################
	###################################################
	if ( length(C_K) < mc_beforeStop) {  # case of exhausted clusters, only 1 one split remains for them,
	# blockPCA entries are later marked as exhausted=TRUE
        	S_Eg.C_K <- cor ( Eg [ ,C_K ] )^2 #building similarity of residuals E1 in cluster C_i
		clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant, list(n_group=2) )
		clustering_vector2list(do.call(clfun2, clf2Args))-> C_K_subdivision
		cks_attrs<- attributes(C_K_subdivision)
		lapply(C_K_subdivision, function(C_K.j)  C_K[ C_K.j ] ) -> C_K_subdivision #translate local to global indexes
		attributes(C_K_subdivision)<- cks_attrs

	} else { sizeDistrState= 9999; n_tries= length(clfun1OtherArgs_ranges[[1]]); i_try=1
		 clfun1args= clfun1OtherArgs_functionalDependence(g=g, cl_tuple=K )
		 clfun1args<- c(clfun1args, clfun1OtherArgs_constant) 
		while ((i_try < n_tries + 1) && (sizeDistrState>0)) {
			for (pname in names(clfun1OtherArgs_ranges))
				{
				clfun1args[[pname]] = clfun1OtherArgs_ranges[[pname]][[i_try]]
				}
			print(sprintf("trying cluster tuple %s with settings:", K))
			print(clfun1args[ names(clfun1args) %in% names(clfun1OtherArgs_ranges) ] )	
			cluster_and_separate_noise( similarity_matrix= S_Eg.C_K, 
					   clustering_function=clfun1, 
					   min_cl_size=min_cl_size, 
					   clustering_functionOtherArgs= clfun1args) -> C_K_subdivision; cks_attrs<- attributes(C_K_subdivision)
			
			lapply(C_K_subdivision, function(C_K.j)  C_K[ C_K.j ] ) -> C_K_subdivision #translate local to global indexes
			attributes(C_K_subdivision)<- cks_attrs
			assess_size_distribution( clustering_list=C_K_subdivision,
						varLeft2explain=left2explain, residuals=Eg)-> sizeDistrState 
			i_try=i_try+1
			print(sprintf("sizeDistrState= %d", sizeDistrState))
			print(sprintf("itry = %d", i_try))
			print(sprintf("n_tries = %d", n_tries))
			print((i_try < n_tries + 1) && (sizeDistrState>0) )
		}

		if (sizeDistrState==0) {   # if TRUE, then  C_K_subdivision has one element and 
					   # and also, this one cluster has enough variance to meet the left2explain target value
					   # therefore we have to split it using clfun2 instead of clfun1
					   # we switch the domain to subdivide from C_K to C_{K=/=K0}
					   C_K<- C_K_subdivision[[1]]  
        				   S_Eg.C_K <- cor ( Eg [ ,C_K ] )^2  #have to recalculate similarity matrix
				#	    S_Eg.C_K[ X_similarity_matrix[C_K, C_K] == 0 ]= 0 #carry zeros from the similarity of X
					   # check: maybe C_{i=/=0} has less that mc_beforeStop elements...
					   if ( length(C_K) < mc_beforeStop) { # if yes, use clfun2 as a last splitting step and mark the subdivision of C_K as exhausted later 
						    clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant, list(n_group=2) )
						    clustering_vector2list(do.call(clfun2, clf2Args))-> C_K_subdivision
						    cks_attrs<- attributes(C_K_subdivision)
						    lapply(C_K_subdivision, function(C_K.j)  C_K[ C_K.j ] ) -> C_K_subdivision #translate local to global indexes
						    attributes(C_K_subdivision)<- cks_attrs
					    } else {  # otherwise C_{K=/=K0} has enough elements, so we can apply clfun2 with normal parameters, so set sizeDistrState <0     
						     sizeDistrState= 1
					    }
				       }	
		if (sizeDistrState>0) {    #in that case, clfun1 did not work (too small clusters/ too big noise part)
					    #have to use clfun2 ("percollator") instead (without noise / non-noise distinction)
					    #either on the whole initial C_K or substituted for C_{K=/=K0} in previous sizeDistrState==0 clause.
			print(sprintf('cluster tuple %s getting processed by clfun2', K))
			clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant )
			settings_scores<-vector()
			candidate_cl_lists<-list()
			settings_m<- length(clfun2OtherArgs_ranges[[1]])
			for (i_settings in 1:settings_m) {
				for (pname in names(clfun2OtherArgs_ranges))
					clf2Args[[pname]] = clfun2OtherArgs_ranges[[pname]][[i_settings]]	
				clustering_vector2list(do.call(clfun2, clf2Args))-> candidate_cl
				can_attrs<- attributes(candidate_cl)
				lapply(candidate_cl, function(C_K.j)  C_K[ C_K.j ] ) -> candidate_cl #translate local to global indexes
				attributes(candidate_cl)<- can_attrs
				candidate_cl_lists[[ i_settings ]] = candidate_cl
				settings_scores [[ i_settings ]] =  normalized_cl_entropy( candidate_cl )     # "best" precollating clustering is the one which proportions are closest
													      # ... to groups of even size
				}
			print('picked best')
			C_K_subdivision<- candidate_cl_lists [[ which.max( settings_scores) ]]	
				       }
	} # end else from if (size(C_K) < minimalsizeBeforeStop)

	#######################################################	
	##Clustering choice END################################
	#######################################################
	#at this point, whatever subdivision happened, the result is sensible

	noiseless_C_K_var<- sum(apply(  Eg[,unlist(C_K_subdivision)], 2, var ))
	# rename each cluster 'j' in subdivision from 'j' to 'i1,i2..ig,j where K= i1,i2...ig' 
	names(C_K_subdivision) <- unlist(lapply( names(C_K_subdivision), function(j) .p( K, j ) ) )
	print(lapply(C_K_subdivision, length))
	subdivide_cluster( hclustering_list=hcl_expanded,
			      cluster_tuple=K,
                                subdivision= C_K_subdivision) -> hcl_expanded #add subclusters to the hclust list

	f_K = left2explain/noiseless_C_K_var #target fraction of variance for the reconstruction of E1 in C_i
	blockwise_PCA_reconstruction(Eg,Eg_variances, clustering_list= C_K_subdivision,
				      f=f_K, 
					k=k) -> C_K.blockPCA	
	if ( length(C_K) < mc_beforeStop )
		for (bl in 1:length(C_K.blockPCA))
			C_K.blockPCA[[bl]]$exhausted=TRUE
 	hcl_expanded_blockPCA<- c( hcl_expanded_blockPCA, C_K.blockPCA)	
	}
return( list( hclustering_list= hcl_expanded,
	      hcluster_blockPCA= hcl_expanded_blockPCA) )
}

get_reconstruction_at_g<- function(X,hcluster_blockPCA, hclustering_list, g) {

Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_blockPCA [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
}
return(Xg)
}




subcluster_and_blockwisePCA2<- function(X, X_similarity_matrix, 
				X_variances, hclustering_list, hcluster_PCA_blocks, k, f,
				clfun2, min_cl_size, mc_beforeStop=1.5*min_cl_size,
				 clfun2OtherArgs_constant=NULL,
				  clfun2OtherArgs_ranges=NULL){


Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order
g<- get_max_g(hclustering_list)

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_PCA_blocks [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
}
Xg_variances<- apply(Xg, 2, var)
Eg = X - Xg #residuals
Eg_variances<- apply(Eg, 2, var)
 g_partition<- get_partition_at_g( hclustering_list, g)
 g_blockPCA<- hcluster_PCA_blocks[ names(g_partition) ] 
for (bl in 1:length(g_blockPCA)) #if 'exhausted' status was not set for a block, declare it to be FALSE
	if (is.null(g_blockPCA[[bl]]$exhausted)) g_blockPCA[[bl]]$exhausted = FALSE
#theese lists will be expanded by additional blocks/subdivisions
hcl_expanded<- hclustering_list
hcl_expanded_blockPCA<- hcluster_PCA_blocks


unsaturated_tuples<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) !x$saturated) )) ]

exhausted_clusters<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) x$exhausted) )) ]
unsaturated_tuples<- setdiff ( unsaturated_tuples, exhausted_clusters)    
#touched_tuples<- vector()
print(sprintf("in P_%d, number of unsaturated tuples = %d ",g, length(unsaturated_tuples) ) )
if (!length(unsaturated_tuples)) stop(' all clusters are either exhausted, or unsaturated. stopping')
for (K in unsaturated_tuples) {
	C_K <-  g_partition[[K]]
        left2explain= f* sum(X_variances[ C_K ] ) - sum(Xg_variances[ C_K] )
	if (left2explain>0){ 
#	touched_tuples=c(touched_tuples, K)
	C_K.Eg_varOrder<-order(Eg_variances[ C_K ], decreasing=TRUE)
	cumul_vars<- cumsum(Eg_variances[C_K][ C_K.Eg_varOrder ])
	min(which(cumul_vars > left2explain))-> up2which
	C_K [ C_K.Eg_varOrder ] [1:up2which ] -> sufficient_subset	
#	print(length(C_K))
	C_K <- sufficient_subset
#	print(length(sufficient_subset))
#	print( sum(apply(Eg[,C_K, drop=FALSE],2,var)))
#	print( left2explain)
        S_Eg.C_K <- cor ( Eg [ ,C_K, drop=FALSE ] )^2 #building similarity of residuals E1 in cluster C_i
	###################################################	
	##Clustering choice################################
	###################################################
	if ( length(C_K) < mc_beforeStop) {  # case of exhausted clusters, only 1 one split remains for them,
	# blockPCA entries are later marked as exhausted=TRUE
		if (length(C_K) > k) {
		clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant, list(n_group=2) )
		clustering_vector2list(do.call(clfun2, clf2Args))-> C_K_subdivision
		cks_attrs<- attributes(C_K_subdivision)
		lapply(C_K_subdivision, function(C_K.j)  C_K[ C_K.j ] ) -> C_K_subdivision #translate local to global indexes
		attributes(C_K_subdivision)<- cks_attrs} else { C_K_subdivision<- list(C_K);
								attributes(C_K_subdivision)<- list(params=list(n_group=1)) 
								names(C_K_subdivision) = c( '1')
							       }

	} else { 
		print(sprintf('cluster tuple %s getting processed by clfun2', K))
			clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant )
			settings_scores<-vector()
			candidate_cl_lists<-list()
			settings_m<- length(clfun2OtherArgs_ranges[[1]])
			for (i_settings in 1:settings_m) {
				for (pname in names(clfun2OtherArgs_ranges))
					clf2Args[[pname]] = clfun2OtherArgs_ranges[[pname]][[i_settings]]	
				clustering_vector2list(do.call(clfun2, clf2Args))-> candidate_cl
				can_attrs<- attributes(candidate_cl)
				lapply(candidate_cl, function(C_K.j)  C_K[ C_K.j ] ) -> candidate_cl #translate local to global indexes
				attributes(candidate_cl)<- can_attrs
				candidate_cl_lists[[ i_settings ]] = candidate_cl
				settings_scores [[ i_settings ]] =  normalized_cl_entropy( candidate_cl )     # "best" precollating clustering is the one which proportions are closest
													      # ... to groups of even size
				}
			print('picked best')
			C_K_subdivision<- candidate_cl_lists [[ which.max( settings_scores) ]]	
		}

	#######################################################	
	##Clustering choice END################################
	#######################################################
	#at this point, whatever subdivision happened, the result is sensible

	if (is.null(names(C_K_subdivision))) stop('null names before change!')
	noiseless_C_K_var<- sum(apply(  Eg[,unlist(C_K_subdivision), drop=FALSE], 2, var ))
	# rename each cluster 'j' in subdivision from 'j' to 'i1,i2..ig,j where K= i1,i2...ig' 
	names(C_K_subdivision) <- unlist(lapply( names(C_K_subdivision), function(j) .p( K, j ) ) )
	#print(lapply(C_K_subdivision, length))
	subdivide_cluster( hclustering_list=hcl_expanded,
			      cluster_tuple=K,
                                subdivision= C_K_subdivision) -> hcl_expanded #add subclusters to the hclust list

	f_K = left2explain/noiseless_C_K_var #target fraction of variance for the reconstruction of E1 in C_i
	if (left2explain > noiseless_C_K_var) stop('something is wrong, left over variance of residuals is smaller than the left2explain value')
	if (is.null(names(C_K_subdivision))) stop('null names!')
	blockwise_PCA_reconstruction(Eg,Eg_variances, clustering_list= C_K_subdivision,
				      f=f_K, 
					k=k) -> C_K.blockPCA	
	if ( length(C_K) < mc_beforeStop )
		for (bl in 1:length(C_K.blockPCA))
			C_K.blockPCA[[bl]]$exhausted=TRUE
 	hcl_expanded_blockPCA<- c( hcl_expanded_blockPCA, C_K.blockPCA)	}
	}
return( list( hclustering_list= hcl_expanded, #touched_tuples=touched_tuples,
	      hcluster_blockPCA= hcl_expanded_blockPCA) )
}

get_reconstruction_at_g<- function(X,hcluster_blockPCA, hclustering_list, g) {

Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_blockPCA [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
}
return(Xg)
}
subcluster_and_blockwisePCA3<- function(X, X_similarity_matrix, 
				X_variances, hclustering_list, hcluster_PCA_blocks, k, f,
				clfun2,# min_cl_size, mc_beforeStop=1.5*min_cl_size,
				 clfun2OtherArgs_constant=NULL,
				  clfun2OtherArgs_ranges=NULL){


Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order
g<- get_max_g(hclustering_list)

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_PCA_blocks [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
}
Xg_variances<- apply(Xg, 2, var)
Eg = X - Xg #residuals
Eg_variances<- apply(Eg, 2, var)
 g_partition<- get_partition_at_g( hclustering_list, g)
 g_blockPCA<- hcluster_PCA_blocks[ names(g_partition) ] 
for (bl in 1:length(g_blockPCA)) #if 'exhausted' status was not set for a block, declare it to be FALSE
	if (is.null(g_blockPCA[[bl]]$exhausted)) g_blockPCA[[bl]]$exhausted = FALSE
#theese lists will be expanded by additional blocks/subdivisions
hcl_expanded<- hclustering_list
hcl_expanded_blockPCA<- hcluster_PCA_blocks


unsaturated_tuples<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) !x$saturated) )) ]

exhausted_clusters<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) x$exhausted) )) ]
unsaturated_tuples<- setdiff ( unsaturated_tuples, exhausted_clusters)    
#touched_tuples<- vector()
print(sprintf("in P_%d, number of unsaturated tuples = %d ",g, length(unsaturated_tuples) ) )
if (!length(unsaturated_tuples)) stop(' all clusters are either exhausted, or unsaturated. stopping')
for (K in unsaturated_tuples) {
	C_K <-  g_partition[[K]]
	print("haveC_K")
        left2explain= f* sum(X_variances[ C_K ] ) - sum(Xg_variances[ C_K] )
	if (left2explain>0){ 
#	touched_tuples=c(touched_tuples, K)
#	print(length(C_K))
#	print(length(sufficient_subset))
#	print( sum(apply(Eg[,C_K, drop=FALSE],2,var)))
#	print( left2explain)
        S_Eg.C_K <- cor ( Eg [ ,C_K, drop=FALSE ] )^2 #building similarity of residuals E1 in cluster C_i
	###################################################	
	##Clustering choice################################
	###################################################
	# blockPCA entries are later marked as exhausted=TRUE
		if (length(C_K) <= k) 
		{ C_K_subdivision<- list(C_K);
		attributes(C_K_subdivision)<- list(params=list(n_group=1)) 
		names(C_K_subdivision) = c( '1')
		} else {

		clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant )
		settings_scores<-vector()
		candidate_cl_lists<-list()
		settings_m<- length(clfun2OtherArgs_ranges[[1]])
		for (i_settings in 1:settings_m) {
			for (pname in names(clfun2OtherArgs_ranges))
				clf2Args[[pname]] = clfun2OtherArgs_ranges[[pname]][[i_settings]]	
			clustering_vector2list(do.call(clfun2, clf2Args))-> candidate_cl
			can_attrs<- attributes(candidate_cl)
			lapply(candidate_cl, function(C_K.j)  C_K[ C_K.j ] ) -> candidate_cl #translate local to global indexes
			attributes(candidate_cl)<- can_attrs
			candidate_cl_lists[[ i_settings ]] = candidate_cl
			settings_scores [[ i_settings ]] =  normalized_cl_entropy( candidate_cl )     # "best" precollating clustering is the one which proportions are closest
												      # ... to groups of even size
			}
		print('picked best')
		C_K_subdivision<- candidate_cl_lists [[ which.max( settings_scores) ]]	
		}
	#######################################################	
	##Clustering choice END################################
	#######################################################
	#at this point, whatever subdivision happened, the result is sensible

	if (is.null(names(C_K_subdivision))) stop('null names before change!')
	noiseless_C_K_var<- sum(apply(  Eg[,unlist(C_K_subdivision), drop=FALSE], 2, var ))
	# rename each cluster 'j' in subdivision from 'j' to 'i1,i2..ig,j where K= i1,i2...ig' 
	names(C_K_subdivision) <- unlist(lapply( names(C_K_subdivision), function(j) .p( K, j ) ) )
	#print(lapply(C_K_subdivision, length))
	subdivide_cluster( hclustering_list=hcl_expanded,
			      cluster_tuple=K,
                                subdivision= C_K_subdivision) -> hcl_expanded #add subclusters to the hclust list

	f_K = left2explain/noiseless_C_K_var #target fraction of variance for the reconstruction of E1 in C_i
	if (left2explain > noiseless_C_K_var) stop('something is wrong, left over variance of residuals is smaller than the left2explain value')
	if (is.null(names(C_K_subdivision))) stop('null names!')
	blockwise_PCA_reconstruction(Eg,Eg_variances, clustering_list= C_K_subdivision,
				      f=f_K, 
					k=k) -> C_K.blockPCA	
	if ( length(C_K) < k )
		for (bl in 1:length(C_K.blockPCA))
			C_K.blockPCA[[bl]]$exhausted=TRUE
 	hcl_expanded_blockPCA<- c( hcl_expanded_blockPCA, C_K.blockPCA)	
	}}
return( list( hclustering_list= hcl_expanded, #touched_tuples=touched_tuples,
	      hcluster_blockPCA= hcl_expanded_blockPCA) )
}


subcluster_and_blockwisePCA4<- function(X, X_similarity_matrix, 
				X_variances, hclustering_list, hcluster_PCA_blocks, k, f,
				clfun2,# min_cl_size, mc_beforeStop=1.5*min_cl_size,
				 clfun2OtherArgs_constant=NULL,
				  clfun2OtherArgs_ranges=NULL){


Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order
g<- get_max_g(hclustering_list)
gj_vars<-rep(0,ncol(X))
for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_PCA_blocks [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
	gj_vars=gj_vars+apply(X_g_j,2,var)
}
Xg_variances<- apply(Xg, 2, var)
print('variances of sums:')
print(sum(Xg_variances))
print('sum of variances:')
print(sum(gj_vars))
Eg = X - Xg #residuals
Eg_variances<- apply(Eg, 2, var)
 g_partition<- get_partition_at_g( hclustering_list, g)
 g_blockPCA<- hcluster_PCA_blocks[ names(g_partition) ] 
for (bl in 1:length(g_blockPCA)) #if 'exhausted' status was not set for a block, declare it to be FALSE
	if (is.null(g_blockPCA[[bl]]$exhausted)) g_blockPCA[[bl]]$exhausted = FALSE
#theese lists will be expanded by additional blocks/subdivisions
hcl_expanded<- hclustering_list
hcl_expanded_blockPCA<- hcluster_PCA_blocks


unsaturated_tuples<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) !x$saturated) )) ]

exhausted_clusters<- names(g_blockPCA) [which( unlist(lapply(g_blockPCA, function(x) x$exhausted) )) ]
unsaturated_tuples<- setdiff ( unsaturated_tuples, exhausted_clusters)    
#touched_tuples<- vector()
print(sprintf("in P_%d, number of unsaturated tuples = %d ",g, length(unsaturated_tuples) ) )
if (!length(unsaturated_tuples)) stop(' all clusters are either exhausted, or unsaturated. stopping')

for (K in unsaturated_tuples) {
	C_K <-  g_partition[[K]]
        left2explain= f* sum(X_variances[ C_K ] ) - sum(Xg_variances[ C_K] )
	stopifnot( sum(Eg_variances[C_K]) > left2explain)
	if (left2explain>0){ 
#	touched_tuples=c(touched_tuples, K)
	Eg.C_K_varOrder<- order(Eg_variances [ C_K ], decreasing= TRUE)
	firstSufficient<-min(which(cumsum(Eg_variances[Eg.C_K_varOrder]	) > left2explain ))
	C_K[Eg.C_K_varOrder][1:firstSufficient]->C_K
	stopifnot( sum(Eg_variances[C_K]) > left2explain)
        S_Eg.C_K <- cor ( Eg [ ,C_K, drop=FALSE ] )^2 #building similarity of residuals E1 in cluster C_i
	###################################################	
	##Clustering choice################################
	###################################################
	# blockPCA entries are later marked as exhausted=TRUE
		if (length(C_K) <= k) 
		{ C_K_subdivision<- list(C_K);
		attributes(C_K_subdivision)<- list(params=list(n_group=1)) 
		names(C_K_subdivision) = c( '1')
		} else {

		clf2Args<- c(list( S_Eg.C_K), clfun2OtherArgs_constant )
		settings_scores<-vector()
		candidate_cl_lists<-list()
		settings_m<- length(clfun2OtherArgs_ranges[[1]])
		for (i_settings in 1:settings_m) {
			for (pname in names(clfun2OtherArgs_ranges))
				clf2Args[[pname]] = clfun2OtherArgs_ranges[[pname]][[i_settings]]	
			clustering_vector2list(do.call(clfun2, clf2Args))-> candidate_cl
			can_attrs<- attributes(candidate_cl)
			lapply(candidate_cl, function(C_K.j)  C_K[ C_K.j ] ) -> candidate_cl #translate local to global indexes
			attributes(candidate_cl)<- can_attrs
			candidate_cl_lists[[ i_settings ]] = candidate_cl
			settings_scores [[ i_settings ]] =  normalized_cl_entropy( candidate_cl )     # "best" precollating clustering is the one which proportions are closest
												      # ... to groups of even size
			}
		print('picked best')
		C_K_subdivision<- candidate_cl_lists [[ which.max( settings_scores) ]]	
		}
	#######################################################	
	##Clustering choice END################################
	#######################################################
	#at this point, whatever subdivision happened, the result is sensible

	if (is.null(names(C_K_subdivision))) stop('null names before change!')
	noiseless_C_K_var<- sum(apply(  Eg[,unlist(C_K_subdivision), drop=FALSE], 2, var ))
	# rename each cluster 'j' in subdivision from 'j' to 'i1,i2..ig,j where K= i1,i2...ig' 
	names(C_K_subdivision) <- unlist(lapply( names(C_K_subdivision), function(j) .p( K, j ) ) )
	#print(lapply(C_K_subdivision, length))
	subdivide_cluster( hclustering_list=hcl_expanded,
			      cluster_tuple=K,
                                subdivision= C_K_subdivision) -> hcl_expanded #add subclusters to the hclust list

	f_K = left2explain/noiseless_C_K_var #target fraction of variance for the reconstruction of E1 in C_i
	if (left2explain > noiseless_C_K_var) stop('something is wrong, left over variance of residuals is smaller than the left2explain value')
	if (is.null(names(C_K_subdivision))) stop('null names!')
	blockwise_PCA_reconstruction(Eg,Eg_variances, clustering_list= C_K_subdivision,
				      f=f_K, 
					k=k) -> C_K.blockPCA
	print( 'left2explain in cluster: ')
	print(left2explain)
	print(' new reconstruction var: ')	
	X_new <- X
	X_new[,]=0
	allocate_blocks( where=X_new, blockwise_PCA= C_K.blockPCA)  -> X_new
	print(sum(apply(X_new,2,var)))
	if ( length(C_K) < k )
		for (bl in 1:length(C_K.blockPCA))
			C_K.blockPCA[[bl]]$exhausted=TRUE
 	hcl_expanded_blockPCA<- c( hcl_expanded_blockPCA, C_K.blockPCA)	
	}}
return( list( hclustering_list= hcl_expanded, #touched_tuples=touched_tuples,
	      hcluster_blockPCA= hcl_expanded_blockPCA) )
}

###############################################################################################################
######FINISHING TOUCHES########################################################################################
###############################################################################################################

add_noise<- function( Xg, X) {

Eg <- X - Xg
Eg_variances<- apply(Eg,2,var)

for (j in 1:length(Eg_variances) )
	Xg[,j] = Xg[,j] + rnorm( nrow(Eg), 0, sqrt(Eg_variances[[j]]) )
return(Xg)
}

get_full_reconstruction<- function(X,hcluster_blockPCA, hclustering_list, add_noise=TRUE, uncenter=TRUE) {

g<-get_max_g(hclustering_list)
print(g)
Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	allocate_blocks( where=X_g_j, blockwise_PCA= hcluster_blockPCA [ names(g_j_clustering) ] )-> X_g_j
	Xg = Xg + X_g_j
}
if (add_noise) Xg<- add_noise(Xg, X)
if (uncenter) { colMeans(X)-> meansX
		for(j in 1:ncol(X)) Xg[,j] = Xg[,j] + meansX[[j]]
	      }
return(Xg)
}


###############################################################################################################
######EVALUATION###############################################################################################
###############################################################################################################

SE<- function(x,y) { (x-y)^2 }
MSE<- function(x,y) { mean(SE(x,y)) }


ccfast<- function(A)    #optimized WGCNA clustering coefficient
      {
      stopifnot(sum(diag(A)==0)==nrow(A))
      N<-nrow(A)
      A.A<-A %*% A
      numer<- vector()
      for(i in 1:N)
          numer[[i]]<- A[i,] %*% A.A[i,]
      denom= colSums(A)^2 - colSums(A^2)
      ifelse(denom==0,0,numer/denom)
      }

ecdf_distance<- function(y, ecdf_x, summarizeFun=NULL) {
ecdf_y<- ecdf(y)
 distt=abs( ecdf_x(y) - ecdf_y(y) )
 if (!is.null(summarizeFun) ) distt=summarizeFun( distt) 
return(distt)
}



###############################################################################################################
######SIMULATION###############################################################################################
###############################################################################################################

PC_generator_matrix<- function(hcluster_blockPCA){

all_PC<- do.call( cbind, lapply(hcluster_blockPCA, function(block) block$PC) )
colnames(all_PC)<-do.call(c, lapply( seq_along(hcluster_blockPCA), function(i) 
					rep( names(hcluster_blockPCA)[[i]],
					     hcluster_blockPCA[[i]]$k_used)
				    )
			)
all_PC
}

multivariate_normal_cholesky<- function( sigma, n_samples){
U = chol(sigma)
t(U) -> L
Z<- matrix(rnorm( ncol(sigma)*n_samples), nrow= n_samples )
t( L %*% t(Z) )
}

#this function simulates original PC_i as Z_i from normal distribution, perserving cov(PC_i,PC_j)
cloned_normalPC_matrix<- function(reference_PC_gen_matrix) {
	PCsig<- cov(reference_PC_gen_matrix)
	PCnorm<-multivariate_normal_cholesky(sigma=PCsig, n_samples=nrow(reference_PC_gen_matrix) ) 
	colnames(PCnorm)<- colnames(reference_PC_gen_matrix)
	PCnorm
}

#below : requires rmetalog

#this function maps PC_i -> H_i, where H_i= (qnorm o F_i) ( PC_i),
# where F_i is the cdf of the metalog distribution fitted to PC_i, 
# and qnorm is quantile function of standard normal distribution

hidden_normalPC<- function(ref_PC, target_metalog, term) {

qnorm(  pmetalog( m=target_metalog, q=ref_PC, term=term)  , 0,1)

}

hidden_normalPC_matrix<- function(reference_PC_gen_matrix, target_metalogs, target_terms) {
stopifnot( ncol(PCnormal_matrix) == length(target_metalogs) )
stopifnot( length(target_metalogs) == length(target_terms) )
PChidden_normal<- reference_PC_gen_matrix

for (j in 1:ncol(reference_PC_gen_matrix))
	PChidden_normal[,j] <- hidden_normalPC( ref_PC= reference_PC_gen_matrix[,j],
						target_metalog= target_metalogs[[j]],
						term= target_terms[[j]] )

PChidden_normal
}

pull2targetMetalog<- function(x_normal, target_metalog_distr, term) {
 u<- pnorm(x_normal, mean= 0, sd=1 )
qmetalog( target_metalog_distr, y= u, term=term) 
}

cloned_metalogPC_matrix<- function( PChidden_normal, target_metalogs, target_terms){

stopifnot( ncol(PChidden_normal) == length(target_metalogs) )
stopifnot( length(target_metalogs) == length(target_terms) )
H_cov= cov(PChidden_normal)
H_clones<-multivariate_normal_cholesky(sigma=H_cov, n_samples=nrow(PChidden_normal) ) 
H_clones<- scale(H_clones, center=TRUE, scale=TRUE)
	colnames(H_clones)<- colnames(PChidden_normal)
PCmeta<- H_clones
for (j in 1:ncol(PCmeta))
	PCmeta[,j] = pull2targetMetalog(x_normal=H_clones[,j],target_metalog_distr=target_metalogs[[j]],
					term=target_terms[[j]]
					)
PCmeta 
}

add_noise2<- function( Xg, Eg_variances) {


for (j in 1:length(Eg_variances) )
	Xg[,j] = Xg[,j] + rnorm( nrow(Xg), 0, sqrt(Eg_variances[[j]]) )
return(Xg)
}

get_cloned_reconstruction<- function(X,hcluster_blockPCA, hclustering_list, cloned_PCmatrix, noise_variances=NULL, uncenter=TRUE) {

g<-get_max_g(hclustering_list)
print(g)
Xg<- X
Xg[,]=0  #placeholder for the reconstruction of gth order

for (g_j in 1:g) {
	get_partition_at_g( hclustering_list, g_j) -> g_j_clustering
	X_g_j <- X
	X_g_j[,]=0
	hcluster_blockPCA[ names(g_j_clustering) ]-> g_j_blocks
	names(g_j_blocks) <- names(g_j_clustering)
	for ( K in names(g_j_blocks) )
		g_j_blocks[[K]]$reconstruction <- cloned_PCmatrix[, colnames(cloned_PCmatrix)==K,drop=FALSE ] %*% g_j_blocks[[K]]$Vt
	allocate_blocks( where=X_g_j, blockwise_PCA= g_j_blocks )-> X_g_j
	Xg = Xg + X_g_j
}
if (!is.null(noise_variances)) Xg<- add_noise2(Xg, noise_variances)
if (uncenter) { colMeans(X)-> meansX
		for(j in 1:ncol(X)) Xg[,j] = Xg[,j] + meansX[[j]]
	      }
return(Xg)
}




