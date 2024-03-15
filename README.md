# hcrs_omics

## License
The software is licensed under the GNU General Public License, version 3.

## How to run using docker

Run the docker daemon in seprate terminal window if its not set up

	sudo dockerd

Then run following commands

	git clone https://github.com/p100mma/hcrs_omics
	cd hcrs_omics
	sudo docker build -t hcrs .
	sudo docker run -it -v .:/home/hcrs_omics hcrs


## Uses:

MCL software (ver 14-137):

https://micans.org/mcl/src/mcl-14-137.tar.gz

R version 4.3.1

R packages used:
- Hmisc
- Matrix
- data.table
- rmetalog
- DynamicTreeCut

## Input data

`gene_expr_data.rds` - gene expression profiles of 1394 breast cancer patients, 8673 genes.
Original source:
- Pereira B et al. The somatic mutation profiles of 2, 433 breast cancers refine their genomic and transcriptomic landscapes. Nature Communications, 2016; 7.1

Preprocessed and prefiltered as in:
- Polewko-Klim A, Mnich K, Rudnicki W, Robust Data Integration Method for Classification of Biomedical Data. Journal ofMedical Systems, 2021; 45.

## Code

functions doing most of the work
- `blockwisePCA_R_engine.R`

Experiment I (run in that order):
- `BRCA_decomposition.R`
- `metalogs.R`*
- `BRCA_simulation.R`
- `WGCNA_clustering.R`
- `WGCNA_simulation.R`
- `plain_SVD.R`

*`metalogs.R` should be run with command line argument ranging from 1 to 70:

	Rscript --no-save metalogs.R 1
 	#(...)
  	Rscript --no-save metalogs.R 70
   
We have used an HPC system to run this part in parallel.


Experiment II:
- `BRCA_nPC.R` - should be run with command line argument ranging from 2 to 5:
  
        Rscript --no-save BRCA_nPC.R 2
        #(...)
        Rscript --no-save BRCA_nPC.R 5

We have used an HPC system to run this part in parallel.

TOM hierarchical clustering plots (run after executing Experiment I scripts):
- `tom.R`
- `tomPlots.R`

Summary of results (run after Experiment I & II):
- `computation_heavy_metrics.R`
- `nPC_comparison_table.R`
- `sim_comparison_table.R`
- `concat_tables.R`

Topology characteristics plots & KS distances (latter in the supplementary material):
- `topology_plots_KSdistances.R`
 
