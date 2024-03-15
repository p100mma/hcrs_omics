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

MCL software:

https://micans.org/mcl/src/mcl-14-137.tar.gz



R packages used:
- Hmisc
- Matrix
- data.table
- rmetalog
- DynamicTreeCut

functions doing most of the work
- `blockwisePCA_R_engine.R`

Experiment I:
- `BRCA_decomposition.R`
- `metalogs.R`
- `BRCA_simulation.R`
- `WGCNA_clustering.R`
- `WGCNA_simulation.R`
- `plain_SVD.R`

Experiment II:
- `BRCA_nPC.R`

TOM plots:
- `tom.R`
- `tomPlots.R`

Summary of results:
- `computation_heavy_metrics.R`
- `nPC_comparison_table.R`
- `sim_comparison_table.R`
- `concat_tables.R`
 
