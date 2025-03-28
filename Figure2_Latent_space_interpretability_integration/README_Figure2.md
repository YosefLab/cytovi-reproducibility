# Figure 2: Latent space interpretability and batch integration
This folder contains the notebooks and scripts to perform the inspection of the latent space with regards to cell type and cell state variability and performs the benchmarking of the integration performance across the two technical replicates of flow cytometry data from Nunez et al. The folder is structured as following:

Notebooks:
- 2_1_Latent_space_interpretability_integration.ipynb: Imports the annodated flow cytometry, mass cytometry and CITE-seq data from Fig 1 and the models trained for the PPCs. Benchmarking of batch correction performance is performed using the two technical replicates from Nunez et al.
- 2_2_Cycombine_batch_correction: R notebook to perform the batch integration of the Nunez et al. data using cyCombine.

