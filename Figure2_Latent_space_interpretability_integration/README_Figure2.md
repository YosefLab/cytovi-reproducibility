# Figure 2: Latent space interpretability and batch integration
This folder contains the notebooks and scripts to perform the inspection of the latent space with regards to cell type and cell state variability and performs the benchmarking of the integration performance across the two technical replicates of flow cytometry data from Nunez et al. The folder is structured as following:

Notebooks:
- 2_1_Latent_space_interpretability_integration.ipynb: Imports the annodated flow cytometry, mass cytometry and CITE-seq data from Fig 1 and the models trained for the PPCs. Benchmarking of batch correction performance is performed using the two technical replicates from Nunez et al. This notebook also contains the evaluation of feature dropout and runtime analysis.
- 2_2_Cycombine_batch_correction.ipynb: R notebook to perform the batch integration of the Nunez et al. data using cyCombine.
- 2_3_FastMNN_batch_correction.ipynb: R notebook to perform the batch integration of the Nunez et al. data using FastMNN.

Scripts:
- 2_1_Nunez_batch_correction_CytoVI.py: Code to train the CytoVI models for batch correction of the Nunez et al. dataset.
- 2_2_Nunez_Runtime_CytoVI_Harmony.py: Code to perform the runtime analysis of the python-based methods (CytoVI, Harmony).
- 2_3_Nunez_Runtime_FastMNN_Cycombine.R: Code to perform the runtime analysis of the R-based methods (FastMNN, cyCombine).
- 2_4_Nunez_feature_ablation.py: Python code to evaluate CytoVI's sensitivity to feature ablation.
- 2_5_Nunez_hyperparam_grid.py: Python code to train CytoVI models with different configurations for latent dimensions and hidden nodes.

