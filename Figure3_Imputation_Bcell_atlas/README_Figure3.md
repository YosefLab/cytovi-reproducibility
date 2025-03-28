# Figure 3: Benchmarking of imputation performance and generation of the B cell maturation atlas
This folder contains the notebooks and scripts to evaluate CytoVIs marker imputation performance, perform cross-study integration and imputation and assemble the B cell maturation atlas from Glass et al. The folder is structured as following:

Notebooks:
- 3_1_Imputation_benchmarking.ipynb: Performs the evaluation of imputation performance of CytoVI including competitor tools.
- 3_2_cyCombine_imputation.ipynb: R notebook to perform the missing marker imputation using cyCombine.
- 3_3_Cross_study_integration_flow_cytometry.ipynb: Notebook to perform the cross-study integration of flow cytometry data from Nunez et al. and Kreutmair et al., including imputation of missing markers.
- 3_4_Bcell_maturation_atlas.ipynb: Imports the mass cytometry data from Glass et al., trains a CytoVI model and imputes missing markers.
- 3_5_Trajectory_inference_Bcells.ipynb: Trajectory inference of the integrated B cell maturation atlas using CellRank.

Scripts:
-3_1_Masking_model_training_imputation.py: Script used to mask markers in the semi-synthetic batches, train CytoVI models and impute missing markers.

