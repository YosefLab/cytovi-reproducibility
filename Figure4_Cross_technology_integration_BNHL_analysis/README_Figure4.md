# Figure 4: Cross-technology integration and analysis of the BNHL cohort
This folder contains the notebooks and scripts to integrate flow and mass cytometry data and impute scatter features in the mass cytometry dataset. It also contains code to reproduce the analysis of the BNHL cohort from Roider et al. The folder is structured as following:

Notebooks:
- 4_1_Scatter_imputation.ipynb: Notebook to reproduce the imputation of scatter features in mass cytometry data.
- 4_2_BNHL_CITEseq_preprocessing.ipynb: Notebook to read and preprocess the CITE-seq BNHL data from Roider et al.
- 4_3_BNHL_CITEseq_TotalVI.ipynb: Notebook used to run CytoVI and TotalVI on the CITEseq data and obtain batch-corrected expression estimates of RNA and protein expression
- 4_4_BNHL_analysis.ipynb: Notebook to reproduce the analysis of the BNHL cohort and generate relevant figures.

Scripts:
- 4_1_BNHL_flow_model_training_DA.py: Script train the flow-flow model of the BNHL cohort and perform the label-free differential abundance analysis.
- 4_2_BNHL_flow_CITEseq_model_training.py: Script to train the flow-CITEseq model of the BNHL cohort and impute transcript expression in flow cytometry data.

