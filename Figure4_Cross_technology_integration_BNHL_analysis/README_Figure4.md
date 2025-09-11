# Figure 4: Cross-technology integration and analysis of the BNHL cohort
This folder contains the notebooks and scripts to integrate mass cytometry and CITE-seq data from Levine et al. and Stuart et al. and benhmark CytoVi in relation to competitor methods. Additionally, it contains code to integrate flow and mass cytometry data and impute scatter features in the mass cytometry dataset. It also contains code to reproduce the analysis of the BNHL cohort from Roider et al. The folder is structured as following:

Notebooks:
- 4_1_Cross_technology_integration_benchmarking.ipynb: Notebook to reproduce the benchmarking of the cross-technology integration between Stuart et al. and Levine et al.
- 4_2_Scatter_imputation.ipynb: Notebook to reproduce the imputation of scatter features in mass cytometry data.
- 4_3_BNHL_CITEseq_preprocessing.ipynb: Notebook to read and preprocess the CITE-seq BNHL data from Roider et al.
- 4_4_BNHL_CITEseq_TotalVI.ipynb: Notebook used to run CytoVI and TotalVI on the CITEseq data and obtain batch-corrected expression estimates of RNA and protein expression
- 4_5_BNHL_analysis.ipynb: Notebook to reproduce the analysis of the BNHL cohort and generate relevant figures.

Scripts:
- 4_1_Crossint_Levine_Stuart_CytoVI_Harmony.py: Code to run the CytoVi and Harmony cross-technology integration.
- 4_2_Crossint_Levine_Stuart_FastMNN_Cycombine.R: Code to perform the FastMNN cross-technology integration.
- 4_3_Crossint_Levine_Stuart_CytoVI_Harmony_Noise.py: Script to run CytoVI and Harmony for the noise injection experiment.
- 4_4_Crossint_Levine_Stuart_FastMNN_Noise.R: Script to run FastMNN for the noise injection experiment.
- 4_5_Crossint_Levine_Stuart_Celltype_ablation_CytoVI_Harmony.py: Script to run CytoVI and Harmony for the cell type ablation experiment.
- 4_6_Crossint_Levine_Stuart_Celltype_ablation_FastMNN.R: Script to run FastMNN for the cell type ablation benchmarking.
- 4_7_BNHL_flow_model_training_DA.py: Script to train the flow-flow model of the BNHL cohort and perform the label-free differential abundance analysis.
- 4_8_BNHL_flow_CITEseq_model_training.py: Script to train the flow-CITEseq model of the BNHL cohort and impute transcript expression in flow cytometry data.

