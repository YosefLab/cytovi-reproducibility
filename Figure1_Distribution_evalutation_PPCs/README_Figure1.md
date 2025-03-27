# Figure 1: Distribution evaluation and Posterior predictive checks (PPCs)
This folder contains the notebooks and scripts to perform the evaluation of distributions of antibody-based single cell data and the evaluation of the generative model using PPCs. The folder is structured as following:

Notebooks:
- 1_1_Preprocessing_annotation_flow_cytometry.ipynb: Imports the flow cytometry data of healthy PBMCs from Nunez et al., performs preprocessing and annotates the normalization samples from two distinct batches separately.
- 1_2_Preprocessing_annotation_mass_cytometry.ipynb: Imports the mass cytometry data of healthy PBMCs from Ingelfinger et al., performs preprocessing and annotates the samples from one individual batch.
- 1_3_Preprocessing_annotation_CITEseq.ipynb: Imports the CITEseq data of PBMCs from Hao et al. and performs preprocessing of one individual sample from one batch.
- 1_4_Distribution_evaluation_PPCs.ipynb: This notebook was utilized to generate the figures and perform the posterior predictive checks of the different technologies and models.

Scripts:
- 1_1_Model_training_flow_cytometry.py: Trains the CytoVI models for flow cytometry data that was used for the PPCs.
- 1_2_Model_training_mass_cytometry.py: Trains the CytoVI models for mass cytometry data that was used for the PPCs.
- 1_3_Model_training_CITEseq.py: Trains the CytoVI models for CITEseq data that was used for the PPCs.
