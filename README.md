TAPAS Evaluation
This repository contains data and R scripts used for evaluating the uncertainty of LLMs in cell type annotation across multiple tissues.

Repository Structure
📁 /data/
This folder contains input data files used for simulation and evaluation across multiple tissue types.
| File                                                                                                                    | Description                                                                                                                                                                                                            |
| ----------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `GPT4_paper_supp.xlsx`                                                                                                  | Supplementary data from the GPTCelltype paper containing marker gene and cell type annotation information collected from ten datasets.                                                                                 |
| `heart_df.rds`<br>`kidney_df.rds`<br>`lung_df.rds`<br>`pbmc_df.RDS`                                                     | Tissue-specific cell type hierarchy data matrices. Each file represents a cell type relationship tree for a specific tissue.                                                                                           |
| `heart_distance_mtx_na.rda`<br>`kidney_distance_mtx_na.rda`<br>`lung_distance_mtx_na.rda`<br>`pbmc_distance_mtx_na.rda` | Pairwise cell type distance matrices corresponding to the tissue-specific `.rds` files. These matrices also include a synthetic “Unknown” cell type, whose distance is set to 1.5 times the maximum observed distance. |

📁 /R_code/
This folder contains all the R scripts used for simulation, assignment, visualization, and statistical evaluation.
| File                                                                                                        | Description                                                                                                                                               |
| ----------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Heart_predict_assign.R`<br>`Kidney_predict_assign.R`<br>`Lung_predict_assign.R`<br>`PBMC_predict_assign.R` | Scripts to run the prediction and assignment simulations under four different scenarios for each tissue type.                                             |
| `Plot_res.R`                                                                                                | Script to generate the final evaluation figures presented in the paper.                                                                                   |
| `bootstrap.R`                                                                                               | Code to generate the null distribution of heterogeneity scores using a bootstrap method. For details, refer to the **Methods** section of the manuscript. |
| `Utility_func.R`                                                                                            | Helper functions used across the simulation and evaluation pipeline.                                                                                      |
