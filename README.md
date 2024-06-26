# Modeling Type 1 Diabetes progression using machine-learning and single-cell transcriptomic measurements in human islets

This repository contains codes for the manuscript titled Modeling Type 1 Diabetes progression using machine-learning and single-cell transcriptomic measurements in human islets.

A schematic workflow of ML-based XGBoost model was built for gene selection and classification.

![pipeline](https://github.com/AbhijeetRPatil/ML_Islets/assets/33159736/7aa0dfc3-5279-4569-b612-f1b77ec92e12)



## R version and packages
R v 4.1.2 was used to develop the models. All the package versions can be found in the manuscript.
Some of the key packages include- 
  ```install.packages("Seurat")```
  ```install.packages("xgboost")```
  ```install.packages("e1071")```
  ```install.packages("foreach")```
  ```install.packages("doParallel")``` and others.

## Data
The processed scRNA-seq dataset (Seurat .RDS object) is uploaded here: https://hpap.pmacs.upenn.edu/analysis.  

## Machine Learning model scripts
The Machine Learning model scripts folder contains the scripts for cell types, unannotated, and LOOCV analysis.

## Libraries
The libraries folder contains the dependencies for the above scripts for generating the data either at cell type level or all cells (unannotated).

## Differential Expression Analysis scripts
The Differential Expression Analysis folder contains the scripts for performing differential expression analysis on either all cells (unannotated) or Pseudobulk.

## Figures and Tables
This folder contains all the scripts used to generate the figures and tables in the manuscript.

## Publication to cite
Patil et al. (2024)  Modeling Type 1 Diabetes progression using machine-learning and single-cell transcriptomic measurements in human islets. Cell Reports Medicine 2024 (Accepted in principle). Link/DOI will be updated soon.
