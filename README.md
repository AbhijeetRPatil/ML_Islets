# Modeling Type 1 Diabetes progression using machine-learning and single-cell transcriptomic measurements in human islets

## R version and packages
R v 4.1.2 was used to develop the models. All the package versions can be found in the manuscript.
Some of the key packages used are mentioned below- 
  ```install.packages("Seurat")```
  ```install.packages("xgboost")```
  ```install.packages("e1071")```
  ```install.packages("foreach")```
  ```install.packages("doParallel")```

## Data
The processed scRNA-seq dataset (Seurat .RDS object) can be found here: https://hpap.pmacs.upenn.edu/analysis.  

## Differential Expression Analysis scripts
The Differential Expression Analysis folder contains the scripts for performing differential expression analysis on either all cells (unannotated) or Pseudobulk.

## Machine Learning model scripts
The Machine Learning model scripts folder contains the scripts for cell types, unannotated, and LOOCV analysis.

## Libraries
The libraries folder contains the dependencies for the above scripts for generating the data either at cell type level or all cells (unannotated).

## Figures and Tables
This folder contains all the scripts used to generate the figures and tables in the manuscript.

## Publication to cite
Patil et al. (2024)  Modeling Type 1 Diabetes progression using machine-learning and single-cell transcriptomic measurements in human islets. Cell Reports Medicine 2024 (Accepted in principle). Link/DOI will be updated soon.
