# NDHS PAF Analysis Script 

This repository contains the R script used to reproduce the NDHS population attributable fraction analysis used for our paper: [Population-attributable burden of modifiable risk factors for depression and anxiety among reproductive-age women in Nepal](https://www.nature.com/articles/s41598-026-43908-8).

## File
- `NDHS PAF analysis.R`

## Input data
The script expects the input file at:

`Data/MergedHHwomen.dta`

The dataset is not included in this repository.

## Important note

This script **does not include variable definitions, data dictionaries, or preprocessing logic explanations**. It assumes prior knowledge of the dataset structure and variable construction.

Researchers attempting to reproduce or adapt this analysis **must define their own variables and modelling logic based on their specific hypotheses**. The script should be treated as a structural reference for the analytical workflow, not a fully self-contained or plug-and-play pipeline.

## Run
Open R in the project folder and run:

```r
source("NDHS PAF analysis.R")
