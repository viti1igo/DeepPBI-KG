---
trigger: always_on
---

# Implementation Plan - WGS vs Key Genes Feature Visualization

This plan details the approach for comparing and visualizing the similarity between whole-genome sequence (WGS) feature vectors and key-gene (KG) feature vectors extracted by the DeepPBI-KG pipeline for phages and bacteria.

## User Review Required

> [!NOTE]
> Since we only have prokka annotations and extracted features for the example dataset (9 phages and 8 bacteria) pre-computed in `template/feature_file/`, the initial analysis will focus on these organisms. 
> To make the analysis robust, we will compute multiple independent metrics:
> 1. Cosine similarity and Pearson correlation of raw and scaled feature vectors for each individual organism.
> 2. PCA dimensionality reduction to 1D and 2D to check for linear relationships between WGS and KG representations.
> 3. Cosine similarity heatmaps of the representations.

## Proposed Changes

We will create a python script `code/visualize_features.py` that loads the pre-computed feature files, performs the dimensionality reduction and correlation analysis, generates high-quality plots, and saves them in the `figure/` directory.

### Visualization Code & Output

#### [NEW] [visualize_features.py](file:///home/axnguyen/Data/DeepPBI-KG/code/visualize_features.py)
A python script to run all analysis and generate figures.

#### [MODIFY] [AGENT.md](file:///home/axnguyen/Data/DeepPBI-KG/.agent/rules/AGENT.md)
The final report containing all visualization methodology, plots, source code, and key insights.

## Verification Plan

### Automated Verification
We will run `visualize_features.py` using python within the activated `DeepPBI-KG` conda environment to produce the figures:
- `figure/phage_wgs_vs_kg_comparison.png`
- `figure/bacterium_wgs_vs_kg_comparison.png`

We will check that the files are created and verify their contents.
