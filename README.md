# spatial-data-integration

A workflow for integrating and comparing spatial transcriptomics data from different platforms (e.g., Visium and Xenium) using R.

This repository contains an R-based analysis that demonstrates how to align and integrate spatial transcriptomics datasets, perform cell‑to‑spot mapping, and visualize multi‑modal spatial molecular data.

## 📌 Overview

Spatial transcriptomics technologies measure gene expression in space, capturing both molecular profiles and tissue architecture. This project shows how to:

- Load spatial data from **10x Genomics Visium** and **Xenium** platforms
- Apply an affine transformation to align the datasets
- Identify neighboring cells per spot using nearest neighbor search
- Aggregate single‑cell data into pseudo‑spots
- Visualize the spatial distributions and compare gene expression across datasets

## 🧠 Key Features

- Imports and preprocesses Visium and Xenium spatial data
- Harmonizes coordinate spaces between platforms
- Uses fixed‑radius neighbor queries for cell mapping
- Generates scatter and overlay plots to assess integration
- Produces integration summary plots (e.g., gene expression comparisons)

## 🧰 Dependencies

This project requires R (≥3.6+) and uses the following packages:

```r
library(RANN)
library(scater)
library(scuttle)
library(harmony)
library(ggspavis)
library(VisiumIO)
library(patchwork)
library(cli)
library(RPostgreSQL)
library(OSTA.data)
library(BayesSpace)
library(Seurat)
library(SpatialExperiment)
library(SpatialExperimentIO)
