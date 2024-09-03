# Overview
mol_data_prep is comprised of R and Python scripts for streamlining the preparation of molecular datasets for machine learning tasks. It provides methods for descriptor calculation, dataset standardization, response variable binarization, clustering for test set selection, and synthetic data generation using SMOTE.

The main script provides a class with handy methods for preparing datasets of molecules for Machine Learning procedures

- Descriptor calculation using MORDRED (peptide descriptors and fingerprints soon to be added)
- Standardization
- Binarization of response variable
- Allows clustering strategies to build more representative test sets
- Easy introduction of synthetic datapoints through SMOTE (SMOGN soon to be added)

# Usage
In the "example" folder you can find a step by step example of a standard procedure for a classification problem

# Citation

A similar pipeline was originally used in:

Machado, L. A., Krempser, E. and Guimarães, A. C. R - A machine learning-based virtual screening for natural compounds capable of inhibiting the HIV-1 integrase. Frontiers in Drug Discovery (2022)
