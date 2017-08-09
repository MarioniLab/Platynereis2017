# Platynereis2017
Code for the single cell RNA seq project on identifying cell types in Platynereis dumerilii

This repository contains R scripts to perform quality control, normalization and downstream analysis of scRNAseq data generated from the developing larva of Platynereis dumerilii.

## Preprocessing

This folder contains scripts to annotate the dataset (Gene_annotation.R), to perform quality control (Quality_control) and to test for biological and technical noise in each batch (Batch_variance.R). It also contains a script to normalize the data (Normalization.R) and to remove potential doublets after finding marker genes (Doublet_removal.R).

## Scripts

Contains analysis scripts to perfrom clustering, differential expression and spatial mapping of scRNAseq data as well as the scripts to generate the figures. 

