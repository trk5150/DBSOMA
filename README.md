# SOMSCAN
This project is a java based tool for training a Self Organizing Map based on single cell RNA-seq data which is then used to evaluate gene lists for the degree of cluster formation by the DBSCAN algorithm. 
The goal is to allow users to in silico screen transcriptional perturbation data to make more accurate predictions 

Code herein was a written by Tim Kunz. The SOM implementation is based on a previous implementation which can be found ___

**In progress**

## Table of Contents
- [Workflow](#Workflow)
- [Executable_Jars](#Executable_Jars)
- [Sample_Files](#Sample_files)
- [Preprocessing](#preprocessing)
- [Training](#Training)
- [Scanning](#Scanning)
- [Viewing](#Viewing)
- [Use](#Use)
- [Supportive_scripts](#Supportive_scripts)
- [License](#license)

## Workflow
The general workflow for use is:
  -Start with single cell RNA seq data
  -Preprocess the single cell data into 2 files required for SOM training:
    -a pairwise gene x gene correlation matrix based on single cell sequencing data
    -a list of gene names associated with each row/column of the square correlation matrix
  -Train a self organizing map with the pairwise correlation matrix
  -Scan the self organizing map (or maps) for pertubation gene lists which overlap with a set of desired states
  -Interact with the resulting scan outputs to generate a list of perturbations predicted to effect genes of the desired state

  The SOM can also be used to view various gene lists before scanning, which can help to decide on states

## Executable_Jars
-Provided are a set of executable jar files which allow for use of the various functionalities described below
-Source code for each is contained in this repository, but interaction with the code is not necessary for use
Files:
-SingleCellParserTSV.jar
-SparesParse.jar
-FilePreviewer.jar
-transposeTSV.jar
-transposeNFixCSV.jar


## Sample_Files
This directory contains a few files which can be used to test the various functionalities of SOMSCAN repository
Included files:
-A randomly generated mock of single cell data to test the preprocessing 
-A randomly generated mock of a correlation matrix and sample gene names to test training
-A sample bash script for implenting the training
-The output SOM from actual training on the single cell sequencing data set
  -Single cell data aquired from GSE114412
-Several transcriptional state targets for scanning
-Several selected perturbation response data sets which overlap with various target states
  -Datasets aquired from https://maayanlab.cloud/sigcom-lincs/#/Download
-An empty directory 
-A sample bash script for running the in silico scan method

## Preprocessing
The SOM training requires two files
Code in the package "correlationMatrixMaker" can be used to produce files required for SOM Training:
-SingleCellParserTSV takes a tab delimeted file containing read count data with cell identifiers on columns and genes in rows a pairwise correlation matrix and associated list of gene names
-SparseParse takes read count data in sparse matrix notation to produce a pairwise correlation matrix and associated list of gene names

Code in the package "file_interaction" can be used to pre-organize single cell data to transpose matricies:
-transposeNFixCSV = convert a CSV file to tsv file
-transposeTSV = transpose a TSV file to have cell identifiers on columns and genes in rows
-FilePreviewer = functionality to read a first 2 lines of a file to ensure the matrix looks as expected


Code in the "CorMatAnalyzer"
-provides functionality to look up a gene or set of genes in a produced pairwise correlation matrix to ensure QC has not removed desired genes

## Training


## Scanning


## Viewing


## Use
A guide to using the provided sample files:


## Supportive_scripts
This repository also contains several python scripts that were used for various data processing and endpoint analyses. However, their general usability is not optimized. I would recommend personalized analyses users based on their exact needs, using these scripts only as guides.
Scripts:
-
-
-

This Repository also contains R scripts used to generate GoTerm analysis plots 
Script:
-

## License

