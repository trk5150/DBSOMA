# SOMSCAN
This project is a java based tool for training a Self Organizing Map based on single cell RNA-seq data which is then used to evaluate gene lists for the degree of cluster formation by the DBSCAN algorithm. 
The goal is to allow users to in silico screen transcriptional perturbation data to make more accurate predictions 

Code herein was a written by Timothy Kunz. The SOM training implementation is a modified from a previous implementation which can be found https://github.com/seqcode/chromosom, with an associate publication https://doi.org/10.1016/j.ymeth.2020.07.002. Other code, including DBSCAN implementation and scripts are original code.

Below is the description of the overall workflow, a guide for using the sample data provided in this repository, and description of the methods used at each step.

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
  
  1) Start with single cell RNA seq data
  
  2) Preprocess the single cell data into 2 files required for SOM training:
  
    a) a pairwise gene x gene correlation matrix based on single cell sequencing data
    
    b) a list of gene names associated with each row/column of the square correlation matrix
  
  3) Train a self organizing map with the pairwise correlation matrix
  
  4) Scan the self organizing map (or maps) for pertubation gene lists which overlap with a set of desired states
  
  5) Interact with the resulting scan outputs to generate a list of perturbations predicted to effect genes of the desired state

  The SOM can also be used to view various gene lists before scanning, which can help to decide on states

## Executable_Jars
Provided are a set of executable jar files which allow for use of the various functionalities described below
Source code for each is contained in this repository, but interaction with the code is not necessary for use
Files:
1) SingleCellParserTSV.jar
2) SparesParse.jar
3) FilePreviewer.jar
4) transposeTSV.jar
5) transposeNFixCSV.jar


## Sample_Files
This directory contains a few files which can be used to test the various functionalities of SOMSCAN repository

Included files:
1) A randomly generated mock of single cell data to test the preprocessing 
2) A randomly generated mock of a correlation matrix and sample gene names to test training
3) A sample bash script for implenting the training
4) The output SOM from actual training on the single cell sequencing data set
   Single cell data aquired from GSE114412
5) Several transcriptional state targets for scanning
6) Several selected perturbation response data sets which overlap with various target states
   Datasets aquired from https://maayanlab.cloud/sigcom-lincs/#/Download
7) An empty directory which will be used to output the scan results
8) A sample bash script for running the in silico scan method

## Preprocessing
The SOM training requires two files
Code in the package "correlationMatrixMaker" can be used to produce files required for SOM Training:
1) SingleCellParserTSV takes a tab delimeted file containing read count data with cell identifiers on columns and genes in rows a pairwise correlation matrix and associated list of gene names
2) SparseParse takes read count data in sparse matrix notation to produce a pairwise correlation matrix and associated list of gene names

Code in the package "file_interaction" can be used to pre-organize single cell data to transpose matricies:
1) transposeNFixCSV = convert a CSV file to tsv file
2) transposeTSV = transpose a TSV file to have cell identifiers on columns and genes in rows
3) FilePreviewer = functionality to read a first 2 lines of a file to ensure the matrix looks as expected


Code in the "CorMatAnalyzer" package can be used to interact with the matrix before beginning training
1) provides functionality to look up a gene or set of genes in a produced pairwise correlation matrix to ensure QC has not removed desired genes

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
MIT License

Copyright (c) 2024 Timothy Kunz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
