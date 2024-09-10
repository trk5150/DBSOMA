# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 14:43:16 2024

@author: tik105
"""

import pandas as pd

file1 = "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//combined_matricies.tsv"
file2 = "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//GSM3722100_YX1L_10X.txt"
file3 = "C://Users//tik105//Desktop//mRNA//Stage 6//GSE114412_Stage_6.all.processed_counts.tsv"
file4 =  "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//GSM3722110_OX3X_10X.txt"
file5 = "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//combined_matricies_and-genes.tsv"
file6 = "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//transposed_from_combined_matricies_and-genes.txt"

# Read the first 5 rows of the tab-delimited file
df = pd.read_csv(file6, sep='\t', nrows=5)
df = df.rename_axis('Barcode', axis=1)

# Select the first 5 columns
subset = df.iloc[:, :10]

print(subset)

#df1 = pd.read_csv(file5, sep='\t')

# Read the second TSV file into another DataFrame
#df2 = pd.read_csv(4, sep='\t')

# Check if the number of rows in both DataFrames is equal
#if df1.shape[0] == df2.shape[0]:
#    print("Both TSV files have the same number of rows.")
#else:
#    print("The number of rows in the TSV files is different.")