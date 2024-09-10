# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 15:46:10 2024

@author: tik105
"""
import pandas as pd
file1 = "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//combined_matricies.tsv"
file2 =  "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//GSM3722110_OX3X_10X.txt"

# Read the first column of the first TSV file
column_to_insert = pd.read_csv(file2, sep='\t', usecols=[0], header=None)

# Remove the first value from the column
column_to_insert = column_to_insert.iloc[1:]

# Reindex the column
column_to_insert.reset_index(drop=True, inplace=True)

# Read the second TSV file into a DataFrame
df = pd.read_csv(file1, sep='\t')

# Insert the column as the first column of the DataFrame

df.insert(0, "Genes", column_to_insert)




# Save the modified DataFrame to a new TSV file
df.to_csv("C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//combined_matricies_and-genes.tsv", sep='\t', index=False)
print("saved to C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies//combined_matricies_and-genes.tsv")

