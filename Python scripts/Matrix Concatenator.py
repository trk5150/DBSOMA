# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:39:08 2024

@author: tik105

Takes all files from directory, assumes tab delimited data matrix, concatentates them into 1 matrix, saves as tab delimited file
"""

import os
import pandas as pd

# Directory containing the files
directory = "C://Users//tik105//Desktop//mRNA//Old v young brain//Matricies"

# List all files in the directory
files = os.listdir(directory)

# Initialize an empty list to store DataFrames
dataframes = []

# Loop through each file
for file in files:
    # Check if the file is a CSV file (you can modify this condition for other file types)
    if file.endswith('.txt'):
        # Construct the full file path
        file_path = os.path.join(directory, file)
        
        # Read the file into a DataFrame
        df = pd.read_csv(file_path, sep='\t')
        
        # Append the DataFrame to the list
        dataframes.append(df)

# Concatenate all DataFrames into a single DataFrame
combined_df = pd.concat(dataframes, ignore_index=True, axis=1)


#Might be necessary in some cases, if files don't all have the same list of genes
#combined_df.fillna(0, inplace=True)


# Now you have a single DataFrame containing data from all files in the directory


# File name for the combined file
output_file_name = 'combined_matricies.tsv'

# Full path for the output file
output_file_path = os.path.join(directory, output_file_name)

# Save the DataFrame as tab-delimited text file
combined_df.to_csv(output_file_path, sep='\t', index=False)

print(combined_df.iloc[:5, :5])
