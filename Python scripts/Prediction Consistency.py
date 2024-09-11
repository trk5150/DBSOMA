# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:22:53 2024

@author: tik105
"""

import pandas as pd

# Specify the file path of the tab-delimited text file
file_path = "C://Users//tik105//Desktop//BigScanTest//Maturation//Maturation Drug Out//5050 full sc-b//SOM_Analysis.txt"

# Load the tab-delimited file into a DataFrame
df = pd.read_csv(file_path, delimiter='\t')

# Define thresholds for 'Structure' and 'quality'
structure_threshold = .04  # Replace with your desired threshold value
quality_threshold = .03    # Replace with your desired threshold value
beta_mat_threshold = 0.1

# Filter out rows where 'Structure' or 'quality' are below the specified thresholds
df_filtered = df[(df['Structure'] >= structure_threshold) & (df['quality'] >= quality_threshold) & (df['1_Beta_mat'] >= beta_mat_threshold)]

# Calculate the 'product' column by multiplying the specified columns
df_filtered['product'] = df_filtered['Structure'] * df_filtered['quality'] * df_filtered['1_Beta_mat']

# Select the top 501 rows based on the 'product' column
top_500_df = df_filtered.nlargest(501, 'product')

# Remove the first row (the one with the highest 'product' value) and reset the index
top_500_df = top_500_df.iloc[1:].reset_index(drop=True)

# Create a list from the entries in the first column of the modified DataFrame
top_list = top_500_df.iloc[:, 0].tolist()

