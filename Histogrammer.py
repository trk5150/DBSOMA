# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:24:43 2024

@author: tik105
"""

import matplotlib.pyplot as plt

import pandas as pd

# Load data from TSV file
data = pd.read_csv("C://Users//tik105//Desktop//BigScanTest//Small Molecules Only output//5050 full sc-b//SOM_Analysis.txt", sep='\t')




column_data = data["quality"]

# Create histogram
plt.hist(column_data, bins=10, color='skyblue', edgecolor='black')

# Add labels and title
plt.xlabel('Quality Value')
plt.ylabel('Count')
plt.title('Quality')


# Set y-axis to logarithmic scale
plt.yscale('log')
plt.xlim(min(column_data), max(column_data))

# Show plot
plt.show()

# Define bin edges
bins = [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, .275, 0.3]  # Adjust as per your requirement

# Create bins and count occurrences
bin_counts = pd.cut(data['quality'], bins=bins).value_counts().sort_index()

# Get bin centers
bin_centers = [(bins[i] + bins[i+1]) / 2 for i in range(len(bins)-1)]

# Display bin centers and bin counts
print("Bin Centers:", bin_centers)
print("Bin Counts:", bin_counts.tolist())
