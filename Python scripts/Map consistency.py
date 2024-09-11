# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:54:15 2024

@author: tik105
"""

import pandas as pd
from collections import Counter
# Load the Excel file
file_path = "C://Users//tik105//Desktop//BigScanTest//Maturation//Map consistency.xlsx"
df = pd.read_excel(file_path)

# Assuming the columns are named 'List1', 'List2', 'List3', 'List4'
drug = df['Drugs'].dropna().tolist()
list1 = df['Map 1'].dropna().tolist()
list2 = df['Map 2'].dropna().tolist()
list3 = df['Map 3&4'].dropna().tolist()
list4 = df['Map 5'].dropna().tolist()

# Combine all lists into one
all_filenames = list1 + list2 + list3 + list4

# Count occurrences of each filename
filename_counts = Counter(all_filenames)

# Separate filenames based on their counts
in_all_4_lists = [file for file, count in filename_counts.items() if count == 4]
in_3_lists = [file for file, count in filename_counts.items() if count == 3]
in_2_lists = [file for file, count in filename_counts.items() if count == 2]
in_1_list = [file for file, count in filename_counts.items() if count == 1]

# Count the occurrences of each filename in all_filenames
filename_counts = Counter(all_filenames)

# For each filename in list1, get the count from the counter
drug_counts = {filename: filename_counts[filename] for filename in drug}

# Convert the dictionary to a DataFrame
drug_counts_df = pd.DataFrame(drug_counts.items(), columns=['Filename', 'Count'])

# Specify the output file path
output_file_path = "C://Users//tik105//Desktop//BigScanTest//Maturation//list1_counts.csv"

# Export the DataFrame to a CSV file
drug_counts_df.to_csv(output_file_path, index=False)