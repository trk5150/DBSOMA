
"""
Created on Mon Jan  8 10:33:45 2024

@author: tik105
"""

import pandas as pd
import os
import csv

printFile = True
directory_path = "C://Users//tik105//Desktop//BigScanTest//Small Molecules Only output//5050 full sc-b//Images"

# Get all file names in the directory
file_names = os.listdir(directory_path)
trimmed_file_names = []
for element in file_names:
    splitter = element.split(" ")[0]
    if len(splitter) > 5:
        trimmed_file_names.append(splitter)

# Create a DataFrame with the file names
df = pd.DataFrame(file_names, columns=['File Names'])

file_info_list = []
header = ["Study","Cell line", "treatment time", "unknown", "drug", "concentration", "up/down"]
file_info_list.append(header)
for file_name in file_names:
    # Split the file name using underscores
    if '_' in file_name:
        parts = file_name.split('_')
        if len(parts) > 4 and len(parts) < 7:
            file_info_list.append(parts)
            parts[len(parts)-1] = parts[len(parts)-1].split(".txt")[0]
            appender = parts[len(parts)-1].split(" ")[-1]
            parts[len(parts)-1] = parts[len(parts)-1].split(" ")[0]
            parts.append(appender)
            if len(parts) == 6:
#                parts[4] = parts[4]+"-gene"
                parts.append(parts[5])    
                parts[5] = "NA"
#        else:
#            print(file_name)

if printFile:
    file_path = "C://Users//tik105//Desktop//BigScanTest//Small Molecules Only output//5050 full sc-b//hits_Analysis_parsed.txt"
    print("printing summary file to " + file_path)
    # Open the file in write mode with a tab delimiter
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
    
        # Write the list as a row in the CSV file
        writer.writerows(file_info_list)
