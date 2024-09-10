# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 11:49:36 2024

@author: tik105
"""
import csv
import os
import pandas as pd

printFile = True

directory_path = "C://Users//tik105//Desktop//mRNA//single_cell_read_counts//Stage 6//L1000 out"

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
            
##** To print a summary of the identified hits above          
if printFile:
    file_path = "C://Users//tik105//Desktop//mRNA//single_cell_read_counts//Stage 6//L1000 out//hits_Analysis.txt"
    print("printing summary file to " + file_path)
    # Open the file in write mode with a tab delimiter
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
    
        # Write the list as a row in the CSV file
        writer.writerows(file_info_list)


##** to fill a directory with only those files of interest

df = pd.DataFrame(file_info_list)
drug_list = df[4].tolist()


def process_line(line, output_dir):
    data = line.split('\t')
    
    checker = data[0].split('_')
#    check = checker[4]
    check = data[0]
    printAFile = False 
#    for element in drug_list:
#        if element in check: 
    for element in trimmed_file_names:
        if element in check:
            printAFile = True
#            print(element)
#            print(check)
#            print("\n")
    
    if printAFile and data and len(data) > 30:
        output_filename = os.path.join(output_dir, f'{data[0][:45]}.txt')
        with open(output_filename, 'w') as output_file:
            output_file.write('\n'.join(data) + '\n')
     

def process_input_file(input_file, output_dir):
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            process_line(line, output_dir)

    except FileNotFoundError:
        print(f"Input file {input_file} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main(input_directory):
    try:
        output_dir = "C://Users//tik105//Desktop//mRNA//single_cell_read_counts//Stage 6//L1000 out"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for filename in os.listdir(input_directory):
            if filename.endswith(".gmt"):
                input_file = os.path.join(input_directory, filename)
                process_input_file(input_file, output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    input_directory = "C://Users//tik105//Desktop//mRNA//single_cell_read_counts//Stage 6//L1000 out"  # Change this to your input directory's path
    main(input_directory)