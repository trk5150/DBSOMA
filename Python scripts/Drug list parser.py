# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:48:00 2024

@author: tik105

Code from chatGPT
"""
import csv

# Open the file containing the list of line-broken strings
file_path = "C://Users//tik105//Desktop//BigScanTest//Maturation//Maturation OE Out//5050 full sc-b//Hit files.txt"  # Replace with the path to your file

# Read the lines from the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Parse each line by underscore
parsed_strings = []
for line in lines:
    # Split the line by underscore
    parts = line.strip().split('_')  # strip() removes leading/trailing whitespaces
    
    # Trim the ".txt" extension from the last element
    last_element = parts[-1]
    last_element_without_extension = last_element.replace('.txt', '')

    # Split the last element further at the space
    last_element_parts = last_element_without_extension.split()

    # Append the parsed parts to the list
    parsed_strings.append(parts[:-1] + last_element_parts)


# Print the parsed strings
for parsed_string in parsed_strings:
    print(parsed_string)

# Write the parsed strings to a tab-delimited file
output_file_path = "C://Users//tik105//Desktop//BigScanTest//Maturation//Maturation OE Out//5050 full sc-b//Parsed.txt"
with open(output_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerows(parsed_strings)