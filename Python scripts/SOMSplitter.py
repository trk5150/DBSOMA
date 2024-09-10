# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:19:48 2024

@author: tik105
"""

import pandas as pd
import os

# Replace 'your_file_path' with the path to your tab-delimited file
file_path = "C://Users//tik105//Desktop//BigScanTest//soms//5050 full sc-b.som"

# Read the tab-delimited file into a DataFrame
df = pd.read_csv(file_path, delimiter='\t', skiprows=2, header=None)

# Parse out the x and y coordinates from the first column
df[[ 'x', 'y']] = df[0].str.extract(r'\[(.*?),(.*?)\]').astype(float)

# Drop the original column containing the combined coordinates
df.drop(0, axis=1, inplace=True)





grid_width = 50
grid_height = 50

# Define the number of rows and columns for the split
num_rows = 2
num_cols = 2

def process_string(input_string):
    # Remove spaces, commas, and parentheses from the input string
    cleaned_string = input_string.replace(' ', '').replace(',', '').replace('(', '').replace(')', '')

    # Split the string into a list of elements
    elements = cleaned_string.split(' ')

    return elements

# Define a function to determine which square a coordinate falls into
def get_square_coordinate(x, y):
    # Calculate the width and height of each square
    square_width = grid_width / num_cols
    square_height = grid_height / num_rows

    # Calculate the column index of the square
    col_index = min(int(x / square_width), num_cols - 1)

    # Calculate the row index of the square
    row_index = min(int(y / square_height), num_rows - 1)

    # Return the row and column indices of the square
    return row_index, col_index


# Apply the function to each row of the DataFrame
df['square_coordinate'] = df.apply(lambda row: get_square_coordinate(row['x'], row['y']), axis=1)


# Group the DataFrame by square coordinate
grouped = df.groupby('square_coordinate')

# Create a directory to store the output files
output_dir = "C://Users//tik105//Desktop//Coding//singlecellworkspace//Search files and SOMS//4 density"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Define a function to process the elements in the column
def process_element(element):
    # Remove spaces, commas, and parentheses from the element
    cleaned_element = element.replace(' ', '').replace(',', '\n').replace('(', '').replace(')', '\n')
    return cleaned_element

# Write separate files for each group
for group_name, group_df in grouped:
    # Extract the square coordinate from the group name
    square_coord = '_'.join(str(x) for x in group_name)
    
    # Create a filename for the output file
    output_filename = os.path.join(output_dir, f'square_{square_coord}.txt')
    
    # Apply the processing function to the elements in the column
    processed_column = group_df[1].apply(process_element)
    
    # Write the processed elements to the output file
    with open(output_filename, 'w') as f:
        for item in processed_column:
            f.write(item + '\n')
        
        
        
        
        
        
        