# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 11:36:53 2024

@author: tik105
"""

import shutil
import os

def read_key_strings_from_file(file_path):
    with open(file_path, 'r') as file:
        key_strings = [line.strip() for line in file if len(line.strip()) > 3]
    return key_strings

def copy_files(source_dir, destination_dir, key_strings, key_strings_2):
    # Get a list of all files in the source directory
    files = os.listdir(source_dir)
    
    # Iterate through each file and copy it to the destination directory
    for file in files:
        source_file = os.path.join(source_dir, file)
        # Check if the filename contains any of the key strings
        for key in key_strings:
            if key.lower() in file.lower():
                for key2 in key_strings_2:
                    if(key2 in file):
                        destination_file = os.path.join(destination_dir, file)
                        # Copy the file to the destination directory
                        shutil.copy(source_file, destination_file)
                        # Print a message
                        print(f"Copied {file} to {destination_dir} (Key: {key})")

# Example usage:
source_directory = "C://Users//tik105//Desktop//LINCs sets//GMT parsed only small molecules"
destination_directory = "C://Users//tik105//Desktop//LINCs sets//MCF10 and HepG2"

key_strings_file = "C://Users//tik105//Desktop//SOMSCAN//Figures//Figure 3 Structure overlap DBSCAN//HepG2MCF10.txt"
key_strings_file_2 = "C://Users//tik105//Desktop//SOMSCAN//Figures//Figure 3 Structure overlap DBSCAN//drugs.txt"

# Read key strings from the file
key_strings = read_key_strings_from_file(key_strings_file)
key_strings_2 = read_key_strings_from_file(key_strings_file_2)

copy_files(source_directory, destination_directory, key_strings, key_strings_2)