# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 11:39:08 2024

@author: tik105
"""

import random

# Open the file in read mode and read all lines into a list
with open('C://Users//tik105//Desktop//LINCs sets//Controls//All L1000 genes (inferred + landmark).txt', 'r') as file:
    lines = file.readlines()

# Process the lines as needed
for line in lines:
    print(line.strip())



# Assuming you already have the 'lines' list from the previous code snippet

n = 100000

# Loop n times
for i in range(n):
    # Generate a random subset of 250 elements
    random_subset = random.sample(lines, 100)
    
    # Specify the filename for the new file
    output_directory = 'C://Users//tik105//Desktop//LINCs sets//Controls//100 genes random outputs'
    current = i
    output_file = output_directory + '//random_subset' + str(current) +'.txt'
    
    # Open the file in write mode
    with open(output_file, 'w') as file:
        # Write each line of the random subset to the file
        for line in random_subset:
            file.write(line)
    if(i%1000 == 0):
        print(f"Random subset saved to '{output_file}'.")

