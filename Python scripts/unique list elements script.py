# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import csv

def load_csv(file_path):
    with open(file_path, 'r') as file:
        names = file.read()
    return names

file1_path = 'C:/Users/tik105/Desktop/Coding/singlecellworkspace/Single Cell Correlation Matrix SOM/Expanding parameters test/best mat 2 shells, 5 count.txt'
file2_path = 'C:/Users/tik105/Desktop/Coding/singlecellworkspace/Single Cell Correlation Matrix SOM/millman paper sets/SC-islets/eec 2 5 pruned on full map.txt'


names_from_file1 = load_csv(file1_path)
names_from_file2 = load_csv(file2_path)

list1 = names_from_file1.split(',')
list2 = names_from_file2.split(',')

#print("Names from file1.csv:", names_from_file1)
#print("Names from file2.csv:", names_from_file2)

unique_names_in_file1 = set(list1) - set(list2)

print("Unique names in file1.txt:", unique_names_in_file1)
