# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 10:43:23 2024

@author: tik105
"""

import csv
import os
import pandas as pd
import json
import warnings

# Suppress all FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

printFile = False

state_file = "C://Users//tik105//Desktop//Coding//singlecellworkspace//Single Cell Correlation Matrix SOM//millman paper sets//primary islet//beta.txt"

# Open the file in read mode
with open(state_file, 'r') as file:
    # Load the content of the file as a JSON-formatted string
    state_as_str = file.read()

    # Parse the JSON string into a Python list
    state_as_list = set(state_as_str.split("\n"))


state_tuned = "C://Users//tik105//Desktop//Coding//singlecellworkspace//Single Cell Correlation Matrix SOM//Expanding parameters test//best mat 2 shells, 5 count - line breaks.txt"

with open(state_tuned, 'r') as file:
    # Load the content of the file as a JSON-formatted string
    tuned_as_str = file.read()
    #tuned_as_str.replace("\t", "")
    # Parse the JSON string into a Python list
    tuned_as_list = set(tuned_as_str.split(" \n"))


#tuned_as_list = state_as_list.union(tuned_as_list)

common_elements = tuned_as_list.intersection(state_as_list)


columns = ['File', 'State Overlap', 'State Overlap %', 'Tuned State Overlap', 'Tuned State Overlap %']
df = pd.DataFrame(columns=columns)

state_over = len(common_elements)
state_over_perc = state_over/len(state_as_list)
tuned_over_perc = state_over/len(tuned_as_list)


#state_data = ["Millman Primary Beta", len(state_as_list), len(state_as_list), 1.0, state_over, state_over_perc]

state_data = {'File':"Millman Primary Beta", 'State Overlap': len(state_as_list), 'State Overlap %' : 1.0, 'Tuned State Overlap':state_over, 'Tuned State Overlap %' : state_over_perc}
df = df.append(pd.Series(state_data), ignore_index=True)
df = df.append({'File':"Primary Beta Tuned", 'State Overlap': state_over, 'State Overlap %' : tuned_over_perc, 'Tuned State Overlap':len(tuned_as_list), 'Tuned State Overlap %' : 1.0}, ignore_index=True)



def process_file(file_path):
    #determine the short name for the file
    namer = file_path.split('_')
    namerr = namer[len(namer)-1].split('.txt')
    namerrr = namer[len(namer)-2]+"_" + namerr[0]
    cell_line = namer[1]
    name = (cell_line + "_" +namer[2] +"_" +namerrr)

    #read in the list of genes
    with open(file_path, 'r') as file:
        # Load the content of the file as a JSON-formatted string
        gene_str = file.read()
        #tuned_as_str.replace("\t", "")
        # Parse the JSON string into a Python list
        genes = set(gene_str.split("\n"))
        common_state = genes.intersection(state_as_list)
        common_tuned = genes.intersection(tuned_as_list)
        state_over = len(common_state)
        tuned_over = len(common_tuned)
        state_over_perc = state_over/len(genes)
        tuned_over_perc = tuned_over/len(genes)
    
    gene_data = {'File':name, 'State Overlap': state_over, 'State Overlap %' : state_over_perc, 'Tuned State Overlap':tuned_over, 'Tuned State Overlap %' : tuned_over_perc}
    return gene_data  


directory_path = 'C://Users//tik105//Desktop//LINCs sets//Hits'
file_list = os.listdir(directory_path)
for file_name in file_list:
    file_path = os.path.join(directory_path, file_name)
    gene_data = process_file(file_path)
    df=df.append(gene_data, ignore_index=True)

