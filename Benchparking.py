# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 10:27:02 2024

@author: tik105
"""

import pandas as pd

# Replace 'your_file.tsv' with the path to your tab-delimited file
SOMSCAN_path = "C://Users//tik105//Desktop//BigScanTest//SOM, CLue and enrichr//Maturation SOM compared to primary islet list//SOMSCAN.txt"
enrichr_path = "C://Users//tik105//Desktop//BigScanTest//SOM, CLue and enrichr//Maturation SOM compared to primary islet list//Enrichr.txt"
clue_io_path = "C://Users//tik105//Desktop//BigScanTest//SOM, CLue and enrichr//Maturation SOM compared to primary islet list//Clue_io.txt"
target_list_path = "C://Users//tik105//Desktop//BigScanTest//SOM, CLue and enrichr//Maturation SOM compared to primary islet list//Maturation Chemicals list.txt"

# Read the tab-delimited file into a DataFrame
SOMSCAN = pd.read_csv(SOMSCAN_path, delimiter='\t')
enrichr = pd.read_csv(enrichr_path, delimiter='\t')
clue_io = pd.read_csv(clue_io_path, delimiter='\t')
target_list = pd.read_csv(target_list_path, delimiter='\t', header=None)




#parse the pert names from SOMSCAN
# Split the column and create new columns
split_columnsE = enrichr['Term'].str.split(' ', expand=True)
split_columnsE.columns = ['drug', 'Up_down', 'blank', 'blank']

# Assign the split columns to new columns in the DataFrame
enrichr = pd.concat([enrichr, split_columnsE], axis=1)
enrichr['rank'] = enrichr['Adjusted P-value'].rank(method='first', ascending=True)


#rename the clue_io column to be consistent
clue_io.rename(columns={'pert_iname': 'drug'}, inplace=True)
clue_io['rank'] = clue_io['norm_cs'].rank(method='first', ascending=False)

#rename the SOM consistency column to be consistent
SOMSCAN.rename(columns={'Compound': 'drug'}, inplace=True)
SOMSCAN['rank'] = SOMSCAN['Product'].rank(method='first', ascending=False)


# Create a new DataFrame to collect ranks
result = pd.DataFrame(columns=['drug', 'SOMSCAN_rank', 'Clue_io_rank', 'enrichr_rank'])

# Function to get the first rank of an identifier in a dataframe
def get_first_rank(SOMSCAN, identifier_col, rank_col, identifier):
    SOMSCAN[identifier_col] = SOMSCAN[identifier_col].str.lower().str.strip()
    identifier = identifier.lower().strip()
    row = SOMSCAN[SOMSCAN[identifier_col] == identifier]
    if not row.empty:
        return row.iloc[0][rank_col]
    return -1
"""
# Iterate through each identifier in df1
for index, row in target_list.iterrows():
    identifier = row[0]
    rank1 = get_first_rank(SOMSCAN, 'drug', 'rank', identifier) 
    rank2 = get_first_rank(clue_io, 'drug', 'rank', identifier)
    rank3 = get_first_rank(enrichr, 'drug', 'rank', identifier)
    
    # Append the collected ranks to the result dataframe
    result = result.append({'drug': identifier, 'SOMSCAN_rank': rank1, 'Clue_io_rank': rank2, 'enrichr_rank': rank3}, ignore_index=True)
"""


# Iterate through each identifier in df1
for index, row in SOMSCAN.iterrows():
    identifier = row['drug']
    rank1 = row['rank']
    
    rank2 = get_first_rank(clue_io, 'drug', 'rank', identifier)
    rank3 = get_first_rank(enrichr, 'drug', 'rank', identifier)
    
    # Append the collected ranks to the result dataframe
    result = result.append({'drug': identifier, 'SOMSCAN_rank': rank1, 'Clue_io_rank': rank2, 'enrichr_rank': rank3}, ignore_index=True)


#output
# Define the file path where you want to save the TSV file
file_path = "C://Users//tik105//Desktop//BigScanTest//SOM, CLue and enrichr//Maturation SOM compared to primary islet list//benching.txt"

# Export DataFrame to TSV
result.to_csv(file_path, sep='\t', index=False)