# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 15:07:19 2023

@author: tik105

reads a given input file, assumes each line is a new data set, names the file based on the first 2 elements, then saves a line-broken text file with the remaining elements
"""
import os 
def process_line(line, output_dir):
    data = line.split('\t')

    if data and len(data) > 30:
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
        output_dir = "C://Users//tik105//Desktop//LINCs sets//GMT Parsed Genes//siRNA//Output"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for filename in os.listdir(input_directory):
            if filename.endswith(".gmt"):
                input_file = os.path.join(input_directory, filename)
                process_input_file(input_file, output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    input_directory = "C://Users//tik105//Desktop//LINCs sets//GMT Parsed Genes//siRNA"  # Change this to your input directory's path
    main(input_directory)    
    
#f'{data[0]}_{data[1]}_{inName}.txt')
        
#output_dir = "C:/Users/tik105/Desktop/Coding/singlecellworkspace/Single Cell Correlation Matrix SOM/search files/crisprKO/outputs/"
#input_directory = "C:/Users/tik105/Desktop/Coding/singlecellworkspace/Single Cell Correlation Matrix SOM/search files/crisprKO" # Change this to your input directory's path
