import os

# Specify the file path to your tab-delimited file
file_path = "C:/Users/tik105/Downloads/differential expression of genes in all datasets(Human).txt"
output_dir = "C:/Users/tik105/Desktop/KnockTF/Human/250"

with open(file_path, 'r') as f:
    lines = f.readlines()
data_dict = {}  # Dictionary to store data for each unique value

limit = 250

for line in lines:
    columns = line.strip().split('\t')
    if len(columns) >= 10 and columns[10].isdigit() and int(columns[10])>0:
        result = '_up' if int(columns[10]) >1.5 else '_down'
        key = columns[1] + result + '_' + str(limit)
        value = columns[2] 

        if key not in data_dict:
            data_dict[key] = []
            count = 0
           # print(key)
        if value not in data_dict[key] and count < limit:
            data_dict[key].append(value)
            count+=1;
#print ("data dict done" + key)



for key, values in data_dict.items():
    
    output_filename = os.path.join(output_dir, f'{key.replace(":", "=")}.txt')
    with open(output_filename, 'w') as output_file:
        output_file.write('\n'.join(values) + '\n')
   # print(f"Created {output_filename}")