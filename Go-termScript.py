import os

def main(input_file):
    try:
        output_dir = "C:/Users/tik105/Desktop/GO_genes"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(input_file, 'r') as f:
            lines = f.readlines()

        data_dict = {}  # Dictionary to store data for each unique value

        for line in lines:
            columns = line.strip().split('\t')
            if len(columns) >= 2:
                key = columns[1]
                value = columns[0]

                if key not in data_dict:
                    data_dict[key] = []
                if value not in data_dict[key]:
                    data_dict[key].append(value)
        #print ("data dict done" + key)
        
        
        
        for key, values in data_dict.items():
            
            output_filename = os.path.join(output_dir, f'{key.replace(":", "=")}.txt')
            with open(output_filename, 'w') as output_file:
                output_file.write('\n'.join(values) + '\n')
            print(f"Created {output_filename}")

    except FileNotFoundError:
        print("Input file not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    input_file = "C:/Users/tik105/Desktop/goa_human_gosorted.txt"  # Change this to your input file's name
    main(input_file)