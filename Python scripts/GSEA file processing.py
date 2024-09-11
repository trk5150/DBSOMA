import os

def process_line(line, output_dir):
    data = line.strip().split('\t')
    if data:
        output_filename = os.path.join(output_dir, f'{data[0]}_{data[1]}_up.txt')
        with open(output_filename, 'w') as output_file:
            output_file.write('\n'.join(data) + '\n')
        print(f"Created {output_filename}")

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
        output_dir = "C:/Users/tik105/Desktop/Coding/singlecellworkspace/Single Cell Correlation Matrix SOM/search files/msig output Dir"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for filename in os.listdir(input_directory):
            if filename.endswith(".symbols.gmt"):
                input_file = os.path.join(input_directory, filename)
                process_input_file(input_file, output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    input_directory = "C:/Users/tik105/Desktop/Coding/singlecellworkspace/Single Cell Correlation Matrix SOM/search files/msigdb_v2023.1.Hs_GMTs"  # Change this to your input directory's path
    main(input_directory)
