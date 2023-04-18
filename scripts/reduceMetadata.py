from datetime import datetime

#script to reduce size of metadata file - only include lines with dates prior to 2021-08-31

# Input and output file paths
input_file = 'metadata.tsv'
output_file = 'metadata-reduced.tsv'

# Predefined date for comparison
compare_date = datetime(2021, 8, 31)

# Open input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    header = infile.readline().strip().split('\t')  # Read and parse the header row
    date_index = header.index('date')  # Find the index of the "date" column
    
    # Write the header row to the output file
    outfile.write('\t'.join(header) + '\n')
    
    for line in infile:
        # Split the line by tab character to extract columns
        columns = line.strip().split('\t')
        
        # Extract the "date" column and parse it as a datetime object
        date_column = columns[date_index]
        if len(date_column) == 10 and "X" not in date_column:
            date_object = datetime.strptime(date_column, '%Y-%m-%d')
        else:
            date_object = datetime.strptime("2023-01-01", '%Y-%m-%d')
        
        # Compare the date object with the predefined date
        if date_object < compare_date:
            # Write the line to the output file
            outfile.write(line)
