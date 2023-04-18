from datetime import datetime

# Input and output file paths
input_file = 'metadata-reduced.tsv'

# set max and min
min_date = datetime(2020, 1, 1)
max_date = datetime(2020, 1, 1)

# Open input and output files
with open(input_file, 'r') as infile:
    header = infile.readline().strip().split('\t')  # Read and parse the header row
    date_index = header.index('date')  # Find the index of the "date" column
    
    for line in infile:
        # Split the line by tab character to extract columns
        columns = line.strip().split('\t')
        
        # Extract the "date" column and parse it as a datetime object
        date_column = columns[date_index]
        if len(date_column) == 10 and "X" not in date_column:
            date_object = datetime.strptime(date_column, '%Y-%m-%d')
            if date_object < min_date:
                min_date = date_object
            if date_object > max_date:
                max_date = date_object

print(f"For {input_file} the maxdate is {max_date} and the mindate is {min_date}")

#For metadata.tsv the maxdate is 2022-02-11 00:00:00 and the mindate is 2010-12-06 00:00:00
#For metadata-reduced.tsv the maxdate is 2021-08-30 00:00:00 and the mindate is 2010-12-06 00:00:00
