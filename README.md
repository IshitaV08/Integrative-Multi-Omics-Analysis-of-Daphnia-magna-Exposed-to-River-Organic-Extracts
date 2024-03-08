#To merge metabolome dataset
import pandas as pd

# Read the datasets
df1 = pd.read_csv('sample_sheet.csv')
df2 = pd.read_csv('Metabolome Trasposed.csv')

# Merge the datasets based on the first column
merged_df = pd.merge(df1, df2, on='SampleID', how='inner')
# Define the directory path where you want to save the CSV file
directory_path = '/rds/homes/i/ixv315/Module 5/'

# Define the file name
file_name = 'metabolomic_merged_df.csv'

# Combine directory path and file name
file_path = directory_path + file_name

# Output the transposed DataFrame to a CSV file in the specified directory
merged_df.to_csv (file_path, index=False)
