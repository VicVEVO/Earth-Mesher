import netCDF4                        
import xarray as xr                                           
import os
import sys
import pandas as pd
import numpy as np

### Parameters analysis

args = sys.argv
stop = True

if len(args) == 1:
    print("Please provide the correct data field names when running the script.\n\nFor further information, run:\n >> python3 script.py --help")
else:    
    option = args[1]

    if option == "-h" or option == "--help":
        print("Usage: python3 script.py <option>\nSupported options include:\n\n--data <field1> <field2> ...\n     Specity the data to keep after the conversion.\n     To keep only the temperature and latitude measurements >> python3 script.py temperature latitude.\n     To keep all fields >> python3 script.py *\n\n--fields\n     Show the different fields of the database.\n\n--help\n     Show publicly available flags.\n")
    else:
        parent_dir = os.getcwd()
        file_path = os.path.join(parent_dir, 'data/nc_files/')
        output_path = os.path.join(parent_dir, 'data/csv_files/')
        os.makedirs(output_path, exist_ok=True)
        files = [item for item in os.listdir(file_path) if item.endswith('.nc')]
        datasets = [xr.open_dataset(os.path.join(file_path, file)) for file in files]

        if option == "--fields":
            print(datasets)
        elif option == "--data":
            stop = False

if stop:
    sys.exit()

## Main program
    
print("Files found : ", files)
print("\nProcessing to csv...\n")

for i in range(len(datasets)):
    output_filename = f"{output_path}{files[i][:-3]}.csv"

    # Empty DataFrame creation to stock all the variables
    all_variables_df = pd.DataFrame()

    # List to stock temporaries dataframes of each variable
    variable_dfs = []

    max_length = max(variable_data.size for variable_data in datasets[i].variables.values())
    
    for variable_name, variable_data in datasets[i].variables.items():
        # We pass over not requested or non-numerical data
        if not(variable_name in args) or variable_data.dtype.kind not in 'iufc':
            continue
        
        # Variable's values' length adjustement
        flattened_data = variable_data.values.ravel()
        flattened_data = np.pad(flattened_data, (0, max_length - len(flattened_data)))

        # Addition of the variable the the temporary dataframe
        variable_df = pd.DataFrame({variable_name: flattened_data})
        variable_dfs.append(variable_df)

    # Concatenation of every temporary dataframes
    all_variables_df = pd.concat(variable_dfs, axis=1)

    # Writing the final dataframe in a .csv file
    all_variables_df.to_csv(output_filename, index=False)
    print(f"{files[i][:-3]}.csv is rendered.")

