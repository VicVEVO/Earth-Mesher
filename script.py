import netCDF4
import xarray as xr
import os
import pandas as pd
import numpy as np

parent_dir = os.getcwd()
file_path = os.path.join(parent_dir, 'data/nc_files/')
output_path = os.path.join(parent_dir, 'data/nc_files_netcdf/')
output_path2 = os.path.join(parent_dir, 'data/csv_files/')
os.makedirs(output_path, exist_ok=True)

files = [item for item in os.listdir(file_path) if item.endswith('.nc')]

print("Files found : ", files)
print("\nProcessing to NetCDF...\n")


for file in files:
    output_filename = os.path.join(output_path, file)
    ds = xr.open_dataset(os.path.join(file_path, file))
    ds.to_netcdf(output_filename)
    print(f"{file} is converted to NetCDF format.")

print("\nProcessing to CSV...\n")

files = [item for item in os.listdir(output_path) if item.endswith('.nc')]
datasets = [xr.open_dataset(os.path.join(output_path, file)) for file in files]

for i in range(len(datasets)):
    output_filename = f"{output_path2}{files[i]}"

    # Créer un DataFrame vide pour stocker toutes les variables
    all_variables_df = pd.DataFrame()

    # Liste pour stocker les DataFrames temporaires de chaque variable
    variable_dfs = []

    # Trouver la longueur maximale parmi toutes les dimensions
    L = []
    for variable_name, variable_data in datasets[i].variables.items():
        # Ignorer les variables non numériques (coordonnées, etc.)
        if variable_data.dtype.kind not in 'iufc':
            continue

        # Trouver la longueur de chaque dimension de la variable
        lengths = [len(dim) for dim in variable_data.dims]
        L.append(max(lengths))

    max_length = max(L)

    for variable_name, variable_data in datasets[i].variables.items():

        # Ignorer les variables non numériques (coordonnées, etc.)
        if variable_data.dtype.kind not in 'iufc':
            continue

        # Ajuster la longueur des valeurs de la variable
        flattened_data = variable_data.values.ravel()
        flattened_data = np.pad(flattened_data, ((0, max_length - len(flattened_data)) if max_length > len(flattened_data) else (0,0)), mode='constant')

        # Ajouter la variable au DataFrame temporaire
        variable_df = pd.DataFrame({variable_name: flattened_data})
        variable_dfs.append(variable_df)

    # Concaténer tous les DataFrames temporaires en un seul DataFrame
    all_variables_df = pd.concat(variable_dfs, axis=1)

    # Écrire le DataFrame dans un fichier CSV
    all_variables_df.to_csv(output_filename.replace('.nc', '.csv'), index=False)
    print(f"{files[i]} is rendered to CSV.")

