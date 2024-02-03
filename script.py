import netCDF4                        
import xarray as xr                                           
import os
import pandas as pd
import numpy as np

parent_dir = os.getcwd()
file_path = os.path.join(parent_dir, 'data/nc_files/')
output_path = os.path.join(parent_dir, 'data/csv_files/')
os.makedirs(output_path, exist_ok=True)

files = [item for item in os.listdir(file_path) if item.endswith('.nc')]
datasets = [xr.open_dataset(os.path.join(file_path, file)) for file in files]

print("Files found : ", files)
print("\nProcessing to csv...\n")

for i in range(len(datasets)):
    output_filename = f"{output_path}{files[i][:-3]}.csv"

    # Créer un DataFrame vide pour stocker toutes les variables
    all_variables_df = pd.DataFrame()

    # Liste pour stocker les DataFrames temporaires de chaque variable
    variable_dfs = []

    # Trouver la longueur maximale parmi toutes les dimensions
    max_length = max(variable_data.size for variable_data in datasets[i].variables.values())

    for variable_name, variable_data in datasets[i].variables.items():
        # Ignorer les variables non numériques (coordonnées, etc.)
        if variable_data.dtype.kind not in 'iufc':
            continue

        # Ajuster la longueur des valeurs de la variable
        flattened_data = variable_data.values.ravel()
        flattened_data = np.pad(flattened_data, (0, max_length - len(flattened_data)))

        # Ajouter la variable au DataFrame temporaire
        variable_df = pd.DataFrame({variable_name: flattened_data})
        variable_dfs.append(variable_df)

    # Concaténer tous les DataFrames temporaires en un seul DataFrame
    all_variables_df = pd.concat(variable_dfs, axis=1)

    # Écrire le DataFrame dans un fichier CSV
    all_variables_df.to_csv(output_filename, index=False)
    print(f"{files[i][:-3]}.csv is rendered.")

