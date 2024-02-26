# How to convert .nc files to .csv

This project aims at converting **.nc** files -used for storing and exchanging data from satellite observations, climate models, and simulation outputs- to **.csv** files. To do so, run `script.py` as described below:

Create a virtual environment and activate it with the following commands:

        python3 -m venv virtual_env && source virtual_env/bin/activate   
        
‎Install -if needed- the needed packages with:
    
         pip3 install netCDF4 xarray pandas numpy

‎Run the script that will convert the **.csv** files in data/csv_files from **.nc** files in data/nc_files, you can find writing conventions with the following line:
    
        python3 script.py --help

<p align="center">
	<a href="https://github.com/VicVEVO/NC-Converter"><img src="https://github.com/VicVEVO/Earth-Mesher/blob/19fbf9c9a4bb75ff870b4a6860909bcef8cb7bfd/submodules/NC-Converter/resources/screen2.png" width="1200"></a>
</p>

Then, you can access and select the fields you need as follows:

<p align="center">
	<a href="https://github.com/VicVEVO/NC-Converter"><img src="https://github.com/VicVEVO/Earth-Mesher/blob/19fbf9c9a4bb75ff870b4a6860909bcef8cb7bfd/submodules/NC-Converter/resources/screen1.png" width="1200"></a>
</p>

‎Deactivate -if needed- the virtual environment:

        deactivate
