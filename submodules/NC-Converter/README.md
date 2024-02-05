# How to convert .nc files to .csv

This project aims at converting **.nc** files -used for storing and exchanging data from satellite observations, climate models, and simulation outputs- to **.csv** files. To do so, run `script.py` as described below:

Create a virtual environment and activate it with the following commands:

        python3 -m venv virtual_env
        source virtual_env/bin/activate   
        
‎Install -if needed- the needed packages with:
    
         pip3 install netCDF4 xarray pandas numpy

‎Run the script that will convert the **.csv** files in data/csv_files from **.nc** files in data/nc_files, you can find writing conventions with the following line:
    
        python3 script.py --help

<p align="center">
	<a href="https://github.com/VicVEVO/NC-Converter"><img src="https://github.com/VicVEVO/NC-Converter/blob/093bdae822a1aa869fb05acab18adc19a1eccd68/resources/help_print.png" width="1200"></a>
</p>

Then, you can access and select the fields you need as follows:

<p align="center">
	<a href="https://github.com/VicVEVO/NC-Converter"><img src="https://github.com/VicVEVO/NC-Converter/blob/e51cb12e796db9b59f0b634f434696504db78c32/resources/fields_print.png" width="1200"></a>
</p>

‎Deactivate -if needed- the virtual environment:

        deactivate
