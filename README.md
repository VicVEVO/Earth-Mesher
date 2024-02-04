# Earth Mesher

This project involves the study of geodesy, i.e. the shape of the Earth, using measurements of its gravity field. Lemme cook

## TODO:
[ ] Prendre en compte les .nc en HDF5 pour les convertir en NetCDF: tout traduire ou déterminer si un .nc est en HDF5 ? (plus dur)
[ ] Ne garder que latitude, longitude, altitude pour les .csv pour qu'ils soient moins lourds: prendre en compte commande utilisateur disant * (on prend tout pour le .csv) ou juste certains champs. 


## How to convert .nc files to .csv

This project uses NASA's TOPEX/POSEIDON geophysical data in **.nc** files, to convert them to **.csv**, use `script.py` as follows:
Create a virtual environment and activate it with the following commands:

        python3 -m venv virtual_env
        source virtual_env/bin/activate   
        
‎Install -if needed- the needed packages with:
    
         pip3 install netCDF4, xarray, os, pandas, numpy

‎Run the script that will convert the **.csv** files in data/csv_files from **.nc** files in data/nc_files:
    
        python3 script.py

‎Deactivate -if needed- the virtual environment:

        deactivate
        
## How to run the project

**Compile** this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `render` binary. Then, **run** it with:

    ./render

## Dependencies

The dependencies are STL, Eigen, libigl and the dependencies
of the `igl::opengl::glfw::Viewer` (OpenGL, glad and GLFW).

The latter requires the glfw module of libigl, which shows up in the CMakeLists.txt 

```cmake
igl_include(glfw)
…
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
```

The CMake build system will automatically download libigl and its dependencies using
[CMake FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html),
thus requiring no setup on your part.

### Use a local copy of libigl
You can use the CMake cache variable `FETCHCONTENT_SOURCE_DIR_LIBIGL` when configuring your CMake project for
the first time to aim it at a local copy of libigl instead.
```
cmake -DFETCHCONTENT_SOURCE_DIR_LIBIGL=<path-to-libigl> ..
```
When changing this value, do not forget to clear your `CMakeCache.txt`, or to update the cache variable
via `cmake-gui` or `ccmake`.
