# Earth Mesher

Using satellite data, create a comprehensive 3D visualization of Earth's measured characteristics (such as altitude, temperature,..). Then, feel free to plot -for instance- Earth's geodesy with the trajectory of the satellite that captured the data.
The data I took for thee example are available there:

## Example by plotting altitude (with points)
<p align="center">
	<a href="https://github.com/VicVEVO/Earth-Mesher/blob/029848795746c3770675454244a5780ebcb422fb/resources/"><img src="https://github.com/VicVEVO/Earth-Mesher/blob/029848795746c3770675454244a5780ebcb422fb/resources/example-1.gif" width="900"></a>
</p>

## Example by plotting altitude (with lines)
<p align="center">
	<a href="https://github.com/VicVEVO/Earth-Mesher/blob/17d05df77ab0ff932cb8bccba06b40b95c17a717/resources/"><img src="https://github.com/VicVEVO/Earth-Mesher/blob/17d05df77ab0ff932cb8bccba06b40b95c17a717/resources/example-2.gif" width="900"></a>
</p>

## How to run the project

**Compile** this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `render` binary. Then, **run** it for a render with points (-P) or a mesh (-S):

    ./render -P (or -L, -S) <FILEPATH/file.csv>

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
