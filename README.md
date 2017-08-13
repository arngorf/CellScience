# CellScience
This software was written to accompany my thesis in examining a FIB-SEM dataset. Preprocessing, semi-automatic segmentation, annotation and storing of organelles, distance function generation, 3D visualization and data outputting is implemented. The software is very rough at the edges in the current form and needs both restructuring, and cleaning up.

# How to use
Dependencies:
1) cmake
2) Qt5
3) VTK
4) Eigen 3 (Removal in progress)

The program has been used on Arch Linux, and has not been tested or used on any other system.

The CMakeLists.txt file is written a bit sketchy with my local paths visible since a general lookup for the required files were not completed. Modifying these manually is thus required.

Assuming the libraries are found by cmake in the system, the program can be run as follows
1) Open a terminal
2) Enter the root directory of the project
3) #mkdir build
4) #cd build
5) #cmake ..
6) #make
7) #./mt

The program will here scream at you for not supplying a 3D tiff image file. You will have to manually change the directory in NumericToolControl class.

Feel free to send questions my way to arngorf@gmail.com
