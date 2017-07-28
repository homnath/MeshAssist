# Getting started                                                        
                
\section pack_struct Package structure
                                                                
 The MeshAssist pacakge can be obtained using Git. Use the following command in the terminal:

 git clone \-\-recursive https://github.com/homnath/MeshAssist.git

 The package has the following structure:

## MeshAssist
### doc/: Documentation files including this one. 
### src/ : Contains all source files.
### input/: Contains example input files.
### LICENSE: License.
### Makefile: GNU make file.
 
\section prereq Prerequisites

Package requires Make utility, latest C and Fortran compilers. For matlab files, Matlab is necessary.

\section config Configuration

Open src/Makefile and modify the C and Fortran compilers if necessary.

\section compile Compile

Type the following command in the terminal

make all

Matlab files can be opened in and run from Matlab.

\section run Run

[command] [input_file]

Example:

./bin/xyz2jou ./input/xyz2jou_example.utm

See Chapter "File Documentation" for all available commands. 
