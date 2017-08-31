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

The package requires Make utility, latest C and Fortran compilers. For matlab files, Matlab is necessary.

\section config Configuration

Open src/Makefile and modify the C and Fortran compilers if necessary.

\section compile Compile

Type the following command in the terminal

make all

Matlab files can be opened in and run from Matlab.

\section run Run

\em command \em input_file [\em Options]

Example:

./bin/xyz2jou ./input/xyz2jou_example.utm

See Chapter "File Documentation" for all available commands. 

\section bug Bug Report

hgharti_AT_princeton_DOT_edu
