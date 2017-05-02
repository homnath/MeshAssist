# Makefile for MeshAssist
# REVISION
#   HNG, Feb 07,2013
# source and binary directory
SRC = .
BIN = ../bin

# C compiler
CC = gcc
CFLAGS = -lm -fbounds-check -O5 -fno-trapping-math -Wunused

# Fortran compiler
FC = gfortran
FFLAGS = -std=f2003 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wall -Wline-truncation -fno-trapping-math -fno-backslash -Wunused
#FFLAGS = -std=f2008 -O5
#FFLAGS = std03 -O3
CC_COMPILE = $(CC) $(CFLAGS)
FC_COMPILE = $(FC) $(FFLAGS)

# main executable for 1D FEM
exodus2specfem2d = $(BIN)/exodus2specfem2d
exodus2specfem3d = $(BIN)/exodus2specfem3d
exodusold2specfem3d = $(BIN)/exodusold2specfem3d
gocad2vtu = $(BIN)/gocad2vtu
vtk1d2jou = $(BIN)/vtk1d2jou
vtk2d2jou = $(BIN)/vtk2d2jou
xyz2jou = $(BIN)/xyz2jou

default: all

clean:
	rm -f *.mod $(BIN)/*.o core

all: exodus2specfem2d \
     exodus2specfem3d \
     exodusold2specfem3d \
     gocad2vtu \
     vtk1d2jou \
     vtk2d2jou \
     xyz2jou

exodus2specfem2d:
	$(CC_COMPILE) -o $(exodus2specfem2d) exodus2specfem2d.c

exodus2specfem3d:
	$(CC_COMPILE) -o $(exodus2specfem3d) exodus2specfem3d.c

exodusold2specfem3d:
	$(CC_COMPILE) -o $(exodusold2specfem3d) exodusold2specfem3d.c

gocad2vtu:
	$(CC_COMPILE) -o $(gocad2vtu) gocad2vtu.c

vtk1d2jou:
	$(FC_COMPILE) -o $(vtk1d2jou) vtk1d2jou.f90

vtk2d2jou:
	$(CC_COMPILE) -o $(vtk2d2jou) vtk2d2jou.c

xyz2jou:
	$(FC_COMPILE) -o $(xyz2jou) xyz2jou.f90
