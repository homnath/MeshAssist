# Makefile for MeshAssist
# REVISION
#   HNG, Feb 07,2013
# source and binary directory
SRC = .
BIN = ../bin

# C compiler
CC = gcc
CFLAGS = -fbounds-check -O5 -fno-trapping-math -Wall -Wunused

# Fortran compiler
FC = gfortran
FFLAGS = -std=f2003 -fbounds-check -fimplicit-none -O5 -pedantic -pedantic-errors -Wall -Wline-truncation -fno-trapping-math -fno-backslash -Wunused
#FFLAGS = -std=f2008 -O5
#FFLAGS = std03 -O3
CC_COMPILE = $(CC) $(CFLAGS)
FC_COMPILE = $(FC) $(FFLAGS)

# main executable for 1D FEM
gocad2vtu = $(BIN)/gocad2vtu
xyz2jou = $(BIN)/xyz2jou

default: all

clean:
	rm -f *.mod $(BIN)/*.o core

all: gocad2vtu xyz2jou

gocad2vtu:
	$(CC_COMPILE) -o $(gocad2vtu) gocad2vtu.c

xyz2jou:
	$(FC_COMPILE) -o $(xyz2jou) xyz2jou.f90
