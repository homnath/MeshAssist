# Makefile for MeshAssist
# DEVELOPER
#   Hom Nath Gharti
#   hngharti_AT_gmail_DOT_com 
#
# HISTORY
#   Feb 07,2013, HNG
#
# source and binary directory

default: all

all: createdir gocad2vtu vtk2jou xyz2jou

clean:
	(cd src; make clean)

createdir:
	(mkdir -p bin; mkdir -p output)

gocad2vtu:
	(cd src; make $@)

vtk2jou:
	(cd src; make $@)

xyz2jou:
	(cd src; make $@)
