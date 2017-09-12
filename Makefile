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

all: createdir \
     exodus2semgeotech \
     exodusold2semgeotech \
     exodus2specfem2d \
     exodus2specfem3d \
     exodusold2specfem3d \
     gid2semgeotech \
     gocad2vtu \
     vtk1d2jou \
     vtk2d2jou \
     xyz2jou

clean:
	(cd src; make clean)

createdir:
	(mkdir -p bin; mkdir -p output)

exodusall2semgeotech:
	(cd src; make $@)

exodus2semgeotech:
	(cd src; make $@)

exodusold2semgeotech:
	(cd src; make $@)

exodus2specfem2d:
	(cd src; make $@)

exodus2specfem3d:
	(cd src; make $@)

exodusold2specfem3d:
	(cd src; make $@)

gid2semgeotech:
	(cd src; make $@)

gocad2vtu:
	(cd src; make $@)

vtk1d2jou:
	(cd src; make $@)

vtk2d2jou:
	(cd src; make $@)

xyz2jou:
	(cd src; make $@)
