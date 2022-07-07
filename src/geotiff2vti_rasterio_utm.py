# Convert a geotiff (.tiff/.tif) file into a vti (.vti) file of data points in 
# cartesian coordinate form. Can be opened in ParaView.
#
# Accesses the metadata stored in the geotiff to define the size of the data, 
# the UTM coordinates, and spacing of points.
#
# Rasterio module is easiest to access through Anaconda, or it can be installed 
# through pip
#
# utm module can can be accessed through anaconda or it can be installed with pip
#
# Call the code in Linux terminal using: python geotiff2vti_rasterio_utm.py
#==============================================================================
# Input:
# Only input is the geotiff file as 'inFile'. Be careful with the input file 
# extension; geotiffs can end with '.tif' or '.tiff' which will affect reading.
# Provide file to convert:
inFile = 'bathymetry_sumatra.tiff'
#
# Output:
# User must provide a name for the vti file, as well as the variable of interest
# Provide name for output file
outFile = 'bathymetry_sumatraTest.vti'
# Output variable name
outVar = 'Elevation'
#==============================================================================
import rasterio as rio
import math as mt
import numpy as np
import utm

print('File to convert: %s' % inFile)
print('File to create: %s' % outFile)
print('Variable: %s' % outVar)

# Set constants
rEarth = 6371 #km
degToRad = mt.pi/180

# Open the file as a geotiff and as an array
img = rio.open(inFile)
data = img.read()
d = np.array(np.fliplr(data),dtype='float32')
dataFull = img.read()

# Get origin data
ox = img.bounds.left
oy = img.bounds.bottom
oz = data[0,0]

# Convert the long & lat to UTM
utm = utm.from_latlon(oy,ox)
oxUTM = utm[0]
oyUTM = utm[1]
ozUTM = oz[0]

print('Origin: [UTM Easting: %f,  UTM Northing: %f,  Elev: %f]' % (oxUTM,oyUTM,ozUTM))

# Get grid dimensions in degrees
widthLong = abs(img.bounds.right - img.bounds.left)
heigthLat = abs(img.bounds.top - img.bounds.bottom)

# Convert Long & Lat to distance on the surface
widthKm = widthLong*rEarth*degToRad
heightKm = heigthLat*rEarth*degToRad

print('Data width: %f km' % widthKm)
print('Data height: %f km' % heightKm)
print('Min elevation: %f km' % np.min(data))
print('Max elevation: %f km' % np.max(data))

# Get spacing of the points 
dx = widthKm/(dataFull.shape[2]-1)
dy = heightKm/(dataFull.shape[1]-1)
dz = 1 # May need to modify to run 3D maps

# Get the number of points in the data
nx = dataFull.shape[2]
ny = dataFull.shape[1]
nz = dataFull.shape[0]

# Get the number of bytes of the actual vti data (excludes metadata)
offset=0
ncomp=1
byte_int=4
byte_float=4
nbyte=nx*ny*nz*ncomp*byte_float
nbyte=nbyte+byte_int
bnbyte = np.int32(nbyte)

# Create a new file and start writing to it
# 'wb' opens the file in binary and then use '.encode(ascii)' to write ascii
nfile = open(outFile, 'wb')

# Write ascii
nfile.write('<?xml version=\"1.0\"?>\n'.encode('ascii'))
nfile.write('<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n'.encode('ascii'))
nfile.write('<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%.6f %.6f %.6f\" Spacing=\"%.6f %.6f %.6f\">\n'.encode('ascii') % (nx-1,ny-1,nz-1,oxUTM,oyUTM,ozUTM,dx,dy,dz))
nfile.write('<Piece Extent=\"0 %d 0 %d 0 %d\">\n'.encode('ascii') % (nx-1,ny-1,nz-1))
nfile.write(f'<PointData Scalars=\"{outVar}\">\n'.encode('ascii'))
nfile.write(f'<DataArray type="Float32" Name=\"{outVar}\" format=\"appended\" offset=\"{offset}\" />\n'.encode('ascii'))
nfile.write('</PointData>\n'.encode('ascii'))
nfile.write('<CellData>\n'.encode('ascii'))
nfile.write('</CellData>\n'.encode('ascii'))
nfile.write('</Piece>\n'.encode('ascii'))
nfile.write('</ImageData>\n'.encode('ascii'))
nfile.write('<AppendedData encoding=\"raw\">\n'.encode('ascii'))
nfile.write('_'.encode('ascii'))

# Write binary
nfile.write(bnbyte)
nfile.write(d)

# Write ascii
nfile.write('\n'.encode('ascii'))
nfile.write('</AppendedData>\n'.encode('ascii'))
nfile.write('</VTKFile>\n'.encode('ascii'))
nfile.close()
#==============================================================================