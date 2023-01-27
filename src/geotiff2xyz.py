# Convert a geotiff (.tiff/.tif) into a vertical list of xyz coordinates
# Can be used with the xyz2jpu.f90 to turn .geotiff into .jou for CUBIT use
#
# Rasterio module is easiest to access through Anaconda, or it can be installed 
# through pip
#
# utm module can can be accessed through anaconda or it can be installed with pip
#------------------------------------------------------------------------------
# Input:
# Only input is the geotiff file as 'inFile'. Be careful with the input file 
# extension; geotiffs can end with '.tif' of '.tiff' which will effect reading.
# Provide file to convert:
inFile = 'bathymetry_sumatra.tiff'
#
# Output:
# User must provide a name for the vti file, as well as the variable of interest
# Provide name for output file
outFile = 'bathymetry_sumatra_xyz.utm'
#==============================================================================
import rasterio as rio
import math as mt
import numpy as np
import utm

print('File to convert: %s' % inFile)
print('File to create: %s' % outFile)

# Set constants
rEarth = 6371 
degToRad = mt.pi/180
m2km = 0.001 # *May need to adujst this parameter* #

# Open the file as a geotiff and as an array
img = rio.open(inFile)
data = img.read(1)
d = np.array(np.fliplr(data))
dataFull = img.read()

# Get origin data
ox = img.bounds.left
oy = img.bounds.bottom

# Convert the long & lat into UTM
utm = utm.from_latlon(oy,ox)
oxUTM = utm[0]
oyUTM = utm[1]

print('Origin: [UTM Easting: %f,  UTM Northing: %f]' % (oxUTM,oyUTM))

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
maxElev = np.max(data)
minElev = np.min(data)

# Get spacing of the points 
dx = widthKm/(dataFull.shape[2]-1)
dy = heightKm/(dataFull.shape[1]-1)

# Get the number of points in the data
nx = dataFull.shape[2]
ny = dataFull.shape[1]
nz = dataFull.shape[0]

# Create a new file and start writing to it 
nfile = open(outFile, 'w')

for j in range(ny):
    x = oxUTM
    y = oyUTM + (dy*j)
    for i in range(nx):
        nfile.write('%f      %f      %f' % (x+(dx*i),y,d[j,i]*m2km))
        if i!= nx-1:
            nfile.write('\n')
    if j != ny-1:
        nfile.write('\n')
nfile.close()
#==============================================================================