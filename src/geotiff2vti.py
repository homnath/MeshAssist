# Install geotiff using the command:
# pip install geotiff
from geotiff import GeoTiff as gt
import numpy as np

# Set the GEO TIFF file.
# Coordinate system. Top left corner is the origin.
# (0,0)
# ------------------------->
# |
# |
# |
# |
# |
# |
# |
# |
# V
tiff_file = '../input/bathymetry_sumatra.tiff'
gtiff = gt(tiff_file)
shape = gtiff.tif_shape
nx = shape[1]
ny = shape[0]
nz = 1
# Bounding box.
bbox = gtiff.tif_bBox
# Put origin on the bottom left corner.
ox = bbox[0][0]
oy = bbox[1][1]
oz = 0.0
name = 'Elevation'
# Spacing in degrees
hd = 1.0/60.0
# Read data.
zarr = gtiff.read()
# Convert data to float32 array.
z = np.array(np.flipud(zarr),dtype='float32')

# Open file to write.
vti_file = 'bathymetry_sumatra.vti'
vf = open(vti_file, "wb")

# Calculate file layout information.
offset=0
ncomp=1
byte_int=4
byte_float=4
nbyte=nx*ny*nz*ncomp*byte_float
nbyte=nbyte+byte_int
bnbyte = np.int32(nbyte)

# Write VTI header.
vf.write("<?xml version=\"1.0\"?>\n".encode('ascii'))
vf.write("<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n".encode('ascii'))
vf.write("<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%.6f %.6f %.6f\" Spacing=\"%.6f %.6f %.6f\">\n".encode('ascii') % (nx-1,ny-1,nz-1,ox,oy,oz,hd,hd,hd))
vf.write("<Piece Extent=\"0 %d 0 %d 0 %d\">\n".encode('ascii') % (nx-1,ny-1,nz-1))
vf.write(f"<PointData Scalars=\"{name}\">\n".encode('ascii'))
vf.write(f"<DataArray type=\"Float32\" Name=\"{name}\" format=\"appended\" offset=\"{offset}\" />\n".encode('ascii'))
vf.write("</PointData>\n".encode('ascii'))
vf.write("<CellData>\n".encode('ascii'))
vf.write("</CellData>\n".encode('ascii'))
vf.write("</Piece>\n".encode('ascii'))
vf.write("</ImageData>\n".encode('ascii'))
vf.write("<AppendedData encoding=\"raw\">\n".encode('ascii'))
vf.write("_".encode('ascii'))

# Write total bytes (Binary).
vf.write(bnbyte)

# Write data (Binary).
vf.write(z)

# Write VTI footer.
vf.write("\n".encode('ascii'))
vf.write("</AppendedData>\n".encode('ascii'))
vf.write("</VTKFile>\n".encode('ascii'))
vf.close()
