# Convert a geotiff (.tiff/.tif) into a vertical list of xyz coordinates
# Can be used with the xyz2jou.f90 to turn .geotiff into .jou for CUBIT use
#
# Rasterio module is easiest to access through Anaconda, or it can be installed 
# through pip
#
# utm module can can be accessed through anaconda or it can be installed with pip
#------------------------------------------------------------------------------
# Input:
# Only input is the geotiff file as 'infname'. Be careful with the input file 
# extension; geotiffs can end with '.tif' of '.tiff' which will effect reading.
# Provide file to convert:
infname = 'bathymetry_sumatra.tiff'
#
# Output:
# User must provide a name for the vti file, as well as the variable of interest
# Provide name for output file
#outFile = 'bathymetry_sumatra_xyz.utm'
#==============================================================================
import rasterio as rio
import math as mt
import numpy as np
import utm
import os
import argparse

def process_file(infname, write_vti, write_xyz):

    # Check if the input file has a .txt extension
    if not infname.endswith('.tiff'):
        print("Error: Input file must have a .tiff extension.")
        return

    # Check if the input file exists
    if not os.path.exists(infname):
        print(f"Error: File '{infname}' not found.")
        return

    print('File to process: %s' % infname)    
    
    # Set constants
    rEarth = 6371 
    degToRad = mt.pi/180
    m2km = 0.001 # *May need to adujst this parameter* #
    
    # Open the file as a geotiff and as an array
    img = rio.open(infname)
    data = img.read(1)    
    dataFull = img.read()
    
    # Get origin data
    ox = img.bounds.left
    oy = img.bounds.bottom
    oz = 0.0
    
    # Convert the long & lat into UTM
    utms = utm.from_latlon(oy,ox)
    oxUTM = utms[0]
    oyUTM = utms[1]
    ozUTM = oz
        
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
    dz = 1 # May need to modify to run 3D maps
    
    print(f"{dx:.9f}, {dy:.9f}")
    # Get the number of points in the data
    nx = dataFull.shape[2]
    ny = dataFull.shape[1]
    nz = dataFull.shape[0]
    
    # If write_option is 1, write the content to a .out file
    if write_vti == 1 or write_xyz == 1:
        # Get the base name without extension
        # base_fname = os.path.splitext(infname)[0]
        # Get the base name without extension and path
        base_fname = os.path.splitext(os.path.basename(infname))[0]
        
    if write_vti == 1:
        z = np.array(np.flipud(data),dtype='float32')
        # Output variable name
        outVar = 'Elevation'
        # Create a new file and start writing to it 
        outfile = base_fname + '.vti'
        print('Writing a file: %s' % outfile)  
        # Create a new file and start writing to it
        # 'wb' opens the file in binary and then use '.encode(ascii)' to write ascii
        outf = open(outfile, "wb")
    
        # Get the number of bytes of the actual vti data (excludes metadata)
        offset=0
        ncomp=1
        byte_int=4
        byte_float=4
        nbyte=nx*ny*nz*ncomp*byte_float
        nbyte=nbyte+byte_int
        bnbyte = np.int32(nbyte)

        # Write ascii
        outf.write('<?xml version=\"1.0\"?>\n'.encode('ascii'))
        outf.write('<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n'.encode('ascii'))
        outf.write('<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%.6f %.6f %.6f\" Spacing=\"%.6f %.6f %.6f\">\n'.encode('ascii') % (nx-1,ny-1,nz-1,oxUTM,oyUTM,ozUTM,dx,dy,dz))
        outf.write('<Piece Extent=\"0 %d 0 %d 0 %d\">\n'.encode('ascii') % (nx-1,ny-1,nz-1))
        outf.write(f'<PointData Scalars=\"{outVar}\">\n'.encode('ascii'))
        outf.write(f'<DataArray type="Float32" Name=\"{outVar}\" format=\"appended\" offset=\"{offset}\" />\n'.encode('ascii'))
        outf.write('</PointData>\n'.encode('ascii'))
        outf.write('<CellData>\n'.encode('ascii'))
        outf.write('</CellData>\n'.encode('ascii'))
        outf.write('</Piece>\n'.encode('ascii'))
        outf.write('</ImageData>\n'.encode('ascii'))
        outf.write('<AppendedData encoding=\"raw\">\n'.encode('ascii'))
        outf.write('_'.encode('ascii'))

        # Write binary        
        outf.write(bnbyte)
        outf.write(z)

        # Write ascii
        outf.write('\n'.encode('ascii'))
        outf.write('</AppendedData>\n'.encode('ascii'))
        outf.write('</VTKFile>\n'.encode('ascii'))
        outf.close()
    
    if write_xyz == 1:
        z = np.array(np.fliplr(data),dtype='float32')
        # Create a new file and start writing to it 
        outfile = base_fname + '.xyz'
        print('Writing a file: %s' % outfile)
        outf = open(outfile, 'w')
        
        for j in range(ny):
            x = oxUTM
            y = oyUTM + (dy*j)
            for i in range(nx):
                outf.write('%f      %f      %f' % (x+(dx*i),y,z[j,i]*m2km))
                if i!= nx-1:
                    outf.write('\n')
            if j != ny-1:
                outf.write('\n')
        outf.close()        
    
#==============================================================================

# Main routine
if __name__ == "__main__":
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Process a .txt file and optionally write to a .out file.")
    parser.add_argument('infname', type=str, help="Input file name with .txt extension")
    parser.add_argument('-vti', type=int, choices=[0, 1], default=0, help="1 to write VTI file, 0 to skip")
    parser.add_argument('-xyz', type=int, choices=[0, 1], default=0, help="1 to write XYZ file, 0 to skip")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with provided arguments
    process_file(args.infname, args.vti, args.xyz)
#==============================================================================