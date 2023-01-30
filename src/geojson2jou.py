#==============================================================================    
# This program converts the geojson file to CUBIT Journal file.
# DEVELOPER:
#   Hom Nath Gharti
#   Department of Geological Sciences and Geological Engineering.
#   Queen's University, Kingston, ON, Canada.
# USAGE: 
#   python geojson2jou.py -f geojson_file_name
#   EXAMPLE: python geojson2jou.py -f ../input/buildings.geojson
# OPTIONS:
#   -f geojson_file_name, --file geojson_file_name
#   -u convert to UTM coordinate system, --utm geojson_file_name 
# HISTORY:
#   Jan 30,2023: Fixed bugs and added options.    
#   Jan 27,2023: First created.    
#==============================================================================    
import os
import sys
import json
import utm
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str)
parser.add_argument("-u", "--utm", type=int, default=0)
args = parser.parse_args()

if args.file == None:
  exit('ERROR: You must enter the file name using -f or --file option!')

# Open the file.
inpfname = args.file
with open(inpfname) as data:
    geo = json.load(data)

print('Input file: "%s"' % inpfname)
if args.utm > 0:
  print('Converting to UTM coordinates.')

# Split file name to extract the base file name.    
split_fname = os.path.splitext(os.path.basename(inpfname))

# CUBIT Journal file name.
outfname = split_fname[0]+'.jou'

i2=0
npolygon=0   
npoint=0 
outf = open(outfname, "w")
for feature in geo['features']:        
    npolygon+=1
    for coord in feature['geometry']['coordinates']:        
        np=len(coord)
        npoint+=np-1
        i1=i2+1
        n = 0
        #i1 = 0
        for xy in coord:
            n+=1
            if n >= np:
                n=n-1
                break

            if args.utm > 0:
              [utm_east, utm_north, utm_zone, utm_letter] = utm.from_latlon(xy[0],xy[1])
              outf.write(f"create vertex {utm_east} {utm_north} {0}\n")
            else:
              outf.write(f"create vertex {xy[0]} {xy[1]} {0}\n")
        i2 = i1+n-1
        outf.write(f"create surface vertex {i1} to {i2}\n")

print('Output file: "%s"' % outfname)        
print('Total %d polygons and %d points written.\n' % (npolygon,npoint))
#==============================================================================    
