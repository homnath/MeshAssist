import rasterio as rio
from rasterio.plot import show

img = rio.open('bathymetry_sumatra.tiff')
show(img)