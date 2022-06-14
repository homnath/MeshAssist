% This program reads and plots GEOTIFF files, e.g., etopo1 TIFF file.
% Bathymetry and topography can be extraced from: 
% https://www.ncei.noaa.gov/maps/grid-extract/
% AUTHOR:
%   Hom Nath Gharti, Queen's University, Kingston, Ontario, Canada.
% HISTORY:
%   Jan 31, 2022: First created.
% RUN:
%   process_geotiff('tiff_file_name')
%   Example:
%   process_geotiff('bathymetry_sumatra.tiff')
% INPUT:
%   tiff_file: etopo TIFF file name.
%==========================================================================
function process_geotiff(tiff_file)
% Earth radius (km).
R_EARTH = 6371;
DEGTORAD = pi/180.;

% Check if input file entered.
if nargin==0
    error('no input file name!')
end

% Check if the file exists.
if ~exist(tiff_file,'file')
    error('file ''%s'' doesn''t exist!',tiff_file);
end

% Read TIFF file.
info = geotiffinfo(tiff_file);
z = double(imread(tiff_file, 'tif'));
nx = info.Width;
ny = info.Height;

% Plot TIFF file
figure
pcolor(z);
shading flat;
colorbar('horiz');
axis ij, axis equal;
title('Bathymetry');

% Print Grid dimensions
fprintf(1,'Grid dimensions:\n');
fprintf(1,'    nx: %d\n',nx);
fprintf(1,'    ny: %d\n',ny);

% Print origin of the model
fprintf(1,'Grid origin:\n');
fprintf(1,'    longref: %.4f\n',info.BoundingBox(1,1));
fprintf(1,'    latref: %.4f\n',info.BoundingBox(1,2));

% Grid interval:
fprintf(1,'Grid interval:\n');
fprintf(1,'    In Arc minutes: %.4f\n',1);
fprintf(1,'    In km         : %.4f\n',DEGTORAD*R_EARTH/60);

