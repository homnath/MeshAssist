%> @file dem2vti.m
%> @brief Converts DEM image file to ASCII XYZ file.
%>
%>   This program converts DEM map (TIFF image format) to ASCII XYZ 
%>   file and optionally ParaView/VTK VTI file format according to the 
%>   parameters defined in the input file. \n
%>   Note: Choosing a relativel small sampling interval may freeze the 
%>   program due to 'surfl' function.
%>
%> <!-- @author  Hom Nath Gharti (Princeton University), LPG, and 
%>   Michael Roth (NORSAR) -->
%>
%> ## Usage:
%>   dem2vti(\em input_file, [\em output_path]) \n\n
%>   Example: \n 
%>   dem2vti('dem2vti_example.in') \n
%>   OR \n
%>   dem2vti('dem2vti_example.in','../output')
%>
%> ## Input:
%>   input_file: Name of DEM file.      
%>
%> ## Options:
%>   An optional argument which must be a legitimate path can be
%>   provided as the output path. Default path is the current path.
%>
%> ## Output:
%>   All output files will be saved in the output_path provided.
%>   If no output_path is provided, current path is used.
%------------------------------------------------------------------------------

function dem2vti(infname,varargin)

% output path
if nargin == 1
    outpath=''; % current folder
elseif nargin == 2
	outpath=varargin{1};
else
	error('Invalid number of arguments!');
end
if ~isempty(outpath) && outpath(end) ~= '/'
    outpath=strcat(outpath,'/');
end
small=1e-12;

% Air properties (SI Units)
vp_air=300; 
vs_air=small;
ro_air=600; % Heavy air

% Read main input file
fprintf(1,'reaing input file...');
infmain=fopen(infname,'r');
fgetl(infmain);
infdem=fscanf(infmain,'%s\n',1);
fgetl(infmain);
ox=fscanf(infmain,'%f\n',1);
oy=fscanf(infmain,'%f\n',1);
oz=fscanf(infmain,'%f\n',1);
fgetl(infmain);
dh=fscanf(infmain,'%f\n',1);
fgetl(infmain);
nx=fscanf(infmain,'%d\n',1);
ny=fscanf(infmain,'%d\n',1);
nz=fscanf(infmain,'%d\n',1);
fgetl(infmain);
vp_rock=fscanf(infmain,'%f\n',1);
vs_rock=fscanf(infmain,'%f\n',1);
ro_rock=fscanf(infmain,'%f\n',1);
fgetl(infmain);
plotfig=fscanf(infmain,'%d\n',1);
fgetl(infmain);
save_xyz=fscanf(infmain,'%d\n',1);
fgetl(infmain);
save_vti2d=fscanf(infmain,'%d\n',1);
fgetl(infmain);
save_vti3d=fscanf(infmain,'%d\n',1);
fclose(infmain);
fprintf(1,'complete!\n');

[inpath,fheader]=fileparts(infname);
if inpath(end) ~= '/'
    inpath=strcat(inpath,'/');
end

% Compute extent
xmax=ox+(nx-1)*dh;
ymax=oy+(ny-1)*dh;
zmax=oz+(nz-1)*dh;

% Read DEM file
fprintf(1,'reaing DEM file...');
demfname=strcat(inpath,infdem);
info = geotiffinfo(demfname);
Z0 = double(imread(demfname, 'tif'));

% Analyse data
nxb = info.Width;
nyb = info.Height;
%demsize = nxb * nyb;
%effsize = length(find(DEM));
%pct = 100 * effsize / demsize;

if ~isfield(info.CornerCoords,'PCSX')
 oxb=info.CornerCoords.X(1);
 oyb=info.CornerCoords.Y(3);
else
    oxb=info.CornerCoords.PCSX(2);
    oyb=info.CornerCoords.PCSY(2);
end
fprintf(1,'complete!\n');
fprintf(1,'Origin: %f,%f,%f\n',oxb,oyb,min(min(Z0)));
fprintf(1,'Bounding box: %f,%f,%f\n',info.BoundingBox(2,1),...
    info.BoundingBox(2,2),max(max(Z0)));

% Check extent
if ox<oxb || oy<oyb || oz<min(min(Z0)) || ...
xmax>info.BoundingBox(2,1) || ymax>info.BoundingBox(2,2) || zmax>max(max(Z0))
    error('Defined model is outside the range!');
end

fprintf(1,'downsampling data...');
% Downsample data
[X, Y] = meshgrid((0:nxb-1)+oxb,(0:nyb-1)+oyb);
Z0 = flipud(Z0);
[Xi, Yi] = meshgrid((0:dh:(nx-1)*dh)+ox,(0:dh:(ny-1)*dh)+oy);

Z = interp2(X,Y,Z0,Xi,Yi);
clear Z0

%B=flipud(B);
% Avoid out of memory
if numel(Z) > 10e6
    disp('Too many points');
    return;
end
fprintf(1,'complete!\n');

% Extract elevation in 1D array
zvec=reshape(Z',1,nx*ny);

if plotfig==1
    fprintf(1,'plotting figure...');
    figure;
    set(gcf,'Units','pixels','Position',[10 10 1600 1000])
    surfl(Xi, Yi, Z);
    lighting phongpwd
    shading interp;
    set(gca, 'Box', 'on');
    axis image;
    grid on;
    colormap bone;
    view(15, 15);   % Optionnal
    view(0, 90)
    xticks=get(gca,'xtick');
    nxtick=length(xticks);
    xticklabel=cell(1,nxtick);
    for i=1:nxtick
        xticklabel{i}=num2str(xticks(i)-xticks(1));
    end
    set(gca,'xticklabel','')
    xstring=sprintf('UTM %6d - %6d', xticks(1),xticks(end));
    xlabel(xstring)
    yticks=get(gca,'ytick');
    nytick=length(yticks);
    yticklabel=cell(1,nytick);
    for i=1:nytick
        yticklabel{i}=num2str(yticks(i)- yticks(1));
    end
    set(gca,'yticklabel','')
    ystring=sprintf('UTM %7d - %7d', yticks(1),yticks(end));
    ylabel(ystring)
         
    view(2);
    set(gcf, 'Renderer', 'OpenGL');
    alpha(.9)

    xlabel(xstring)
    ylabel(ystring)
    colormap(bone)
    fprintf(1,'complete!\n');
end

% write XYZ ASCII file
if save_xyz ==1    
    fprintf(1,'saving XYZ file...');
    nxy=nx*ny;
    outxyz=strcat('outpath',fheader,'_elevation.xyz');
    outf=fopen(outxyz,'w');
    fprintf(outf,'%f %f %f\n',...
        [reshape(Xi,1,nxy); reshape(Yi,1,nxy); reshape(Z,1,nxy)]);
    fclose(outf);
    fprintf(1,'complete!\n');
end

% Write 2D .vti file to be viewed in vtk/ParaView
if save_vti2d == 1
    fprintf(1,'saving 2D VTI file...');
    outvti=strcat(outpath,fheader,'_elevation.vti');
    write_vti(outvti,[ox oy oz],[dh dh dh],[nx ny 1],zvec,'Elevation');
    fprintf(1,'complete!\n');
end

if save_vti3d == 1
    fprintf(1,'craeating & saving 3D VTI files...');

    % Create model
    vp(1:nx*ny*nz,1)=vp_rock;
    vs(1:nx*ny*nz,1)=vs_rock;
    ro(1:nx*ny*nz,1)=ro_rock;
    %zslope(1:nx*ny,1)=0;

    addtop=0; % Add extra thickness at the top for topography
    % TODO: this loop can be make much faster?
    for iz=1:nz
        z=oz+(nz-iz)*dh; % Start from the top left corner that will be used as the origin
        for iy=1:ny
            for ix=1:nx
                if z > (Z(iy,ix)+addtop)                  
                    ind=(iz-1)*ny*nx+(iy-1)*nx+ix;
                    vp(ind)=vp_air;
                    vs(ind)=vs_air;
                    ro(ind)=ro_air;
                end
            end
        end
    end
    clear Z

    % Shift the origin to top left corner
    oz_top=oz+(nz-1)*dh;

    % Write P velocity model
    outvti=strcat(outpath,fheader,'_vp.vti');
    write_vti(outvti,[ox oy oz_top],[dh dh -dh],[nx ny nz],vp,'Vp');

    % Write S velocity model
    outvti=strcat(outpath,fheader,'_vs.vti');
    write_vti(outvti,[ox oy oz_top],[dh dh -dh],[nx ny nz],vs,'Vs');

    % Write density model
    outvti=strcat(outpath,fheader,'_rho.vti');
    write_vti(outvti,[ox oy oz_top],[dh dh -dh],[nx ny nz],ro,'Rho');
    fprintf(1,'complete!\n');
end
disp('-------------------------------------------------');
%==============================================================================
