% FILE
%   dxf2jou.m
% PURPOSE
%   This function converts AUTOCAD 2000 (other ?) DXF ASCII file to a
%   CUBIT/Trelis journal file, and optionally to a ParaView/VTK VTU ASCII
%   file --- a unstructured mesh files (.vtu). The function extracts only
%   the faces represented by 'AcDbFace' tokens in the DXF file.
% RUN
%   dxf2jou(inpfname, [save_vtu])
%   Example:
%   dxf2jou('../input/dxf2jou_example.dxf')
%   OR
%   dxf2jou('../input/dxf2jou_example.dxf',1)
% INPUT
%   inpfname: DXF input file name
% OPTION
%   An optional argument can be provided to save VTU
%   ASCII file (0: No [DEFAULT] or 1: yes) .
%   
% OUTPUT
%   .vtu file which can be visualized with ParaView/VTK
% AUTHOR
%   Hom Nath Gharti
%   formerly ay Institute of Engineering, Tribhuvan University, Nepal
%   formerly at NORSAR
%   Department of Geosciences, Princeton University, USA
%   hngharti_AT_gmail_DOT_com
% REVISION
%   May 01,2017: Hom Nath Gharti
%   - Clean up
%   Feb 15,2010: Hom Nath Gharti
%   - First created
%--------------------------------------------------------------------------

function dxf2jou(inpf_name)
if ~exist('inpf_name','var')
    error('dxf2vtk accepts exactly 1 argument! Enter dxf file as an argument!');
end
[path,fheader]=fileparts(inpf_name);
if isempty(path)
    path='.';
end

inpf=fopen(inpf_name,'r');
if inpf<=0
	error('File %s cannot be opened!\n',inpf_name);
end

% count number of 3D faces
nface=0; nline=0;
fprintf(1,'counting faces...');
while(~feof(inpf))    
    if strcmp(fgetl(inpf),'AcDbFace')
        nface=nface+1;        
    end
    if nface==0
        nline=nline+1;
    end
end
fclose(inpf);
fprintf(1,'complete!\n');
fprintf(1,'total number of faces: %d\n',nface);

fprintf(1,'extracting coordinates...');
coord=zeros(3*nface,3);
inpf=fopen(inpf_name,'r');
for i_line=1:nline; fgetl(inpf); end
for i_face=1:nface
    inode=[1 2 3]+(i_face-1)*3;
    fgetl(inpf);
    strblk=textscan(inpf,'%s',18);    
    coord(inode,:)=reshape(str2double(strblk{1}(2:2:18)),[3,3])';
    for i_line=1:22; fgetl(inpf); end
end
fprintf(1,'complete!\n');

fprintf(1,'removing duplicate nodes...');
% Remove duplicates and renumber
[ucoord, m, n]=uunique(coord);
clear coord;
fprintf(1,'complete!\n');

fprintf(1,'assigning face connectivity...');
face=int32(1:3*nface)';
for i_m=1:length(m)
    face(n==i_m)=i_m;    
end
fprintf(1,'complete!\n');

outf_name=strcat(path,'/',fheader,'_mesh','.jou');
fprintf(1,'writing CUBIT journal file %s ...',outf_name);
outf_cubit=fopen(outf_name,'w');
fprintf(outf_cubit,'create vertex %.6f %.6f %.6f\n',ucoord');
fprintf(outf_cubit,'create surface vertex %d %d %d\n',face);
fclose(outf_cubit);
fprintf(1,'complete!\n');

outf_name=strcat(path,'/',fheader,'_mesh','.vtu');
fprintf(1,'writing VTK file %s ...',outf_name);
write_vtu(ucoord,face,5,outf_name);
fprintf(1,'complete!\n');
end

% This function finds unsorted unique    
function [b, im, in] = uunique(a)     
    [~, im, in] = unique(a, 'rows','first'); 
    if nargout > 2 
        [ia, tmp] = sort(im); 
        [~, in] = ismember(in, tmp);
        clear tmp
    else 
       im = sort(im); 
    end 
    b = a(ia,:); 
end

% This function writes a VTU ASCII file
function write_vtu(coord,connect,VTK_etype,outf_name)
outf=fopen(outf_name,'w');
if outf<0
    error('file %s cannot be opened!\n',outf_name);
end
nelmt=size(connect,1);
nnode=size(coord,1);
fprintf(outf, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(outf, '<UnstructuredGrid>\n');
fprintf(outf, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n',nnode,nelmt/3);
fprintf(outf, '<Points>\n');
fprintf(outf, '<DataArray type = "Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(outf, '%f %f %f\n',coord(:,[2,1,3])'); %coord');
fprintf(outf, '</DataArray>\n');
fprintf(outf, '</Points>\n');
fprintf(outf, '<Cells>\n');
fprintf(outf, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
fprintf(outf, '%d %d %d\n',connect-1);
fprintf(outf, '</DataArray>\n');
fprintf(outf, '<DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(outf,'%d\n',3*(1:nelmt));
fprintf(outf,'</DataArray>\n');
fprintf(outf, '<DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(outf, '%d\n',VTK_etype*ones(nelmt,1));
fprintf(outf, '</DataArray>\n');
fprintf(outf, '</Cells>\n');
fprintf(outf, '</Piece>\n');
fprintf(outf, '</UnstructuredGrid>\n');
fprintf(outf, '</VTKFile>\n');
fclose(outf);
end
%==============================================================================

