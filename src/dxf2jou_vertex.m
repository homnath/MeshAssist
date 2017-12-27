%> @file dxf2jou_vertex.m
%> @brief Converts DXF AUTOCAD file to CUBIT/Trelis journal file.
%>
%> This function converts AUTOCAD 2000 (other ?) DXF ASCII file to a
%> CUBIT/Trelis journal file. The function extracts only
%> the 'VERTEX' tokens in the DXF file.
%>
%> <!-- @author Hom Nath Gharti (hgharti_AT_princeton_DOT_edu) -->
%>
%> ## Usage:
%>   dxf2jou_vertex(\em input_file) \n\n
%>   Example: \n
%>   dxf2jou('../input/dxf2jou_vertex_example.dxf') \n
%>
%> ## Input:
%>   input_file: DXF input file name.
%>
%> ## Output:
%>   .jou file which can be opened with CUBIT/Trelis.
% TODO: add creating polyline/s in journal file
%--------------------------------------------------------------------------

function dxf2jou_vertex(inpf_name)
if ~exist('inpf_name','var')
    error(strcat('dxf2vtk accepts exactly 1 argument!', ...
    ' Enter dxf file as an argument!'));
end
[path,fheader]=fileparts(inpf_name);
if isempty(path)
    path='.';
end

inpf=fopen(inpf_name,'r');
if inpf<=0
	error('File %s cannot be opened!\n',inpf_name);
end

% count number of 3D vertices
nvert=0;
fprintf(1,'counting vertices...');
while(~feof(inpf))
    if strcmp(fgetl(inpf),'VERTEX')
        nvert=nvert+1;
    end    
end
fclose(inpf);
fprintf(1,'complete!\n');
fprintf(1,'total number of vertices: %d\n',nvert);

fprintf(1,'extracting coordinates...');
inpf=fopen(inpf_name,'r');
coord=zeros(nvert,3);
ivert=0;
tline = fgetl(inpf);
while ischar(tline)
    if strcmp(tline,'VERTEX')
        ivert=ivert+1;
        %skip 3 lines
        fgetl(inpf);
        fgetl(inpf);
        fgetl(inpf);
        coord(ivert,1)=fscanf(inpf,'%f\n',1);
        fgetl(inpf);
        coord(ivert,2)=fscanf(inpf,'%f\n',1);
        fgetl(inpf);
        coord(ivert,3)=fscanf(inpf,'%f\n',1);
        fgetl(inpf);
               
    end
    tline = fgetl(inpf);
end
if nvert ~= ivert
   error('number of vertices mismatch!'); 
end
fprintf(1,'complete!\n');

% Remove duplicates and renumber
fprintf(1,'removing duplicate nodes...');
ucoord=uunique(coord,'rows','stable');
clear coord;
fprintf(1,'complete!\n');

outf_name=strcat(path,'/',fheader,'_vertex','.jou');
fprintf(1,'writing CUBIT journal file %s ...',outf_name);
outf_cubit=fopen(outf_name,'w');
fprintf(outf_cubit,'create vertex %.6f %.6f %.6f\n',ucoord');
fclose(outf_cubit);
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
%==========================================================================