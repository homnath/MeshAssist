% This function converts AUTOCAD 2000 (other ?) DXF files to VTK
% unstructured mesh files (.vtu). The function extracts only the faces
% represented by 'AcDbFace' tokens.
% RUN
% >> dxf2vtk(inpfname)
% INPUT
%   inpfname: DXF input file name
% OUTPUT
%   .vtu file which can be visualized with ParaView/VTK
% AUTHOR
%   Hom Nath Gharti
%   NORSAR
%   homnath_AT_norsar_DOT_com
% REVISION
%   Feb 15,2010
function dxf2vtk_inkscape(inpf_name)
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

% count total number of polylines
npoly=0; nline=0;
fprintf(1,'counting polylines...');
while(~feof(inpf))
    if strcmp(fgetl(inpf),'POLYLINE')
        npoly=npoly+1;
    end
    if npoly==0
        nline=nline+1;
    end
end
fclose(inpf);
fprintf(1,'complete!\n');
fprintf(1,'total number of polylines: %d\n',npoly);

nvert=zeros(npoly,1);
nskip=zeros(npoly,1);
% count total number of vertices
ipoly=0;
fprintf(1,'counting vertices...');
inpf=fopen(inpf_name,'r');
ncount=0;
iscount=0;
while(~feof(inpf))
    sline=fgetl(inpf);
    if strcmp(sline,'POLYLINE')
        ipoly=ipoly+1;
        ncount=0; %
        iscount=1;
    end
    if iscount
        ncount=ncount+1;
    end
    if strcmp(sline,'VERTEX')
        nvert(ipoly)=nvert(ipoly)+1;
        nskip(ipoly)=ncount;
        iscount=0;
    end
end
fclose(inpf);
% correct the number of lines from POLYLINE just before the first 'VERTEX'
nskip=nskip-1;

fprintf(1,'complete!\n');
fprintf(1,'total number of vertices in each polylines: %d\n',nvert);

fprintf(1,'extracting coordinates...');

outf_name=strcat(path,'/',fheader,'.jou');
fprintf(1,'writing CUBIT journal file %s ...',outf_name);
outf_cubit=fopen(outf_name,'w');

iv0=1; % numbering starts from 1
inpf=fopen(inpf_name,'r');
for i_line=1:nline; fgetl(inpf); end % skip until the first POLYLINE
for i_poly=1:npoly
    coord=zeros(nvert(i_poly),3);
%     while ~strcmp(sline,'VERTEX')
%        sline=fgetl(inpf) 
%     end
    % skip until just before the first VERTEX after POLYLINE
    for i_line=1:nskip(i_poly); fgetl(inpf); end    
    for i_vert=1:nvert(i_poly)
        fgetl(inpf); % VERTEX
        fgetl(inpf);
        fgetl(inpf);

        fgetl(inpf);
        coord(i_vert,1)=fscanf(inpf,'%f\n',1);
        fgetl(inpf);
        coord(i_vert,2)=fscanf(inpf,'%f\n',1);
        fgetl(inpf);
        coord(i_vert,3)=fscanf(inpf,'%f\n',1);

        fgetl(inpf);
        fgetl(inpf);
        fgetl(inpf);
        fgetl(inpf);
        fgetl(inpf);
        %------------------------------------------------------------------
    end
    fgetl(inpf);
    fgetl(inpf);
    fgetl(inpf);
    fgetl(inpf);

    fprintf(1,'removing duplicate nodes...');
    % Remove duplicates and renumber
    [ucoord, m, ~]=uunique(coord);
    clear coord;
    fprintf(1,'complete!\n');

    nuvert=length(m); % number of unique vertices


    fprintf(outf_cubit,'create vertex %.6f %.6f %.6f\n',ucoord');
    fprintf(outf_cubit,'create curve spline vertex %d to %d %d\n',iv0,iv0+nuvert-1,iv0);
    iv0=iv0+nuvert;
end
fclose(inpf);
fprintf(1,'complete!\n');



fclose(outf_cubit);
fprintf(1,'complete!\n');
end

% Function unsorted unique
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
