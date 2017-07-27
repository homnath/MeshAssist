%> @file write_vti.m
%> @brief Writes 3D gridded data to VTK VTI file.
%> 
%>   This function writes the VTI binary file for structured grid data, such as 
%>   finite difference data and tomography data. The VTI file can be visualized 
%>   in ParaView (http://www.paraview.org/).
%> 
%> @author Hom Nath Gharti (Princeton University)
%> 
%> ## Input:
%>   fname   : output file name \n
%>   ox      : origin vector [ox oy oz] \n
%>   dh      : sampling interval vector [dx dy dz] \n
%>   nx      : grid number vector [nx ny nz] \n
%>   name    : output variable name
%> 
%> ## Usage:
%>   call this functions with appropriate variables.
%> 
%> ## Notes:
%> 
%> - For a BigEndian architecture, replace "LittleEndian" with "BigEndian"
%>     in the line below.
%------------------------------------------------------------------------------
function write_vti(fname,ox,dh,nx,data,name)

outf=fopen(fname,'w');
fprintf(outf,'<?xml version="1.0"?>\n');
fprintf(outf,'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n');
fprintf(outf,'<ImageData WholeExtent="0 %d 0 %d 0 %d" Origin="%.6f %.6f %.6f" Spacing="%.6f %.6f %.6f">\n'...
    ,nx(1)-1,nx(2)-1,nx(3)-1,ox(1),ox(2),ox(3),dh(1),dh(2),dh(3));
fprintf(outf,'<Piece Extent="0 %d 0 %d 0 %d">\n',nx(1)-1,nx(2)-1,nx(3)-1);
fprintf(outf,'<PointData Scalars="%s">\n',name);

offset=0;
ncomp=1;
byte_int=4;
byte_float=4;
nbyte=nx(1)*nx(2)*nx(3)*ncomp*byte_float;
nbyte=nbyte+byte_int;

fprintf(outf,'<DataArray type="Float32" Name="%s" format="appended" offset="%d" />\n',name,offset);
fprintf(outf,'</PointData>\n');
fprintf(outf,'<CellData>\n');
fprintf(outf,'</CellData>\n');
fprintf(outf,'</Piece>\n');
fprintf(outf,'</ImageData>\n');
fprintf(outf,'<AppendedData encoding="raw">\n');
fprintf(outf,'_');
fwrite(outf,nbyte,'int32');
fwrite(outf,data,'float32');
fprintf(outf,'\n');
fprintf(outf,'</AppendedData>\n');
fprintf(outf,'</VTKFile>\n');
fclose(outf);
