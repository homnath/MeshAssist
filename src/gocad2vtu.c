/** @file gocad2vtu.c 
*   @brief Converts GOCAD ASCII file to VTU file.
*
*  This program converts the GOCAD ASCII file (3-noded triangular meshes)
*  to VTK XML .vtu binary file (unstructured mesh file) which can be 
*  visualized/processed in ParaView or VTK.
*
*  <!-- @author Hom Nath Gharti (hgharti_AT_princeton_DOT_edu) -->
*
* ##Dependencies:
*	stringmanip.c
*
* ## Compile
*  - in parent folder, type:
*  make
*  OR
*  - in src/ folder, type
*	gcc gocad2vtu.c -o gocad2vtu
*
* ## Usage:
*	./bin/gocad2vtu \em input_file [\em Options] \n\n
*	Example: ./bin/gocad2vtu ./input/gocad2vtu_example.ts
*
*  ##Options:
*	- -fac: Use this option to multiply the coordinates by a certain factor, this is 
*		helpful for unit conversion, e.g. for m to km use 0.001, for km to m 
*    use 1000, example: gocad2vtu T2_horizon.ts -fac=0.001
*
* ##Notes:
*	- Output .vtu file is binary, therefore endianness of the processor
*  architechture is important.
*
*	- This program automatically identify the endianness and write the output
*  accordingly. Hence if you run and process/visualize .vtu file in the
*  architecture with different endianness there may be an error. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

#define maxline 150 /* Max number of characters in a line */
#define maxenod 3 /* Max number of nodes per element */
#define maxvar 4 /* Max number of variable sets be plotted in ParaView file */
#define LE 0 /* Little Endian */
#define BE 1 /* Big Endian */

int main(int argc,char **argv)
{
	int i,j,etype,nnod,nelmt; 
	int endian,dumi,enod[maxenod];
	int ndim,nenod,nvar;
	int bytes[maxvar],off[maxvar],temp;
	float tempx,tempy,tempz,fac;
	char byte_order[12],buffer[maxline],dumc[10],fonly[62],outfname[62];
	
	FILE *inf,*outf0,*outf1; /* ,*outf_vtk; */
	
	ndim = 3; 	/* Number of dimensions */
	nenod = 3; 	/* Number of nodes per element */
	nvar = 4;
	etype = 5; /* VTK triangle */	

	printf("Processing...\n");
		
	fac=1.0;	
		
	/* Open input file */
	if(argc<2){
		fprintf(stderr,"ERROR: input file not entered!\n");
		exit(-1);
	}
	inf=fopen(argv[1],"rb");
	if(inf==NULL){
		fprintf(stderr,"ERROR: file \"%s\" not found!",argv[1]);
		exit(-1);
	}
	
	/* Scan options */
	if(argc>2){
		for(i=2;i<argc;i++){
			if(matchfirstword(argv[i],"-fac")){
				getvalue(argv[i],"fac","f",(int *)&fac);
			}
			else{
				printf("ERROR: unrecognized option \"%s\"",argv[i]);
				exit(-1);
			}
		}
	}
	
	/* Processor endianness */
	endian=getEndian();
	if(endian == LE){
		strcpy(byte_order,"LittleEndian");
	}
	else if(endian == BE){
		strcpy(byte_order,"BigEndian");
	}
	else{
		printf("ERROR: illegal endianness!\n");
		exit(-1);
	}	
	
	printf("--------------------------------\n");
	printf("Input file: %s\n",argv[1]);
	printf("Coordinates multiplication factor: %f\n",fac);	
	printf("--------------------------------\n");
	
	outf0=fopen("coord_temp","wb");
	outf1=fopen("connect_temp","wb");
	nnod=0; nelmt=0;
	while(fgets(buffer,maxline,inf) != NULL){
		if(matchfirstword(buffer,"VRTX") || matchfirstword(buffer,"PVRTX")){
			nnod++;			
			sscanf(buffer,"%s %d %f %f %f",dumc,&dumi,&tempx,&tempy,&tempz);
			if(fac != 1.0){
				tempx*=fac; tempy*=fac; tempz*=fac;
			}
			fwrite(&tempx,sizeof(float),1,outf0);
			fwrite(&tempy,sizeof(float),1,outf0);
			fwrite(&tempz,sizeof(float),1,outf0);
		}
		
		if(matchfirstword(buffer,"TRGL")){
			nelmt++;	
			sscanf(buffer,"%s %d %d %d",dumc,&enod[0],&enod[1],&enod[2]);
			fwrite(&enod[0],sizeof(int),1,outf1);
			fwrite(&enod[1],sizeof(int),1,outf1);
			fwrite(&enod[2],sizeof(int),1,outf1);
		}
	}
	fclose(inf);
	fclose(outf0);
	fclose(outf1);
	if(nnod<=0 || nelmt<=0){
		printf("ERROR: zero entities [nodes=%d, elements=%d]!\n",nnod,nelmt);
		exit(-1);
	}
	printf("Number of nodes: %d\n",nnod);
	printf("Number of elements: %d\n",nelmt);
	
	printf("--------------------------------\n");
	
	bytes[0] = (ndim*nnod)*sizeof(float); 	/* Coordinates */
	bytes[1] = (nenod*nelmt)*sizeof(int); 	/* Connectivity */
	bytes[2] = (nelmt)*sizeof(int); 		/* Offsets */
	bytes[3] = (nelmt)*sizeof(int); 		/* Types */
	
	off[0]=0; /* 1st offset */
	for (i=0; i<nvar; i++){		
		if(i<nvar-1)off[i+1]=off[i]+sizeof(int)+bytes[i];
		bytes[i]=bytes[i]+sizeof(int);
	}
	
	removeExtension(argv[1],fonly);	
	sprintf(outfname,"%s_mesh.vtu",fonly);
	
	outf0=fopen(outfname,"wb");	
	
	fprintf(outf0,"<?xml version=\"1.0\"?>\n");
	fprintf(outf0,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",byte_order);
	fprintf(outf0,"<UnstructuredGrid>\n");
	fprintf(outf0,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",nnod,nelmt);
	fprintf(outf0,"<Points>\n");
	fprintf(outf0,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",off[0]);
	fprintf(outf0,"</Points>\n");
	fprintf(outf0,"<PointData>\n");
	fprintf(outf0,"</PointData>\n");		
	fprintf(outf0,"<Cells>\n");
	fprintf(outf0,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\"/>\n",off[1]);
	fprintf(outf0,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\"/>\n",off[2]);
	fprintf(outf0,"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%d\"/>\n",off[3]);
	fprintf(outf0,"</Cells>\n");
	fprintf(outf0,"<CellData>\n");
	fprintf(outf0,"</CellData>\n");
	fprintf(outf0,"</Piece>\n");
	fprintf(outf0,"</UnstructuredGrid>\n");
	fprintf(outf0,"<AppendedData encoding=\"raw\">\n");
	fprintf(outf0,"_");
	
	inf=fopen("coord_temp","rb");
	
	/* outf_vtk=fopen("test.vtk","w");
	fprintf(outf_vtk,"# vtk DataFile Version 2.0\n");
	fprintf(outf_vtk,"Unstructured Grid Example\n");
	fprintf(outf_vtk,"ASCII\n");
	fprintf(outf_vtk,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(outf_vtk,"POINTS %d float\n",nnod);*/

	/* Coordinates */
	fwrite(&bytes[0],sizeof(int),1,outf0);
	for(i=0;i<nnod;i++){
		for(j=0;j<ndim;j++){
			fread(&tempx,sizeof(float),1,inf);
			fwrite(&tempx,sizeof(float),1,outf0);
			/*fprintf(outf_vtk,"%.6f ",tempx);*/
			/*if(j==ndim-1)fprintf(outf_vtk,"\n");*/
		}
	}
	fclose(inf);
	printf("Coordinates: SUCCESS\n");	
	
	inf=fopen("connect_temp","rb");
	/*fprintf(outf_vtk,"CELLS %d %d\n",nelmt,nelmt*(nenod+1));*/
	/* Connectivity */
	fwrite(&bytes[1],sizeof(int),1,outf0);
	for(i=0;i<nelmt;i++){
		for(j=0;j<nenod;j++){
			fread(&enod[j],sizeof(int),1,inf);
			enod[j]--;
			fwrite(&enod[j],sizeof(int),1,outf0);
			/*fprintf(outf_vtk,"4 %d %d %d\n",enod[0],enod[1],enod[2]);*/
		}
	}
	fclose(inf);
	remove("coord_temp");
	remove("connect_temp");
	
	/* Offsets */
	fwrite(&bytes[2],sizeof(int),1,outf0);
	temp=0;
	for(i=0;i<nelmt;i++){
		temp+=3;
		fwrite(&temp,sizeof(int),1,outf0);

	}
	printf("Offsets: SUCCESS\n");
	
	/* Types */
	/*fprintf(outf_vtk,"CELL_TYPES %d\n",nelmt);*/
	fwrite(&bytes[3],sizeof(int),1,outf0);
	for(i=0;i<nelmt;i++){
		fwrite(&etype,sizeof(int),1,outf0);
		/*fprintf(outf_vtk,"%d\n",etype);*/
		
	}
	/*fclose(outf_vtk);*/
	printf("Types: SUCCESS\n");
	
	fprintf(outf0,"\n");
	fprintf(outf0,"</AppendedData>\n");
	fprintf(outf0,"</VTKFile>\n");
	fclose(outf0);	
	
	printf("Status: SUCCESS\n");
	printf("--------------------------------\n");	

	return(0);
}
/*=============================================================================*/
