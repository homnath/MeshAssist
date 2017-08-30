/** @file vti2cell.c
*  @brief This file converts VTI file to VTU file.
*
*  This program converts the 2D/3D Binary VTK XML .vti file to unstructured 
*  mesh files (.vtu). This program also generates the mesh files required by 
*  SPECFEM2D and SPECFEM3D. Note that the file formats in SPECFEM2D and 
*  SPECFEM3D are different. This should be made same format as soon as 
*  possible. For this, source codes within the decompose folder of SPECFEM3D 
*  and cubit2specfem3d.py need to be changed. 
*
*  <!-- @author Hom Nath Gharti (hgharti_AT_princeton_DOT_edu) -->
*
* ## Dependencies:
*  stringmanip.c
*
* ## Compile:
*  gcc vti2cell.c -o vti2cell -lm
*
*  ## Usage: 
*  vti2cell \em input_file [\em Options] \n\n
*  Example: \n
*  vti2cell py_plane_model.vti
*
* ## Options:
* - -fac=factor (real)
*    Use this option to multiply the coordinates by certain factor, this is 
*    helpful for unit conversion, e.g. for m to km use 0.001, for km to m use 1000
*    Example: vti2cell2d py_plane_model.vti -fac=1000
* - -xmat=exclusion material id/s (integer/s)
*    Use this option to exclude certain region of the model, e.g. exclusion of air.
*    Appropriate id/s should be supplied, id s are number orderd according to the value 
*    of corresponding material properties and numbered starting from 1. This way,
*    lowest value will have id 1 and so on.
*    Example: vti2cell py_plane_model.vti -xmat=1,2
*    This command will exclude the regions with material id 1 and 2.
*      
*    Example: vti2cell py_plane_model.vti -fac=1000 -xmat=1
*    This command multiply the coordinates by 1000 and exclude the region with material id 1
* - -step=step size (integer)
*    Use this option to coarsen the mesh. This value represent the number of grids to be used as 1 element,
*    e.g., if you want to make 2 grids as 1 element, use -step=2
* - -zup=z axis direction indicator (integer)
*    Use this option to indicate whether the Z axis direction is up
*
* ## Toto:
* - make uniformity for 2D,3D, e.g., writing and reading coordinates
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

/* Bounds for some parameters */
#define maxline 150 /* Max number of characters in a line */
#define maxvar 5  /* Max number of variable sets to be plotted in ParaView */
#define maxmat 1000000  /* Max of materials */

#define ON 1
#define OFF 0
#define eps 1.0e-16

#define LE 0 /* Little Endian */
#define BE 1 /* Big Endian */

/* inline function to compare two float values a and b*/
int comp_float (const void *a, const void *b) {
   if (*(float *)a < *(float *)b) return -1;
   if (*(float *)a > *(float *)b) return 1;
   return 0;
}

int main(int argc,char **argv)
{
  int ch,i,ielmt,inode,j,k;
  int ii,ji,ki;
  int endian,etype,ienode; 
  int dumi,*enode;
  int ndim,nenode,nvar,nmat;
  int nx,ny,nz; /* number of grid points */
  int nex,ney,nez; /* number of elements */
  int nnx,nny,nnz; /* number of nodes */
  int next,bytes[maxvar],offset[maxvar];
  int nnode,nelmt,new_nnode,new_nelmt,inew_elmt;
  int **elmt_node,*nstat,*estat,*nmir,*elmt_mat; 
  int wx1,wx2,wy1,wy2,wz1,wz2; /* Whole extent */
  int px1,px2,py1,py2,pz1,pz2; /* Piece extent */
  int imat,emat_id,match;
  int i_char,matid,nxmat,nchar,temp,*xmat,*new_matid;
  int step; /* coarsening factor */
  int xup,yup,zup; /* switch to indicate whether the z direction is up */
  int ioff,off0;
  /* offset from off0, offset from beginning of file to where binary grid data
     starts in vti */
  int idatum,isign,jdatum,jsign,kdatum,ksign;
  /* controllers of k index according to the z axis orientation (up or down)*/
  float mat[maxmat];
  float ox,oy,oz; /* Origin */
  float dx,dy,dz; /* Spacig */
  float tempx,tempy,tempz;  
  float fac;  
  char byte_order[12],buffer[maxline],string[maxline],stag[maxline];
  char fonly[62],outfname[62],vname[62],vtype[62];
  FILE *inf,*outf0,*outf1,*outf_abs,*outf_xmin,*outf_xmax,*outf_ymin,         \
  *outf_ymax,*outf_zmin,*outf_zmax;
  
  /* Actual parameters */   
  nvar = 5; /* Number of variable sets to be plotted in ParaView */  
    
  fac=1.0; nxmat=0; step=1; 
  xup=yup=zup=1;
  idatum=jdatum=kdatum=0;
  isign=jsign=ksign=1; /* default is z up */
    
  /* Open input file */
  if(argc<2){
    fprintf(stderr,"ERROR: input file not entered!\n");
    exit(-1);
  }
  inf=fopen(argv[1],"rb");
  if(inf==NULL){
    fprintf(stderr,"ERROR: file \"%s\" not found!\n",argv[1]);
    exit(-1);
  }
  
  if(argc>2){
    for(i=2;i<argc;i++){
      if(matchfirstword(argv[i],"-fac")){
        getvalue(argv[i],"fac","f",(int *)&fac);
      }
      else if(matchfirstword(argv[i],"-xmat")){
        printf("%s\n",argv[i]);
		/* count number of material regions to exclude */
		nxmat=1; /* default is at least 1 */
		nchar=strlen(argv[i]);
		for(i_char=0;i_char<nchar;i_char++){
			if(argv[i][i_char]==','){
			nxmat++;
			}
		}				
		xmat=malloc(nxmat*sizeof(int));
		getintegervect(argv[i],"=",nxmat,xmat);
    }
      else if(matchfirstword(argv[i],"-step")){
        getvalue(argv[i],"step","d",&step);
        if(step<=0){
          printf("WARNING: \"-step\" must be a positive integer! changed to default value 1!\n");
        }
      }
	  else if(matchfirstword(argv[i],"-xup")){
        getvalue(argv[i],"xup","d",&xup);
        if(xup == 0 || xup == 1){
          continue;
        }
        else{
          xup=1;
          printf("WARNING: \"-xup\" must either be 0 or 1! changed to default value 1!\n");
        }
      }
	  else if(matchfirstword(argv[i],"-yup")){
        getvalue(argv[i],"yup","d",&yup);
        if(yup == 0 || yup == 1){
          continue;
        }
        else{
          yup=1;
          printf("WARNING: \"-yup\" must either be 0 or 1! changed to default value 1!\n");
        }
      }
      else if(matchfirstword(argv[i],"-zup")){
        getvalue(argv[i],"zup","d",&zup);
        if(zup == 0 || zup == 1){
          continue;
        }
        else{
          zup=1;
          printf("WARNING: \"-zup\" must either be 0 or 1! changed to default value 1!\n");
        }
      }
      else{
        printf("ERROR: unrecognized option \"%s\"",argv[i]);
        exit(-1);
      }
    }
  }

  /* processor endianness */
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
  printf("input file: %s\n",argv[1]);
  printf("coordinates multiplication factor: %f\n",fac);  
  printf("coursening factor: %d\n",step);
  if(xup==1){
    printf("x direction: UP\n");
  }
  else {
    printf("x direction: DOWN\n");
  }
  if(yup==1){
    printf("y direction: UP\n");
  }
  else {
    printf("y direction: DOWN\n");
  }
  if(zup==1){
    printf("z direction: UP\n");
  }
  else {
    printf("z direction: DOWN\n");
  }
  if(nxmat>0 && nxmat<=maxmat){
    printf("Exclusion of certain region: YES [ ");
    for(i=0;i<nxmat;i++){
      printf("%d ",xmat[i]);
    }
      printf("]\n");
    }
  else {
    printf("exclusion of certain region: NO [*]\n");
  }
  printf("--------------------------------\n");

  printf("reading header...\n");
  /* read header information ASCII part */
  while(fgets(buffer,maxline,inf) != NULL && !matchfirstword(buffer,"<AppendedData")){
    if(matchfirstword(buffer,"<ImageData")){
      /* WholeExtent */
      getvalue(buffer,"WholeExtent","s",(int *)string);
      getfirstquote(string,stag);     
      sscanf(stag,"%d %d %d %d %d %d",&wx1,&wx2,&wy1,&wy2,&wz1,&wz2);     
      printf("WholeExtent: %d %d %d %d %d %d\n",wx1,wx2,wy1,wy2,wz1,wz2);
      
      /* Origin */
      getvalue(buffer,"Origin","s",(int *)string);
      getfirstquote(string,stag);           
      sscanf(stag,"%f %f %f",&ox,&oy,&oz);      
      printf("Origin: %f %f %f\n",ox,oy,oz);
      
      /* Spacing */
      getvalue(buffer,"Spacing","s",(int *)string);
      getfirstquote(string,stag);
      sscanf(stag,"%f %f %f",&dx,&dy,&dz);      
      printf("Spacing: %f %f %f\n",dx,dy,dz);
    }
    
    if(matchfirstword(buffer,"<Piece")){
      /* WholeExtent */
      getvalue(buffer,"Extent","s",(int *)string);
      getfirstquote(string,stag);     
      sscanf(stag,"%d %d %d %d %d %d",&px1,&px2,&py1,&py2,&pz1,&pz2);     
      printf("Piece Extent: %d %d %d %d %d %d\n",px1,px2,py1,py2,pz1,pz2);
    }
    
    if(matchfirstword(buffer,"<DataArray")){
      /* type */
      getvalue(buffer,"type","s",(int *)string);
      getfirstquote(string,stag);     
      strcpy(vtype,stag);     
      printf("type: %s\n",vtype);
      
      /* Name */
      getvalue(buffer,"Name","s",(int *)string);
      getfirstquote(string,stag);     
      strcpy(vname,stag);     
      printf("Name: %s\n",vname);
    }
  }
  /* Transition character */
  ch=fgetc(inf);
  printf("transition char: %c\n",ch);
  if(ch != '_'){
    printf("ERROR: wrong transition character found!\n");
    exit(-1);
  } 
        
  printf("header status: SUCCESS\n");
  printf("--------------------------------\n");
  
  printf("original model...\n");
  
  /* Number of grid points */
  nx=wx2-wx1+1;
  ny=wy2-wy1+1;
  nz=wz2-wz1+1;

  ndim=0;
  if(nx>1)ndim+=1;
  if(ny>1)ndim+=1;
  if(nz>1)ndim+=1;

  if(ndim==3){
	  nenode=8;
	  etype=12; /* VTK hexahedron */
  }
  else if(ndim==2){
	  nenode=4;
	  etype=9; /* VTK quadrilateral */
  }
  else if(ndim==1){
	  nenode=2;
	  etype=3; /* VTK line */
  }
  enode=malloc(nenode*sizeof(int));
  
  /* Number of element */
  nex=floor(nx/step);
  ney=floor(ny/step); 
  nez=floor(nz/step);
  if(nex==0)nex=1;
  if(ney==0)ney=1;
  if(nez==0)nez=1;
  nelmt=nex*ney*nez;
  
  /* Number of nodes */
  nnx=nex+1;
  nny=ney+1;
  nnz=nez+1;
  if(ndim==1){
    nny=nnz=1;
  }
  else if(ndim==2){
    nnz=1;
  }
  nnode=nnx*nny*nnz;  
  
  /* Memory allocation */
  elmt_node=malloc (nelmt * sizeof(int *));
  if(elmt_node == NULL){
    fprintf(stderr, "ERROR: out of memory\n");
    exit(-1);
  }

  for (i=0; i<nelmt; i++){
    elmt_node[i]=malloc (nenode * sizeof(int ));
    if(elmt_node[i] == NULL){
      fprintf(stderr, "ERROR: out of memory\n");
      exit(-1);
    }
  }

  nstat=malloc(nnode * sizeof(int));
  if(nstat == NULL){
    fprintf(stderr, "ERROR: out of memory\n");
    exit(-1);
  }
  estat=malloc(nelmt * sizeof(int));
  if(estat == NULL){
    fprintf(stderr, "ERROR: out of memory\n");
    exit(-1);
  }
  
  elmt_mat=malloc(nelmt * sizeof(int));
  if(elmt_mat == NULL){
    fprintf(stderr, "ERROR: out of memory\n");
    exit(-1);
  }
  
  /* Initialization */
  for(i=0;i<nnode;i++)nstat[i]=OFF;
  for(i=0;i<nelmt;i++)estat[i]=OFF;    
  
  /* New origin */
  
  /* always set origin at the bottom */
  /* Change origin from top to bottom if necessary */
  ox = ox-0.5*dx;
  if(xup==0){
    dx=fabs(dx);    
    ox=ox-nex*(step*dx);
  }
  if(ndim>1){
    oy = oy-0.5*dy;
    if(yup==0){
      dy=fabs(dy);    
      oy=oy-ney*(step*dy);
	  }
  }
  if(ndim>2){
    oz = oz-0.5*dz;
    if(zup==0){
      dz=fabs(dz);    
      oz=oz-nez*(step*dz);
	  }
  }

  /* new sampling nterval */
  dx=step*dx;
  dy=step*dy;
  dz=step*dz;
  
  removeExtension(argv[1],fonly); 
  sprintf(outfname,"%s_mesh.vtu",fonly);
  
  outf0=fopen(outfname,"w");

  /* Bytes and offsets for ParaView file */ 
  bytes[0] = (3*nnode)*sizeof(float);i
  /* Coordinates: in paraview both x, y, and z coordinates should be provided
     for all dimensions */
  bytes[1] = (nenode*nelmt)*sizeof(int);  /* Connectivity */
  bytes[2] = (nelmt)*sizeof(int);     /* Offsets */
  bytes[3] = (nelmt)*sizeof(int);     /* Types */
  bytes[4] = (nelmt)*sizeof(float);     /* Cell data */
  
  offset[0]=0; /* 1st offset */
  for (i=0; i<nvar; i++){   
    if(i<nvar-1)offset[i+1]=offset[i]+sizeof(int)+bytes[i];
    bytes[i]=bytes[i]+sizeof(int);
  } 

  /* Header for VTK XML .vtu file */
  fprintf(outf0,"<?xml version=\"1.0\"?>\n");
  fprintf(outf0,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",byte_order);
  fprintf(outf0,"<UnstructuredGrid>\n");
  fprintf(outf0,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",nnode,nelmt);
  fprintf(outf0,"<Points>\n");
  fprintf(outf0,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset[0]);
  fprintf(outf0,"</Points>\n");
  fprintf(outf0,"<PointData>\n");
  fprintf(outf0,"</PointData>\n");    
  fprintf(outf0,"<Cells>\n");
  fprintf(outf0,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\"/>\n",offset[1]);
  fprintf(outf0,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\"/>\n",offset[2]);
  fprintf(outf0,"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%d\"/>\n",offset[3]);
  fprintf(outf0,"</Cells>\n");
  fprintf(outf0,"<CellData>\n");
  fprintf(outf0,"<DataArray type=\"Float32\" Name=\"Velocity\" format=\"appended\" offset=\"%d\"/>\n",offset[4]);
  fprintf(outf0,"</CellData>\n");
  fprintf(outf0,"</Piece>\n");
  fprintf(outf0,"</UnstructuredGrid>\n");
  fprintf(outf0,"<AppendedData encoding=\"raw\">\n");
  fprintf(outf0,"_");

  /* Coordinates */
  outf1=fopen("nodes_coords_file","w");
  fprintf(outf1,"%d\n",nnode);  

  next=0;

  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  /* loop through all the nodes */
  inode=0;
  for(k=0;k<nnz;k++){
    tempz=oz+k*dz; if(fac!=1.0)tempz*=fac;
    for(j=0;j<nny;j++){
      tempy=oy+j*dy; if(fac!=1.0)tempy*=fac;
      for(i=0;i<nnx;i++){
        tempx=ox+i*dx; if(fac!=1.0)tempx*=fac;
        fwrite(&tempx,sizeof(float),1,outf0);
        fwrite(&tempy,sizeof(float),1,outf0);
        fwrite(&tempz,sizeof(float),1,outf0);
        if(ndim==3){			       
          fprintf(outf1,"%d %.6f %.6f %.6f\n",inode+1,tempx,tempy,tempz);
        }
        else if(ndim==2){       
          fprintf(outf1,"%.6f %.6f\n",tempx,tempy);
        }
        else if(ndim==1){       
          fprintf(outf1,"%.6f\n",tempx);
        }
        inode++;
      }
    }
  }
  fclose(outf1);
  printf("coordinates: SUCCESS\n"); 
  
  /* Connectivity */
  outf1=fopen("mesh_file","w");
  fprintf(outf1,"%d\n",nelmt);

  /* Absorbing surface files */
  if(ndim==3){
	  outf_xmin=fopen("absorbing_surface_file_xmin","w");
	  outf_xmax=fopen("absorbing_surface_file_xmax","w");
	  outf_ymin=fopen("absorbing_surface_file_ymin","w");
	  outf_ymax=fopen("absorbing_surface_file_ymax","w");
	  outf_zmin=fopen("absorbing_surface_file_bottom","w");
    /* for pyhaesalmi I will write both zmin and zmax in same file */
	  printf("WARNING: top surface is saved as absorbing boundary and appended to bottom surface!\n");
	  
	  fprintf(outf_xmin,"%d\n",ney*nez);
	  fprintf(outf_xmax,"%d\n",ney*nez);
	  fprintf(outf_ymin,"%d\n",nex*nez);
	  fprintf(outf_ymax,"%d\n",nex*nez);
	  fprintf(outf_zmin,"%d\n",2*(nex*ney)); 
  }
  else if(ndim==2){
	  outf_abs=fopen("absorbing_surface_file","w");
	  fprintf(outf_abs,"%d\n",2*(nex+ney));
	  printf("WARNING: all periphery is saved as absorbing boundary surface!\n");
  }
  
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  ielmt=0;
  for(k=0;k<nez;k++){
    for(j=0;j<ney;j++){
      for(i=0;i<nex;i++){
        if(ndim==3){
          /* This segment is only valid for 8-noded hexahedron */       
          enode[0]=k*nny*nnx+j*nnx+i;
          enode[1]=enode[0]+1;
	  
          enode[2]=enode[1]+nnx;
          enode[3]=enode[0]+nnx;
	  
          enode[4]=enode[0]+nnx*nny;
          enode[5]=enode[4]+1;
              
          enode[6]=enode[5]+nnx;
          enode[7]=enode[4]+nnx;
          /*---------------------------------------------------*/
        }
        else if(ndim==2){
          /* This segment is only valid for 4-noded quadrilateral */       
          enode[0]=j*nnx+i;
          enode[1]=enode[0]+1;
	  
          enode[2]=enode[1]+nnx;
          enode[3]=enode[0]+nnx;
          /*---------------------------------------------------*/
        }
        
        /* write absorbing boundary at xmin */
        if(ndim==3){
          if(i==0){
            fprintf(outf_xmin,"%d %d %d %d %d\n",ielmt+1,enode[3]+1,enode[0]+1,enode[4]+1,enode[7]+1);
          }
          /* write absorbing boundary at xmax */
          if(i==nex-1){
            fprintf(outf_xmax,"%d %d %d %d %d\n",ielmt+1,enode[2]+1,enode[1]+1,enode[5]+1,enode[6]+1);
          }
          /* write absorbing boundary at ymin */
          if(j==0){
            fprintf(outf_ymin,"%d %d %d %d %d\n",ielmt+1,enode[0]+1,enode[1]+1,enode[5]+1,enode[4]+1);
          }
          /* write absorbing boundary at ymin */
          if(j==ney-1){
            fprintf(outf_ymax,"%d %d %d %d %d\n",ielmt+1,enode[3]+1,enode[2]+1,enode[6]+1,enode[7]+1);
          }
          /* write absorbing boundary at zmin */
          if(k==0){
            fprintf(outf_zmin,"%d %d %d %d %d\n",ielmt+1,enode[0]+1,enode[1]+1,enode[2]+1,enode[3]+1);
          }
          /* write absorbing boundary at zmax */
          if(k==nez-1){
            fprintf(outf_zmin,"%d %d %d %d %d\n",ielmt+1,enode[4]+1,enode[5]+1,enode[6]+1,enode[7]+1);
          }
        }
        else if(ndim==2){
          if(i==0){
            fprintf(outf_abs,"%d %d %d\n",ielmt+1,enode[0]+1,enode[3]+1);
          }
          /* write absorbing boundary at xmax */
          if(i==nex-1){
            fprintf(outf_abs,"%d %d %d\n",ielmt+1,enode[1]+1,enode[2]+1);
          }
          /* write absorbing boundary at ymin */
          if(j==0){
            fprintf(outf_abs,"%d %d %d\n",ielmt+1,enode[0]+1,enode[1]+1);
          }
          /* write absorbing boundary at ymin */
          if(j==ney-1){
            fprintf(outf_abs,"%d %d %d\n",ielmt+1,enode[3]+1,enode[2]+1);
          }
        }

        fprintf(outf1,"%d ",ielmt+1);
        for(ienode=0;ienode<nenode;ienode++){
          elmt_node[ielmt][ienode]=enode[ienode];
          fwrite(&enode[ienode],sizeof(int),1,outf0);
          fprintf(outf1,"%d ",enode[ienode]+1);
          if(ienode==nenode-1)fprintf(outf1,"\n");        
        }
        ielmt++;
      }
    }
  }
  //exit(-1);
  if(ndim==3){
	  fclose(outf_xmin);
	  fclose(outf_xmax);
	  fclose(outf_ymin);
	  fclose(outf_ymax);
	  fclose(outf_zmin);
  }
  else if(ndim==2){
	  fclose(outf_abs);
  }
  fclose(outf1);
  printf("connectivity: SUCCESS\n");
  
  /* Offsets */
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  temp=0;
  for(i=0;i<nelmt;i++){
    temp+=nenode;
    fwrite(&temp,sizeof(int),1,outf0);
  }
  printf("offsets: SUCCESS\n");
  
  /* Types */
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;  
  for(i=0;i<nelmt;i++){
    fwrite(&etype,sizeof(int),1,outf0);
  }
  printf("types: SUCCESS\n");
  
  /* Cell data */ 
  outf1=fopen("materials_val","wb");
  fread(&dumi,sizeof(int),1,inf);
  off0=ftell(inf);  /* from this offset gridded data starts */

  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  
  for(imat=0;imat<maxmat;imat++)mat[imat]=0.0;
  
  if(zup==0){
    kdatum=nz-1;
    ksign=-1;
  }
  nmat=0;
  for(ki=0;ki<nz;ki+=step){
    k=kdatum+ksign*ki;
    for(ji=0;ji<ny;ji+=step){
      j=jdatum+jsign*ji;
      for(ii=0;ii<nx;ii+=step){
        i=idatum+isign*ii;
        ioff=sizeof(float)*(k*ny*nx+j*nx+i);/*printf("%d\n",ioff); exit(-1);*/
        fseek(inf,off0+ioff,SEEK_SET);
        fread(&tempx,sizeof(float),1,inf);          
        fwrite(&tempx,sizeof(float),1,outf0);
        fwrite(&tempx,sizeof(float),1,outf1);

        /* count and designate material regions */        
        if(nmat==0){/* Initialize */
          mat[nmat]=tempx;printf("material property: %.6f\n",tempx);
          nmat++;
        }       
        else{
          match=OFF;
          for(imat=0;imat<nmat;imat++){
            if(fabs(tempx-mat[imat])<eps){
              /* current value matched with stored value */
              match=ON;
              break;
            }
          }

          if(match == OFF){
            mat[nmat]=tempx;printf("material property: %.6f\n",tempx);
            nmat++;
            if(nmat>maxmat){
              printf("ERROR: material number %d > bound %d!\n",nmat,maxmat);
              exit(-1);
            }
          }
        } /* if(nmat */       
      } /* for(i */
    } /* for(j */
  } /* for(ki */
  fclose(inf);
  fclose(outf1);

  /* sort the material designation in ascending order and diplay */ 
  qsort(mat,nmat,sizeof(float),comp_float);
  printf("material properties: num=%d, min=%f, max=%f\n",nmat,mat[0],mat[nmat-1]);

  /* new material ID after exclusion */
  new_matid=malloc(nmat*sizeof(int));
  for(i=0;i<nmat;i++){
    new_matid[i]=i+1;		
  }
  for(i=0;i<nxmat;i++){
    if(xmat[i]<1 || xmat[i]>nmat){
      printf("ERROR: exclusion material ID: %d is outside the range: [%d %d]\n",xmat[i],1,nmat);
      exit(-1);
    }		
    new_matid[xmat[i]-1]=OFF;
  }
  matid=0;
  for(i=0;i<nmat;i++){		
    if(new_matid[i]>OFF){
      matid++;
      new_matid[i]=matid;
    }
  }
  if(matid<=0){
    printf("WARNING: all material regions cannot be excluded! Nothing will be excluded!\n");
    nxmat=0;
  }
  
  inf=fopen("materials_val","rb");
  outf1=fopen("materials_file","w");
  ielmt=0;new_nelmt=0;
  
  ielmt=0;  
  for(k=0;k<nez;k++){
    for(j=0;j<ney;j++){
      for(i=0;i<nex;i++){
        /*ioff=sizeof(float)*(k*ney*nx+j*nx+i);
        fseek(inf,off0+ioff,SEEK_SET);*/
        fread(&tempx,sizeof(float),1,inf);
        emat_id=0;
        for(imat=0;imat<nmat;imat++){
          if(fabs(tempx-mat[imat]) <= eps){
            elmt_mat[ielmt]=imat+1;
            fprintf(outf1,"%d %d\n",ielmt+1,imat+1);
            emat_id=1;
            break;
          }
        }
        if(emat_id==0){
          printf("ERROR: material mismatched [%d %f]!\n",ielmt,tempx);
          exit(-1);
        }
        
        if(nxmat>0){ /* Exclusion of material */
          if(new_matid[elmt_mat[ielmt]-1]>OFF){
            new_nelmt++;
            estat[ielmt]=ON;
            for(inode=0;inode<nenode;inode++)nstat[elmt_node[ielmt][inode]]=ON;
          }
          /*if(tempx != mat[xmat[1]-1]){
              new_nelmt++;
              estat[ielmt]=ON;
              for(inode=0;inode<nenode;inode++)nstat[elmt_node[ielmt][inode]]=ON;
          }*/
        }
        ielmt++;
      }
    }
  }
  free(elmt_node);
  fclose(outf1);    
  fclose(inf);
  printf("cell data: SUCCESS\n"); 
  
  fprintf(outf0,"\n");
  fprintf(outf0,"</AppendedData>\n");
  fprintf(outf0,"</VTKFile>\n");
  fclose(outf0);
  
  printf("nodes: %d\n",nnode);
  printf("elements: %d\n",nelmt);
  printf("status: SUCCESS\n");
  printf("--------------------------------\n");
  
  if(nxmat<=0 || nxmat>maxmat)return(0);

  /* processing to exclude air */

  printf("exclude certain region...\n");
  
  /* count new nodes */
  new_nnode=0;
  for(i=0;i<nnode;i++){
    if(nstat[i]==ON)new_nnode++;
  }

  sprintf(outfname,"%s_mesh_xmat.vtu",fonly);
  
  outf0=fopen(outfname,"wb");
  
  /* bytes and offsets for ParaView file */
  bytes[0] = (3*new_nnode)*sizeof(float);  /* Coordinates: in paraview both x, y, and z coordinates should be provided for all dimensions*/
  bytes[1] = (nenode*new_nelmt)*sizeof(int);  /* Connectivity */
  bytes[2] = (new_nelmt)*sizeof(int);     /* Offsets */
  bytes[3] = (new_nelmt)*sizeof(int);     /* Types */
  bytes[4] = (new_nelmt)*sizeof(float);     /* Cell data */
  
  offset[0]=0; /* 1st offset */
  for (i=0; i<nvar; i++){   
    if(i<nvar-1)offset[i+1]=offset[i]+sizeof(int)+bytes[i];
    bytes[i]=bytes[i]+sizeof(int);
  }
  
  /* header for VTK XML .vtu file */
  fprintf(outf0,"<?xml version=\"1.0\"?>\n");
  fprintf(outf0,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"%s\">\n",byte_order);
  fprintf(outf0,"<UnstructuredGrid>\n");
  fprintf(outf0,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",new_nnode,new_nelmt);
  fprintf(outf0,"<Points>\n");
  fprintf(outf0,"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset[0]);
  fprintf(outf0,"</Points>\n");
  fprintf(outf0,"<PointData>\n");
  fprintf(outf0,"</PointData>\n");    
  fprintf(outf0,"<Cells>\n");
  fprintf(outf0,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%d\"/>\n",offset[1]);
  fprintf(outf0,"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%d\"/>\n",offset[2]);
  fprintf(outf0,"<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%d\"/>\n",offset[3]);
  fprintf(outf0,"</Cells>\n");
  fprintf(outf0,"<CellData>\n");
  fprintf(outf0,"<DataArray type=\"Float32\" Name=\"Velocity\" format=\"appended\" offset=\"%d\"/>\n",offset[4]);
  fprintf(outf0,"</CellData>\n");
  fprintf(outf0,"</Piece>\n");
  fprintf(outf0,"</UnstructuredGrid>\n");
  fprintf(outf0,"<AppendedData encoding=\"raw\">\n");
  fprintf(outf0,"_");


  nmir = malloc (nnode * sizeof(int));
  if(nmir == NULL){
    fprintf(stderr, "ERROR: out of memory!\n");
    exit(-1);
  }

  /* coordinates and mirror */
  inf=fopen("nodes_coords_file","r");
  fscanf(inf,"%d\n",&dumi);
  
  outf1=fopen("nodes_coords_file_xmat","w");
  fprintf(outf1,"%d\n",new_nnode);
  
  next=0;
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  inode=0;
  for(i=0;i<nnode;i++){
    if(ndim==3){
	  fscanf(inf,"%d %f %f %f\n",&dumi,&tempx,&tempy,&tempz);
	}
	else if(ndim==2){
	  fscanf(inf,"%f %f\n",&tempx,&tempy);
	}
    else if(ndim==1){
	  fscanf(inf,"%f\n",&dumi,&tempx);
	}
    if(nstat[i]==ON){
		fwrite(&tempx,sizeof(float),1,outf0);
		fwrite(&tempy,sizeof(float),1,outf0);
		fwrite(&tempz,sizeof(float),1,outf0);
		if(ndim==3){		  
		  fprintf(outf1,"%d %.6f %.6f %.6f\n",inode+1,tempx,tempy,tempz);
		}
		else if(ndim==2){
		  fprintf(outf1,"%.6f %.6f\n",tempx,tempy);
		}
		else if(ndim==1){
		  fprintf(outf1,"%.6f\n",tempx);
		}
		nmir[i]=inode;
		inode++;
    }
  }
  free(nstat);
  fclose(outf1);
  printf("coordinates and mirror: SUCCESS\n");
  
  /* connectivity */
  inf=fopen("mesh_file","r");
  fscanf(inf,"%d\n",&dumi);
  
  outf1=fopen("mesh_file_xmat","w");
  fprintf(outf1,"%d\n",new_nelmt);
  
  /* absorbing surface files */
  if(ndim==3){
	  outf_xmin=fopen("absorbing_surface_file_xmin_xmat","w");
	  outf_xmax=fopen("absorbing_surface_file_xmax_xmat","w");
	  outf_ymin=fopen("absorbing_surface_file_ymin_xmat","w");
	  outf_ymax=fopen("absorbing_surface_file_ymax_xmat","w");
	  outf_zmin=fopen("absorbing_surface_file_bottom_xmat","w"); /* for pyhaesalmi I will write both zmin and zmax in same file */
	  printf("WARNING: top surface is saved as absorbing boundary and appended to bottom surface!\n");  
	  
	  /* this must be modified if the boundaries are not intact after the exclusion of certain region/s */
	  fprintf(outf_xmin,"%d\n",ney*nez);
	  fprintf(outf_xmax,"%d\n",ney*nez);
	  fprintf(outf_ymin,"%d\n",nex*nez);
	  fprintf(outf_ymax,"%d\n",nex*nez);
	  fprintf(outf_zmin,"%d\n",2*(nex*ney));
	  printf("WARNING: boundaries are assumed to be intact!\n");
  }
  else{
	  outf_abs=fopen("absorbing_surface_file_xmat","w");  
    printf("WARNING: all periphery is saved as absorbing boundary surface!\n");

	  /* this must be modified if the boundaries are not intact after the 
    exclusion of certain region/s */
	  fprintf(outf_abs,"%d\n",2*(nex+ney));
	  printf("WARNING: boundaries are assumed to be intact!\n");
  }
  
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  
  ielmt=0;
  inew_elmt=0;
  for(k=0;k<nez;k++){
    for(j=0;j<ney;j++){
      for(i=0;i<nex;i++){
        /*for(i=0;i<nelmt;i++){*/
        /* node number just follows represents the numbering starting from 1 */
        if(ndim==3){
          fscanf(inf,"%d %d %d %d %d %d %d %d %d\n",&dumi,&enode[0],&enode[1],\
          &enode[2],&enode[3],&enode[4],&enode[5],&enode[6],&enode[7]);
        }
        else if(ndim==2){
          fscanf(inf,"%d %d %d %d %d\n",&dumi,&enode[0],&enode[1],&enode[2],  \
          &enode[3]);
        }
        else if(ndim==1){
          fscanf(inf,"%d %d %d\n",&dumi,&enode[0],&enode[1]);
        }
        if(estat[ielmt]==ON){
          /* nmir was determined based on the numbering starting from 0 but 
          SPECFEM uses node number starting from 1 */       
          /* write connectivity considering mirror points to new nodes */
          fprintf(outf1,"%d ",inew_elmt+1);
          for(ienode=0;ienode<nenode;ienode++){
            enode[ienode]=nmir[enode[ienode]-1];
            /* nmir was determined based on the numbering starting from 0 */
            fwrite(&enode[ienode],sizeof(int),1,outf0);
            fprintf(outf1,"%d ",enode[ienode]+1);
            /* SPECFEM uses node number starting from 1 not 0 */
            if(ienode==nenode-1)fprintf(outf1,"\n");
          }
          if(ndim==3){
            /* write absorbing boundary at xmin */
            if(i==0){
              fprintf(outf_xmin,"%d %d %d %d %d\n",inew_elmt+1,enode[3]+1,    \
              enode[0]+1,enode[4]+1,enode[7]+1);
            }
            /* write absorbing boundary at xmax */
            if(i==nex-1){
              fprintf(outf_xmax,"%d %d %d %d %d\n",inew_elmt+1,enode[2]+1,    \
              enode[1]+1,enode[5]+1,enode[6]+1);
            }
            /* write absorbing boundary at ymin */
            if(j==0){
              fprintf(outf_ymin,"%d %d %d %d %d\n",inew_elmt+1,enode[0]+1,    \
              enode[1]+1,enode[5]+1,enode[4]+1);
            }
            /* write absorbing boundary at ymin */
            if(j==ney-1){
              fprintf(outf_ymax,"%d %d %d %d %d\n",inew_elmt+1,enode[3]+1,    \
              enode[2]+1,enode[6]+1,enode[7]+1);
            }
            /* write absorbing boundary at zmin */
            if(k==0){
              fprintf(outf_zmin,"%d %d %d %d %d\n",inew_elmt+1,enode[0]+1,    \
              enode[1]+1,enode[2]+1,enode[3]+1);
            }
            /* write absorbing boundary at zmax */
            if(k==nez-1){
              fprintf(outf_zmin,"%d %d %d %d %d\n",inew_elmt+1,enode[4]+1,    \
              enode[5]+1,enode[6]+1,enode[7]+1);
            }
          }
          else if(ndim==2){
            /* write absorbing boundary at xmin */
            if(i==0){
              fprintf(outf_abs,"%d %d %d\n",inew_elmt+1,enode[0]+1,enode[3]+1);
            }
            /* write absorbing boundary at xmax */
            if(i==nex-1){
              fprintf(outf_abs,"%d %d %d\n",inew_elmt+1,enode[1]+1,enode[2]+1);
            }
            /* write absorbing boundary at ymin */
            if(j==0){
              fprintf(outf_abs,"%d %d %d\n",inew_elmt+1,enode[0]+1,enode[1]+1);
            }
            /* write absorbing boundary at ymin */
            if(j==ney-1){
              fprintf(outf_abs,"%d %d %d\n",inew_elmt+1,enode[3]+1,enode[2]+1);
            }
          }
          inew_elmt++;        
        }
        ielmt++;        
      }
    }   
  }
  if(inew_elmt!=new_nelmt){
    printf("ERROR: number of elements mismatch for exclusion!\n");
    exit(-1);
  }
  free(nmir);
  fclose(inf);
  if(ndim==3){
	  fclose(outf_xmin);
	  fclose(outf_xmax);
	  fclose(outf_ymin);
	  fclose(outf_ymax);
	  fclose(outf_zmin);
  }
  else if(ndim==2){
	  fclose(outf_abs);
  }
  fclose(outf1);
  printf("connectivity: SUCCESS\n");

  /* Offsets */
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  temp=0;
  for(i=0;i<new_nelmt;i++){
    temp+=nenode;
    fwrite(&temp,sizeof(int),1,outf0);
  }
  printf("offsets: SUCCESS\n");
  
  /* Types */
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;
  for(i=0;i<new_nelmt;i++){
    fwrite(&etype,sizeof(int),1,outf0);
  }
  printf("types: SUCCESS\n");
  
  /* Cell data */
  inf=fopen("materials_val","rb");
  
  fwrite(&bytes[next],sizeof(int),1,outf0);next++;

  outf1=fopen("materials_file_xmat","w");
  //ielmt=0;
  for(i=0;i<nelmt;i++){
    fread(&tempx,sizeof(float),1,inf);        
    if(estat[i]==ON){
      fwrite(&tempx,sizeof(float),1,outf0);     
      fprintf(outf1,"%d\n",new_matid[elmt_mat[i]]);
        /*if(xmat[1]<elmt_mat[i]){
            fprintf(outf1,"%d\n",elmt_mat[i]-1);
        }
        else if(xmat[1]>elmt_mat[i]){
            fprintf(outf1,"%d\n",elmt_mat[i]);
        }
        else{
            printf("ERROR: materials ID mismatched!\n");
            exit(-1);
        }*/
      //ielmt++;
    }
  }
  free(xmat);
  free(new_matid);
  free(elmt_mat);
  free(estat);
  fclose(inf);
  fclose(outf1);

  if(remove("materials_val"))printf("WARNING: \"mesh.tmp\" cannot be deleted!\n");
  
  printf("cell data: SUCCESS\n");
    
  fprintf(outf0,"\n");
  fprintf(outf0,"</AppendedData>\n");
  fprintf(outf0,"</VTKFile>\n");
  fclose(outf0);
  
  printf("new nodes: %d\n",new_nnode);
  printf("new elements: %d\n",new_nelmt);
  printf("status: SUCCESS\n");
  printf("--------------------------------\n");     

  return(0);
}
/*===========================================================================*/
