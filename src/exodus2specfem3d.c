/** @file exodus2specfem3d.c
*  @brief Converts ASCII exodus file to SPECFEM3D files.
*
*  This program converts the Binary (provided that ncdump command exists,
*  type ncdump to check whether ncdump command exists.) or ASCII exodus file
*  exported from the TRELIS/CUBIT to several mesh files 
*  required by the SPECFEM3D Cartesian package. The ncdump commad is a part of
*  NetCFD library which is generally installed already in LINUX, 
*  which can be downloaded freely from
*  https://www.unidata.ucar.edu/software/netcdf/
*
*  @author Hom Nath Gharti (hgharti_AT_princeton_DOT_edu)
*
* ## Dependencies:
*  stringmanip.c: string manipulation routines
*
* ## Compile:
*  gcc exodus2specfem3d.c -o trelis2specfem3d
*
* ## Usage: 
*  exodus2specfem3d <inputfile> [OPTIONS]
*  Example: exodus2specfem3d tunnel.e -bin=1
*  or
*  exodus2specfem3d tunnel.txt
*
* ## Options:
*
* - -fac: use this option to multiply coordinates with some factor. this is 
*        important for unit 
*        conversion, e.g., to convert m to km use -fac=0.001 [DEFAULT 1]
*
* - -bin: use this option if you want to convert exodus binary directly, provided
*        that the command "ncdump" is in the path. The command "ncdump" is a 
*        part of netCDF library that can be downloaded freely from 
*        http://www.unidata.ucar.edu/downloads/netcdf/index.jsp.
*        use -bin=1 for binary or -bin=0 for ascii file. [DEFAULT 0]
*
* - -norm: use this option to check the normal of the faces in order to make sure
*        that the surface nodes are in the right order. use -norm=1 for checking
*        or -norm=0 for no checking [DEFAULT 0]. normally this is not necessary.
*
* - -head: use this option to attach head of input file to output file names. 
*        use -head=1 to attach header or -head=0 not to attach  [DEFAULT 0]
*
* - -tomo: use this option for tomography model. Since tomography model uses
*        negative identifiers, this option will write negative block IDs. 
*        use -tomo=1 to make negative block IDs or -tomo=0 not to make
*        [DEFAULT 0]
*
* ##Issues:
* - - This does not work with older verion of CUBIT. For the older version use
*    exodusold2specfem3d.c.
*
* # Basic steps starting from the TRELIS:
*
* ### step 1: prepare mesh in TRELIS/CUBIT
* - define material regions using "Blocks"
*   For example:
*   block 1 add volume 1
*   block 2 add volume 2 3
* 
*   will assign material region 1 to volume 1 and material region 2 to volumes 2
*   and 3. These material regions will be used to define material properties in
*   "nummaterial_velocity_file". this program will NOT generate 
*   "nummaterial_velocity_file". the file "nummaerial_veolicty_file" must be
*   created to run SPECFEM3D!
* 
* - define surface boundary conditions using "Sidesets"
*   
*   For example:
*   sideset 1 add surface 1
*   sideset 1 name 'free_or_absorbing_surface_file_zmax'
* 
*   will define a free or absorbing surface boundary condition on surface 1 which
*   lies at the top of the volume (zmax). similary,  
*   sideset 2 add surface 3
*   sideset 2 name 'absorbing_surface_file_bottom'
*   
*   will define absorbing boundary condition on the surface 3 which lies at the
*   bottom of the volume (zmin).
*   Note: All the above commands can also be executed using TRELIS/CUBIT GUI.
*   "sideset 1 name 'free_or_absorbing_surface_file_zmax'" is equivalent to
*   clicking sideset 1 and renaming.
*
* ### step2: export mesh file as exodus file say "tunnel.e" (use 3D option)
*
* ### step3: convert "tunnel.e" to SPECFEM3D files
*  exodus2specfem3d tunnel.e -bin=1
*
*There will be several output files:
*
* - nodes_coords_file : coordinates file => total number of nodes followed by 
*  nodal coordinate ? (? -> x, y, z)
*
* - mesh_file : element file => total number of elements followed by connectivity
*  list
*
* - materials_file : material file => total number of elements followed by 
*  material IDs
*
* - surface_file* : sourface boundary condition files => total number of elements
*  followed by element ID and surface nodes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

#define OFF 0
#define ON 1
#define MAXBULK 100000
#define MAXLINE 100

/* auxiliary routines */
void removeExtension(char *, char *);
int get_int(int *, char *, char *);
int look_int(int *, char *, char *);
int getfirstquote(char *, char *);
int check_normal(double [3][4],double [3]);


/* main routine */
int main(int argc,char **argv){
int i,ielmt,iface,itmp,j,k;
int bulksize; 
int ndim; /* geometry dimension */ 
int nnode,nelmt; /* number of nodes, number of elements */
int nblk,nss; /* number of blocks, number of side sets */
int elmt_count,node_count; /* element, node count */
int node_countx,node_county,node_countz; /* element, node count */
int blk_count,ss_count; /* block, node set count */
int dim_stat,ss_rstat,ss_stat,con_stat,coord_stat,mat_stat; /* status */
int coordx_stat,coordy_stat,coordz_stat; /* status */
int *blk_nelmt,*blk_nenod; /* number of elements, number of nodes per element in each node */
int *ss_nside; /* number of nodes in each node set */

double fac,dtmp; /* multiplication factor for coordinates, temporary double */
double **coord;
int **elmt_node;
char *bulk,line[100],token[62],stag[62];
char dumc[250];
char **coord_name; /* coordinates name */
char fonly[62],infname[62],outfname[62],outhead[62];
char **ss_name; /* node set names */

/* change this line to look for other BCs in side sets */
double normal[3],p[3][4];
/* int fnod[6][4]={{1,2,6,5},{2,3,7,6},{3,4,8,7},{4,1,5,8},{2,1,4,3},{5,6,7,8}}; */ /* face nodes */
int fnod[6][4]={{1,2,6,5},{2,3,7,6},{3,4,8,7},{1,5,8,4},{1,4,3,2},{5,6,7,8}}; /* face nodes */
int igfnod,gfnod[4]; /* global face nodes */
int nface_change;

int count_sselmt,count_ssside,sumss_nelmt,*ss_elmt,*ss_side; /* total number of element in absorbing BC and free surface */

int isbin; /* test if binary */
int ishead; /* option to add innput file header in the output file names */
int istomo; /* option to convert block IDs to negative for tomography models */
int ndir,isnorm; /* normal direction, test if normal has to checked */

FILE *inf,*outf_mat,*outf_con,*outf_coord,*outf_side;

/* default factor and binary switch*/    
fac=1.0; isbin=OFF; isnorm=OFF; ishead=OFF; istomo=OFF; ndir=1;

if(argc<2){
  fprintf(stderr,"ERROR: input file not entered!\n");
  exit(-1);
}

/* scan command */
if(argc>2){
  for(i=2;i<argc;i++){
    if(look_double(&dtmp,"-fac=",argv[i])==0){
    fac=dtmp;     
    continue;
  }
  else if(look_int(&itmp,"-bin=",argv[i])==0){
    isbin=itmp;
    continue;
  }else if(look_int(&itmp,"-norm=",argv[i])==0){
    isnorm=itmp;
    continue;
  }else if(look_int(&itmp,"-head=",argv[i])==0){
    ishead=itmp;
    continue;
  }else if(look_int(&itmp,"-tomo=",argv[i])==0){
    istomo=itmp;
    continue;
  }else{
    printf("ERROR: unrecognized option \"%s\"",argv[i]);
    exit(-1);
  }
  }
}

printf("input file: %s\n",argv[1]);
printf("fac: %f\n",fac);
printf("isbin: %d\n",isbin);
printf("isnorm: %d\n",isnorm);
printf("ishead: %d\n",ishead);
printf("istomo: %d\n",istomo);

printf("--------------------------------\n");
/* default input file name is argv[1]*/
strcpy(infname,argv[1]);
removeExtension(argv[1],fonly);

if (isbin){
  printf("converting binary to ascii...");
  /* set input file name */
  sprintf(infname,"%s.txt",fonly);

  /* convert binary netCDF file to ascii file */
  sprintf(dumc,"ncdump %s > %s.txt",argv[1],fonly);  
  printf("conversion command: \"%s\" !\n",dumc);
  if (system(dumc)!=0){
  printf("ERROR: command \"%s\" cannot be executed! use -bin=0 or no option for ascii input file! \n",dumc);
  exit(-1);
  }
  printf("complete!\n");
}

/* open input file */
inf=fopen(infname,"r");
if(inf==NULL){
  fprintf(stderr,"ERROR: file \"%s\" not found!",argv[1]);
  exit(-1);
}
/*printf("--------------------------------\n");*/

/* set out header */
if(ishead){
  sprintf(outhead,"%s_",fonly);
}else{
  /* outhead[0]='\0'; */
  sprintf(outhead,"%s","");
}
/* if the size below is not large enough segmentation fault might occur
   elsewhere such as in malloc */
bulk=malloc(MAXBULK); /* bulk string */
      
/* initialize some variables to 0 */
ndim=0; nblk=0; nnode=0; nelmt=0; nss=0;

/* intialize count to 0 */
blk_count=0; ss_count=0; node_count=0; elmt_count=0;  
node_countx=0; node_county=0; node_countz=0;

/* set default status to OFF */
dim_stat=OFF; ss_rstat=OFF; ss_stat=OFF; con_stat=OFF; coord_stat=OFF;
coordx_stat=OFF; coordx_stat=OFF; coordx_stat=OFF;

fscanf(inf,"%s",token);
if(strcmp(token,"netcdf")!=0){
  printf("ERROR: invalid exodus file or wrong -bin value!\n");
  printf("HINT: try correct value for -bin option or use valid exodus file!\n");
  exit(-1);
}
bulksize=0;
count_sselmt=0; count_ssside=0;
while(!feof(inf)){
  fscanf(inf,"%s",token);
  /* read dimensions */
  if(dim_stat!=ON && strcmp(token,"dimensions:")==0){
    printf("reading dimensions...");
    while(strstr(fgets(line,MAXLINE,inf),"variables:") == NULL){    
      bulksize+=MAXLINE;
      if(bulksize>MAXBULK){
        printf("ERROR: not enough bulk size! Set MAXBULK > %d!\n",bulksize);
        exit(-1);
      }
      strncat(bulk,line,strcspn(line,";")+1);              
    }    
    get_int(&ndim,"num_dim =",bulk);
    if(ndim>0){   
      /* allocate memory */
      coord_name=malloc(ndim*sizeof(char *));
      for(i=0;i<ndim;i++){      
      coord_name[i]=malloc(62*sizeof(char)); /* each name has maximum of 62 characters */   
      }   
    }else{
      printf("ERROR: illegal value of dimension!\n");
      exit(-1);
    }
    get_int(&nnode,"num_nodes =",bulk);
    get_int(&nelmt,"num_elem =",bulk);
    get_int(&nblk,"num_el_blk =",bulk);
   
    /* allocate memory */
    blk_nelmt=malloc(nblk*sizeof(int));
    blk_nenod=malloc(nblk*sizeof(int));
    elmt_node=malloc(8*sizeof(int *));
    if(elmt_node==NULL){
      printf("ERROR: not enough memory!");
      exit(-1);
    }
    for(i=0;i<8;i++){
      elmt_node[i]=malloc(nelmt*sizeof(int)); /* need to make general */
      if(elmt_node[i]==NULL){
        printf("ERROR: not enough memory!");
        exit(-1);
      }
      /*printf("%d %d\n",nelmt,i);*/
    }
    coord=malloc(3*sizeof(double *));
    for(i=0;i<3;i++){
      coord[i]=malloc(nnode*sizeof(double));
      if(coord[i]==NULL){
        printf("ERROR: not enough memory!");
        exit(-1);
      }
    }

    /* sideset information */
    if (look_int(&nss,"num_side_sets =",bulk)!=0){
      nss=0;
    }else{
      /* allocate memory */
      ss_name=malloc(nss*sizeof(char *));
      for(i=0;i<nss;i++){
        ss_name[i]=malloc(62*sizeof(char)); /* each name has maximum of 62 characters */
      }
      ss_nside=malloc(nss*sizeof(int));
    }
    
    if(nss>0){ /* This segment has a significance only if nss has legitimate value */
      for(i=0;i<nss;i++){
        sprintf(stag,"num_side_ss%d =",i+1);
        get_int(&ss_nside[i],stag,bulk);          
      }
    }

    /* block information */        
    for(i=0;i<nblk;i++){
      sprintf(stag,"num_el_in_blk%d =",i+1);
      get_int(&blk_nelmt[i],stag,bulk); 
      
      sprintf(stag,"num_nod_per_el%d =",i+1);
      get_int(&blk_nenod[i],stag,bulk);
    }    
      
    dim_stat=ON;
    free(bulk);
    printf("complete!\n");
    printf(" geometry dimension: %d\n",ndim);    
    printf(" number of nodes: %d\n",nnode);
    printf(" number of elements: %d\n",nelmt);
    printf(" number of blocks: %d\n",nblk);
    printf(" number of sidesets: %d\n",nss);
    continue;
  }  
  
  /* Coordinates */
  if(strcmp(token,"coordx")==0){
    printf("reading x coordinates...");
    fscanf(inf,"%s",dumc);
    for(j=0;j<nnode;j++){
      fscanf(inf,"%lf,",&dtmp); /* read comma separated data */ 
      coord[0][j]=fac*dtmp;
      node_countx++;
    }
    coordx_stat=ON;
    printf("complete!\n");
    continue;
  }  
  if(strcmp(token,"coordy")==0){
    printf("reading y coordinates...");
    fscanf(inf,"%s",dumc);
    for(j=0;j<nnode;j++){
      fscanf(inf,"%lf,",&dtmp); /* read comma separated data */ 
      coord[1][j]=fac*dtmp;
      node_county++;
    }
    coordy_stat=ON;
    printf("complete!\n");
    continue;
  }  
  if(strcmp(token,"coordz")==0){
    printf("reading z coordinates...");
    fscanf(inf,"%s",dumc);
    for(j=0;j<nnode;j++){
      fscanf(inf,"%lf,",&dtmp); /* read comma separated data */ 
      coord[2][j]=fac*dtmp;
      node_countz++;
    }
    coordz_stat=ON;
    printf("complete!\n");
    continue;
  }
  
  /* read and store side boundary conditions */
  if(strcmp(token,"ss_names")==0){
    /* printf("saving side BCs..."); */
    fscanf(inf,"%s",dumc); /* = */      
    for (i=0; i<nss; i++){  
      fscanf(inf,"%s",dumc);   
      getfirstquote(dumc,ss_name[i]);
      if(!strcmp(ss_name[i],"")){ 
        sprintf(ss_name[i],"unnamed_sideset%d",i);
      }
    }
    /* total number of sideset elements */
    sumss_nelmt=0;
    for(i=0;i<nss;i++)sumss_nelmt+=ss_nside[i]; 

    ss_elmt=malloc(sumss_nelmt*sizeof(int));
    ss_side=malloc(sumss_nelmt*sizeof(int)); 
    continue;
  }

  /* read coordinate names */
  if(strcmp(token,"coor_names")==0){    
    fscanf(inf,"%s",dumc); /* = */      
    for (i=0; i<ndim; i++){
      fscanf(inf,"%s",dumc);    
      getfirstquote(dumc,coord_name[i]);    
    } 
    continue;
  }  
  
  if(ss_stat!=ON){
    for(i=0;i<nss;i++){
      sprintf(stag,"elem_ss%d",i+1);
      if(strcmp(token,stag)==0){
        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<ss_nside[i]; j++){
          fscanf(inf,"%d,",&ss_elmt[count_sselmt]); /* read comma separated data */
          count_sselmt++; 
        }
        
        fscanf(inf,"%s",dumc); /* ; */
        fscanf(inf,"%s",token);
        sprintf(stag,"side_ss%d",i+1);
        if(strcmp(token,stag)==0){
          fscanf(inf,"%s",dumc); /* = */
          for(j=0;j<ss_nside[i]; j++){
            fscanf(inf,"%d,",&ss_side[count_ssside]); /* read comma separated data */
            count_ssside++;
          }
        } 
        continue;
      }          
    }
  }
   
  /* Connectivity */
  if(nblk>0 && con_stat!=ON){   
   
    /* write connectivity and material id */  
    for(i=0;i<nblk;i++){
      sprintf(stag,"connect%d",i+1);
      if(strcmp(token,stag)==0){
        blk_count++;
        
        /* open connectivity and material files */
        if(blk_count==1){          
          printf("writing connectivity and materials..."); 
          sprintf(outfname,"%smesh_file",outhead);
          outf_con=fopen(outfname,"w");
          fprintf(outf_con,"%d\n",nelmt);
          
          sprintf(outfname,"%smaterials_file",outhead);
          outf_mat=fopen(outfname,"w");
          /* fprintf(outf_mat,"%d\n",nelmt); */
        }
        
        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<blk_nelmt[i];j++){
          fprintf(outf_con,"%d ",elmt_count+1);
          for(k=0;k<blk_nenod[i];k++){
            fscanf(inf,"%d,",&itmp);
            fprintf(outf_con,"%d ",itmp);
            elmt_node[k][elmt_count]=itmp;
          }
          fprintf(outf_con,"\n");/* new line */
          if(istomo){
            fprintf(outf_mat,"%d %d\n",elmt_count+1,-(i+1));
          }else{
            fprintf(outf_mat,"%d %d\n",elmt_count+1,i+1);
          }
          elmt_count++;
        }        
        
        if(blk_count==nblk){
          con_stat=ON;
          mat_stat=ON;
          free(blk_nelmt);
          free(blk_nenod);
          printf("complete!\n");
          fclose(outf_con);
          fclose(outf_mat);
        }
        continue;
      }
    }
  }
  
}
fclose(inf);

/* write coordinates file */
printf("writing coordinates...");
sprintf(outfname,"%snodes_coords_file",outhead);
outf_coord=fopen(outfname,"w");
fprintf(outf_coord,"%d\n",nnode);
for(i=0;i<nnode;i++)fprintf(outf_coord,"%d %.15lf %.15lf %.15lf\n",i+1,coord[0][i],coord[1][i],coord[2][i]);
fclose(outf_coord);
printf("complete!\n");

/* Write sidesets */
if(nss>0 && ss_stat!=ON){
  printf("writing side sets...");
  count_sselmt=0;
  nface_change=0;
  for(i=0;i<nss;i++){
    ndir=1;
    sprintf(outfname,"%s%s",outhead,ss_name[i]);
    outf_side=fopen(outfname,"w");
    fprintf(outf_side,"%d\n",ss_nside[i]);
    /* assign normal */
    if (strstr(ss_name[i],"xmin")!=NULL){
      normal[0]=-1.0; normal[1]=0.0; normal[2]=0.0;
    }else if(strstr(ss_name[i],"xmax")!=NULL){
      normal[0]=1.0; normal[1]=0.0; normal[2]=0.0;
    }else if(strstr(ss_name[i],"ymin")!=NULL){
      normal[0]=0.0; normal[1]=-1.0; normal[2]=0.0;
    }else if(strstr(ss_name[i],"ymax")!=NULL){
      normal[0]=0.0; normal[1]=1.0; normal[2]=0.0;
    }else if(strstr(ss_name[i],"zmin")!=NULL || strstr(ss_name[i],"bottom")!=NULL){
      normal[0]=0.0; normal[1]=0.0; normal[2]=-1.0;
    }else if(strstr(ss_name[i],"zmax")!=NULL || strstr(ss_name[i],"top")!=NULL){
      normal[0]=0.0; normal[1]=0.0; normal[2]=1.0;
    }else{
      /* do not check normal for unknown-normal faces */
      isnorm=0;
      ndir=1;
    }
    for(j=0;j<ss_nside[i];j++){
      fprintf(outf_side,"%d ",ss_elmt[count_sselmt]);
      ielmt=ss_elmt[count_sselmt]-1; /* index starts from 0!*/
      iface=ss_side[count_sselmt]-1; /* index starts from 0!*/      
      for(k=0;k<4;k++){
        gfnod[k]=elmt_node[fnod[iface][k]-1][ielmt]; /* index starts from 0!*/
        igfnod=gfnod[k]-1; /* index starts from 0!*/
        p[0][k]=coord[0][igfnod]; p[1][k]=coord[1][igfnod]; p[2][k]=coord[2][igfnod];
      }
      if(isnorm){
        ndir=check_normal(p,normal);
      }
      /* printf("Hello %d %d %d %d\n",gfnod[0],gfnod[1],gfnod[2],gfnod[3]); */
      if(ndir==1){
        fprintf(outf_side,"%d %d %d %d\n",gfnod[0],gfnod[1],gfnod[2],gfnod[3]);
      }else if(ndir==-1){
        fprintf(outf_side,"%d %d %d %d\n",gfnod[0],gfnod[3],gfnod[2],gfnod[1]);
        nface_change++;
      }else{
        printf("ERROR: absurd value of ndir:%d for normal orientation!\n",ndir);
        exit(-1);
      }
      count_sselmt++;     
    }
    fclose(outf_side);
    ss_stat=ON;
  }
  printf("complete!\n");
  printf(" number of side sets written: %d\n",nss);
  printf(" number of reversed normals: %d\n",nface_change);
}
 
for(i=0;i<3;i++)free(coord[i]);
free(coord);
for(i=0;i<8;i++)free(elmt_node[i]);
free(elmt_node);
if(nss>0){                                                              
  free(ss_elmt);                                                                  
  free(ss_side);
}                                                                  
/* check status */
if(nnode!=node_countx || nnode!=node_county || nnode!=node_countz){
  printf("ERROR: number of nodes inconsistent!\n");
  exit(-1);
}
if(nelmt!=elmt_count){
  printf("ERROR: number of elements inconsistent!\n");
  exit(-1);
}
if(dim_stat!=ON){
  printf("ERROR: dimensions cannot be read!\n");
  exit(-1);
}

if(ss_stat!=ON){
  printf("WARNING: side boundary conditions cannot be read!\n");
  exit(-1);
}
if(con_stat!=ON){
  printf("ERROR: connectivity cannot be read!\n");
  exit(-1);
}
if(mat_stat!=ON){
  printf("ERROR: material IDs cannot be read!\n");
  exit(-1);
}
if(coordx_stat!=ON){
  printf("ERROR: x coordinates cannot be read!\n");
  exit(-1);
}
if(coordy_stat!=ON){
  printf("ERROR: y coordinates cannot be read!\n");
  exit(-1);
}
if(coordz_stat!=ON){
  printf("ERROR: z coordinates cannot be read!\n");
  exit(-1);
}
printf("--------------------------------\n");
return(0);
}
/*======================================*/

/* this function check the direction of normal to the given normal
and returns 1 if the direction matches, returns -1 if does not. */
int check_normal(double p[3][4],double normal[3])
{ 
  int i;
  double a[3],b[3],comp_normal[3],dot;

  for(i=0;i<3;i++){
    /* compute 2 vectors */
    a[i]=p[i][1]-p[i][0];
    b[i]=p[i][2]-p[i][1];
  }
  /* compute normal */
  comp_normal[0]=(a[1]*b[2]-a[2]*b[1]);
  comp_normal[1]=(a[2]*b[0]-a[0]*b[2]);
  comp_normal[2]=(a[0]*b[1]-a[1]*b[0]);
  dot=0.0;
  for(i=0;i<3;i++)dot=dot+comp_normal[i]*normal[i];

  if(dot>0.0){
    /* directions match */
    return(1);
  }else if(dot<0.0){
    /* directions do not match */
    return(-1);
  }else{
    printf("ERROR: perpendicular normal orientation!\n");
    exit(-1);
  }
}
