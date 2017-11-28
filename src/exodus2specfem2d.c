/** @file exodus2specfem2d.c
 *  @brief Convert ASCII exodus file to SPECFEM2D format
 *
 *  This program converts the ASCII exodus file exported from CUBIT to
 *  several input files required by the SPECFEM2D program. Currently, this
 *  program only handles the 2D quadrilateral elements with four nodes.
 *  The binary exodus file (e.g., .e file) needs to be converted into an ASCII file,
 *  generally using a free console application "ncdump" which is a part of the
 *  netCDF library, and can be downloaded from \n
 *  http://www.unidata.ucar.edu/downloads/netcdf/index.jsp. Please see the
 *  detailed steps below.
 *
 * <!-- @author Hom Nath Gharti (hgharti_AT_princeton_DOT_edu) -->
 *
 *  ## Dependencies:
 *    stringmanip.c: string manipulation routines
 *
 *  ## Compile:
 *    gcc exodus2specfem2d.c -o exodus2specfem2d -lm
 *
 *  ## Usage:
 *    exodus2specfem2d \em input_file [\em Options] \n\n
 *    Example: \n
 *    exodus2specfem2d mesh.e -bin=1 \n
 *    or \n
 *    exodus2specfem2d mesh.txt
 *
 *  ## Options:
 *
 *    - -fac: Use this option to multiply coordinates by some factor. This is
 *          important for unit
 *          conversion, e.g., to convert m to km use -fac=0.001 [DEFAULT 1]
 *
 *    - -bin: Use this option if you want to convert exodus binary directly, provided
 *          that the command "ncdump" is in the path. The command "ncdump" is a
 *          part of netCDF library that can be downloaded for free from \n
 *          http://www.unidata.ucar.edu/downloads/netcdf/index.jsp.
 *          Use -bin=1 for binary or -bin=0 for ascii file. [DEFAULT 0]
 *
 *    - -order: Use this option to check the connectivity order and make sure that
 *          the connectivity is in counterclockwise order. Use -order=1 for
 *          checking or -order=0 for no checking [DEFAULT 0].
 *
 *    - -head: Use this option to attach head of input file to output file names.
 *          Use -head=1 to attach header or -head=0 not to attach  [DEFAULT 0]
 *
 *    - -tomo: Use this option for tomography model. Since tomography model uses
 *           negative identifiers, this option will write negative block IDs.
 *          Use -tomo=1 to make negative block IDs or -tomo=0 not to make
 *          [DEFAULT 0]
 *
 *  # Basic steps starting from TRELIS:
 *
 *  ### Step 1: prepare mesh in TRELIS/CUBIT
 *  - Define material regions using "Blocks"
 *
 *    For example:\n
 *    block 1 add surface 1\n
 *    block 2 add surface 2 3
 *
 *    will assign material region 1 to surface 1 and material region 2 to
 *    surfaces 2 and 3. These material regions will be used to define material
 *    properties in "Par_file". This program will NOT generate "Par_file".
 *    The file "Par_file" must be created to run SPECFEM2D!
 *
 *  - Define element type to be QUAD4
 *
 *    For example:\n
 *    block all element type quad4
 *
 *    NOTE: If the element types are SHELL or SHELL4, "Default" or 3D option should
 *    be selected during export. If the element type is QUAD or QUAD4,	3D
 *    option should be selected. With default or 2D data, it saves only
 *    X and Y coordinates which is not always correct.
 *    Make sure that the node ordering is strictly anticlockwise
 *    (no longer necessary!) for all the elements in CUBIT.
 *
 *  - Define surface boundary conditions using "Sidesets"
 *
 *    For example:\n
 *    sideset 1 add curve 1 \n
 *    sideset 1 name 'free_surface_file'
 *
 *    will define a free or absorbing surface boundary condition on surface.
 *    Similary,\n
 *    sideset 2 add curve 3\n
 *    sideset 2 name 'absorbing_surface_file'
 *
 *    will define absorbing boundary condition on the curve 3.
 *    Note: All the above commands can also be executed using TRELIS/CUBIT GUI.
 *    "sideset 1 name 'free_surface_file'" is equivalent to
 *    clicking sideset 1 and renaming.
 *
 *  ### Step 2: export mesh file as exodus file say "mesh.e" (always use 3D option!)
 *
 *  ### Step 3: convert "mesh.e" to SPECFEM2D files
 *    exodus2specfem2d mesh.e -bin=1
 *
 *  There will be several output files:
 *
 *  - coordinates : coordinates file => total number of nodes followed by
 *    nodal coordinate ? (? -> x, y, z)
 *
 *  - connectivity : element file => total number of elements followed by
 *    connectivity list
 *
 *  - materials : material file => total number of elements followed by
 *    material IDs
 *
 *  - surface* : sourface boundary condition files => total number of elements
 *    followed by element ID and surface nodes
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

#define OFF 0
#define ON 1
/* 0 and 1; developer version */
#define TEST_JACOBIAN 0
#define MAXBULK 100000
#define MAXLINE 100
/* auxiliary routines */
void removeExtension(char *, char *);
int get_int(int *, char *, char *);
int look_int(int *, char *, char *);
int getfirstquote(char *, char *);
int shape(double,double,double**);
int check_normal(double [3][4],double [3]);
int isclockwise(int, double [], double []);
double absmaxval(int, double []);

/* main routine */
int main(int argc,char **argv){
int e,i,inode,itmp,j,k;
int i1,i2,nod1,nod2,n1,n2,n3,n4;
int bulksize;
/* geometry dimension */
int ndim;
/* number of nodes, number of elements */
int nnode,nelmt;
/* number of blocks, number of side sets */
int nblk,nss;
/* element, node count */
int elmt_count,node_count;
/* element, node count */
int node_countx,node_county,node_countz;
/* block, node set count */
int blk_count,ss_count;
/* status */
int dim_stat,ss_rstat,ss_stat,con_stat,coord_stat,mat_stat;
int side1,side2,side3,side4;
int inum,idim[3],switch_coord[3];
/* status */
int coordx_stat,coordy_stat,coordz_stat;
double x[4],z[4],s[4],t[4],**lag4;
double dx_ds,dx_dt,dz_ds,dz_dt,detJ;
/* number of elements, number of nodes per element in each node */
int *blk_nelmt,*blk_nenod;
/* number of nodes in each node set */
int *ss_nside;

/* multiplication factor for coordinates, temporary double */
double fac,dtmp;
double **coord,*xp,*zp;
int **elmt_node;
char *bulk,line[MAXLINE],token[62],dumc[250],etype[12],stag[62];
/* coordinates name */
char **coord_name;
char fonly[62],infname[62],outfname[62],outhead[62];
/* node set names */
char **ss_name;

/* total number of element in absorbing BC and free surface */
int count_sselmt,count_ssside,sumss_nelmt,*ss_elmt,*ss_side;

/* test if binary */
int isbin;
/* option to add innput file header in the output file names */
int ishead;
/* option to convert block IDs to negative for tomography models */
int istomo;
/* normal direction, test if normal has to checked */
int ndir,isorder;

FILE *inf,*outf_mat,*outf_con,*outf_coord,*outf_side;
int   norder[4]={0,3,2,1};
int   isdone,isflag,ncount,ncw,nccw;
/* default factor and binary switch*/
fac=1.0; isbin=OFF; isorder=OFF; ishead=OFF; istomo=OFF; ndir=1;
for(i=0;i<3;i++)switch_coord[i]=ON;

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
  }else if(look_int(&itmp,"-order=",argv[i])==0){
    isorder=itmp;
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
printf("isorder: %d\n",isorder);
printf("ishead: %d\n",ishead);

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
  if (system(dumc)!=0){
  printf("ERROR: command \"%s\" cannot be executed! use -bin=0 or no option \
          for ascii input file! \n",dumc);
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
      /* each name has maximum of 62 characters */
      coord_name[i]=malloc(62*sizeof(char));
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
    elmt_node=malloc(4*sizeof(int *));
    if(elmt_node==NULL){
      printf("ERROR: not enough memory!");
      exit(-1);
    }
    for(i=0;i<4;i++){
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
        /* each name has maximum of 62 characters */
        ss_name[i]=malloc(62*sizeof(char));
      }
      ss_nside=malloc(nss*sizeof(int));
    }

    if(nss>0){
      /* This segment has a significance only if nss has legitimate value */
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
    //printf("saving side BCs...");
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
          /* read comma separated data */
          fscanf(inf,"%d,",&ss_elmt[count_sselmt]);
          count_sselmt++;
        }

        fscanf(inf,"%s",dumc); /* ; */
        fscanf(inf,"%s",token);
        sprintf(stag,"side_ss%d",i+1);
        if(strcmp(token,stag)==0){
          fscanf(inf,"%s",dumc); /* = */
          for(j=0;j<ss_nside[i]; j++){
            /* read comma separated data */
            fscanf(inf,"%d,",&ss_side[count_ssside]);
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
          printf("storing connectivity and writing material IDs...");

          sprintf(outfname,"%smaterials",outhead);
          outf_mat=fopen(outfname,"w");
          /*fprintf(outf_mat,"%d\n",nelmt);*/
        }

        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<blk_nelmt[i];j++){
          /*fprintf(outf_con,"%d ",elmt_count+1);*/
          for(k=0;k<blk_nenod[i];k++){
            fscanf(inf,"%d,",&itmp);
            elmt_node[k][elmt_count]=itmp;
          }
          /* 2D format has only 1 column */
          if(istomo){
            fprintf(outf_mat,"%d\n",-(i+1));
          }else{
            fprintf(outf_mat,"%d\n",i+1);
          }
          elmt_count++;
        }

        if(blk_count==nblk){
          con_stat=ON;
          mat_stat=ON;
          printf("complete!\n");
          fclose(outf_mat);
        }
        continue;
      }
    }

  }

}
fclose(inf);

/* free memory */
for(i=0;i<ndim;i++)free(coord_name[i]);
free(coord_name);

/* write coordinates file */
printf("writing coordinates...");
inum=0;
for(i=0; i<ndim; i++)idim[i]=-9999;
for(i=0; i<ndim; i++){
  if(absmaxval(nnode,coord[i])==0.){
    switch_coord[i]=OFF;
    if(i==0)printf("X coordinate switched OFF!\n");
    if(i==1)printf("Y coordinate switched OFF!\n");
    if(i==2)printf("Z coordinate switched OFF!\n");
    continue;
  }
  if(inum>2){
    fprintf(stderr,"ERROR: number of active dimension must be 2!\n");
    fprintf(stderr,"HINT: run exodus2specfem3d instead for 3D models!\n");
    exit(-1);
  }

  idim[inum]=i;
  inum++;
}
/*if(ndim==3)switch_coord[1]=OFF;*/
/* for SHELL element make Y coordinates OFF. TODO: add option for this. */
sprintf(outfname,"%scoordinates",outhead);
outf_coord=fopen(outfname,"w");
fprintf(outf_coord,"%d\n",nnode);
for(j=0;j<nnode;j++){
  for(i=0;i<ndim;i++){
    /* Do not write Y coordinates */
    if(switch_coord[i]!=OFF)fprintf(outf_coord,"%.6f ",fac*coord[i][j]);
  }
  fprintf(outf_coord,"\n");
}
fclose(outf_coord);
printf("complete!\n");

/* Connectivity */
if(nblk>0 && con_stat==ON){
  xp=(double *)malloc(4 * sizeof(double));
  zp=(double *)malloc(4 * sizeof(double));

  if(isorder==0){
    printf("connectivity order: preserve!\n");
  }else if(isorder==1){
    printf("connectivity order: check for clockwise ordering!\n");
  }
  printf("saving connectivity...");
  sprintf(outfname,"%sconnectivity",outhead);
  outf_con=fopen(outfname,"w");
  fprintf(outf_con,"%d\n",nelmt);

  /* WARNING: only for 4 noded elements */
  ncount=0;
  ncw=0; nccw=0;
  for(i=0; i<nelmt; i++){
    /*printf("Before coordX: %6.2f %6.2f %6.2f %6.2f\n",xp[0],xp[1],xp[2],xp[3]);
 *     printf("Before coordZ: %6.2f %6.2f %6.2f %6.2f\n",zp[0],zp[1],zp[2],zp[3]);*/
    n1=elmt_node[0][i]; n2=elmt_node[1][i];
    n3=elmt_node[2][i]; n4=elmt_node[3][i];
    if(isorder==1){
      /* check ordering, and if found clockwise convert them to
         counterclockwise */
      for(j=0; j<4; j++){
        inode=elmt_node[j][i]-1;
        xp[j]=coord[idim[0]][inode];
        zp[j]=coord[idim[1]][inode];
      }
      isflag=isclockwise(4,xp,zp);
      if(isflag>0){
        /* found clockwise */
        ncw+=1;
        /* convert to counterclockwise */
        n1=elmt_node[0][i]; n2=elmt_node[3][i];
        n3=elmt_node[2][i]; n4=elmt_node[1][i];
      }
      if(isflag<0)nccw+=1;
    }
    /* preserve the order */
    fprintf(outf_con,"%d %d %d %d\n",n1,n2,n3,n4);
  }

  con_stat=ON;
  printf("complete!\n");
  fclose(outf_con);
  free(xp);
  free(zp);
  free(blk_nelmt);
  free(blk_nenod);
}
if(isorder==1){
  printf("clockwise-ordered elements: %d\n",ncw);
  printf("anticlockwise-ordered elements: %d\n",nccw);
}

/* Write boundary files */
/* Side numbering:
 * for QUAD elements ndim=2 or 3, and side numbering - 1,2,3,4 but
 * for SHELL element ndim=3, and side numbering - 3,4,5,6 */
side1=1; side2=2; side3=3; side4=4;
/* side1=1; side2=4; side3=3; side4=2; */
if(ndim==3 && strstr(etype,"SHELL")!=NULL){
  /* SHELL elements. Need to check carefully whether true for all SHELL
 *      elements, and how they number. */
  side1+=2; side2+=2; side3+=2; side4+=2;
}
if(nss>0 && ss_stat!=ON){
  printf("writing boundary files...");
  i1=0; i2=0;
  for(i=0;i<nss;i++){
    sprintf(outfname,"%s%s",outhead,ss_name[i]); /*,i+1);*/
    outf_side=fopen(outfname,"w");
    fprintf(outf_side,"%d\n",ss_nside[i]);
    i2+=ss_nside[i];
    for(j=i1;j<i2;j++){
      if(ss_side[j]==side1){
        nod1=1; nod2=2;
      }else if(ss_side[j]==side2){
        nod1=2; nod2=3;
      }else if(ss_side[j]==side3){
        nod1=3; nod2=4;
      }else if(ss_side[j]==side4){
        nod1=4; nod2=1;
      }else{
        fprintf(stderr,"ERROR: wrong side ID: %d for boundary!\n",ss_side[j]);
        fprintf(stderr,"HINT: make sure that all the Blocks have \"QUAD4\" Element type\n");
        exit(-1);
      }

      e=ss_elmt[j];
      if(strstr(ss_name[i],"absorbing")!=NULL){
        fprintf(outf_side,"%d %d %d %d %d\n",e,2,elmt_node[nod1-1][e-1],
        elmt_node[nod2-1][e-1],ss_side[j]); /* This is only for edge */
      }else{
        fprintf(outf_side,"%d %d %d %d\n",e,2,elmt_node[nod1-1][e-1],
        elmt_node[nod2-1][e-1]); /* This is only for edge */
      }
    }
    fclose(outf_side);
    i1+=ss_nside[i];
    ss_stat=ON;
  }
  printf("complete!\n");
  printf(" number of side sets written: %d\n",nss);
}else{
  printf("WARNING: no boundary files written!\n");
}

if(TEST_JACOBIAN==1){
  /* Check for negative jacobian */
  lag4 = (double **)malloc(3 * sizeof(double *));
  for(i = 0; i < 3; i++)lag4[i] = (double *)malloc(4 * sizeof(double));

  s[0]=-1.0; t[0]=-1.0;
  s[1]= 1.0; t[1]=-1.0;
  s[2]= 1.0; t[2]= 1.0;
  s[3]=-1.0; t[3]= 1.0;
  ncount=0;
  for(i=0; i<nelmt; i++){
    n1=elmt_node[0][i]; n2=elmt_node[1][i];
    n3=elmt_node[2][i]; n4=elmt_node[3][i];
    isdone=0;
    for(j=0; j<4; j++){
      shape(s[j],t[j],lag4);
      printf("%f %f %f %f\n",lag4[0][0],lag4[0][1],lag4[0][2],lag4[0][3]);
      printf("%f %f %f %f\n",lag4[1][0],lag4[1][1],lag4[1][2],lag4[1][3]);
      printf("%f %f %f %f\n",lag4[2][0],lag4[2][1],lag4[2][2],lag4[2][3]);

      /* Anticlockwise mapping */
      x[0]=coord[idim[0]][n1]; x[1]=coord[idim[0]][n2];
      x[2]=coord[idim[0]][n3]; x[3]=coord[idim[0]][n4];

      z[0]=coord[idim[1]][n1]; z[1]=coord[idim[1]][n2];
      z[2]=coord[idim[1]][n3]; z[3]=coord[idim[1]][n4];
      /*
      if(ndim==2){
        z[0]=coord[1][n1]; z[1]=coord[1][n2];
        z[2]=coord[1][n3]; z[3]=coord[1][n4];
      }else{
        z[0]=coord[2][n1]; z[1]=coord[2][n2];
        z[2]=coord[2][n3]; z[3]=coord[2][n4];
      }
      */

      dx_ds=0.0; dx_dt=0.0; dz_ds=0.0; dz_dt=0.0;
      for(k=0; k<4; k++){
        dx_ds+=x[k]*lag4[1][k]; dx_dt+=x[k]*lag4[2][k];
        dz_ds+=z[k]*lag4[1][k]; dz_dt+=z[k]*lag4[2][k];
      }
      detJ=dx_ds*dz_dt-dx_dt*dz_ds;
      printf("Before -element: %d j: %d %f\n",i+1,j,detJ);

      if(!isdone && detJ<=0.0){
        ncount+=1;
        isdone=1;
      }

      if(detJ <= 0.0){
        /* printf("Negative Jacobian: %f for element: %d\n",detJ,i+1); */
        /* Clockwise mapping */
        x[0]=coord[idim[0]][n1]; x[1]=coord[idim[0]][n4];
        x[2]=coord[idim[0]][n3]; x[3]=coord[idim[0]][n2];

        z[0]=coord[idim[1]][n1]; z[1]=coord[idim[1]][n2];
        z[2]=coord[idim[1]][n3]; z[3]=coord[idim[1]][n4];
        /*
        if(ndim==2){
          z[0]=coord[1][n1]; z[1]=coord[1][n2];
          z[2]=coord[1][n3]; z[3]=coord[1][n4];
        }else{
          z[0]=coord[2][n1]; z[1]=coord[2][n2];
          z[2]=coord[2][n3]; z[3]=coord[2][n4];
        }
        */

        dx_ds=0.0; dx_dt=0.0; dz_ds=0.0; dz_dt=0.0;

        for(k=0; k<4; k++){
          dx_ds+=x[k]*lag4[1][norder[k]]; dx_dt+=x[k]*lag4[2][norder[k]];
          dz_ds+=z[k]*lag4[1][norder[k]]; dz_dt+=z[k]*lag4[2][norder[k]];
        }
        detJ=dx_ds*dz_dt-dx_dt*dz_ds;
        printf("After - element: %d j: %d %f\n",i+1,j,detJ);
        if(detJ <= 0.0){
          printf("After - element: %d j: %d %f\n",i+1,j,detJ);
        }
      }
    }
    /* exit(-1); */
  } /*for(i=0; i<nelmt; i++){ */


  printf("ncount:%d",ncount);
  free(lag4);
  free(elmt_node);
} /* if(TEST_JACOBIAN==1){ */


for(i=0;i<3;i++)free(coord[i]);
free(coord);
for(i=0;i<4;i++)free(elmt_node[i]);
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

/* shape functions */
#define quart 0.25
int shape(double s, double t, double **lag4)
{
double sm,sp,tm,tp;

printf("In function: %f %f\n",s,t);
sp = s + 1.0; sm = s - 1.0;
tp = t + 1.0; tm = t - 1.0;

/*lag4 = (double **)malloc(3 * sizeof(double *));
 * for(i = 0; i < 3; i++)lag4[i] = (double *)malloc(4 * sizeof(double));*/

/* Shape functions */
lag4[0][0]=quart*sm*tm; lag4[0][1]=-quart*sp*tm;
lag4[0][2]=quart*sp*tp; lag4[0][3]=-quart*sm*tp;

/* Derivatives with respect to s */
lag4[1][0]=quart*tm; lag4[1][1]=-quart*tm;
lag4[1][2]=quart*tp; lag4[1][3]=-quart*tp;

/* Derivatives with respect to t */
lag4[2][0]=quart*sm; lag4[2][1]=-quart*sp;
lag4[2][2]=quart*sp; lag4[2][3]=-quart*sm;

return 0;
}
/*---------------------------------------------------------------------------*/

/* check whether the nodes are clockwise */
#define SUMEDGE 0
#define NORMAL 1
#define APPROACH NORMAL /* SUMEDGE or NORMAL */
int isclockwise(int n, double x[n], double z[n])
{
int i,iflag;
double sum_edge;
double norm,nvec[3],v1[3],v2[3];

if(APPROACH == SUMEDGE){
  /* sum over edges approach */
  sum_edge=0.0;
  for(i=0; i<n-1; i++){
    sum_edge+=(x[i+1]-x[i])*(z[i+1]+z[i]);
  }
  sum_edge+=(x[0]-x[n-1])*(z[0]+z[n-1]);
  if(sum_edge>0){
    iflag=1;
  }else{
    iflag=-1;
  }
}else{
  /* normal approach */
  v1[0]=  x[1]-x[0]; v1[1]=0.0; v1[2]=  z[1]-z[0];
  v2[0]=x[n-1]-x[0]; v2[1]=0.0; v2[2]=z[n-1]-z[0];

  nvec[0]=v1[1]*v2[2]-v1[2]*v2[1];
  nvec[1]=v1[2]*v2[0]-v1[0]*v2[2];
  nvec[2]=v1[0]*v2[1]-v1[1]*v2[0];
  norm=sqrt(nvec[0]*nvec[0]+nvec[1]*nvec[1]+nvec[2]*nvec[2]);
  nvec[0]/=norm;
  nvec[1]/=norm;
  nvec[2]/=norm;
  if(nvec[1]>0.0){
    iflag=1;
  }else if(nvec[1]<0.0){
    iflag=-1;
  }else{
    fprintf(stderr,"ERROR: Zero jacobian! %f %f\n",nvec[1],norm);
    fprintf(stderr,"coordX: %6.2f %6.2f %6.2f %6.2f\n",x[0],x[1],x[2],x[3]);
    fprintf(stderr,"coordZ: %6.2f %6.2f %6.2f %6.2f\n",z[0],z[1],z[2],z[3]);
    exit(-1);
  }
}
return(iflag);
}
/*---------------------------------------------------------------------------*/

/* Get absolute maximum value of an array */
double absmaxval(int n, double x[n])
{
int i;
double absx,xmax;

xmax=abs(x[0]);
for(i=1; i<n-1; i++){
  absx=abs(x[i]);
  if(absx>xmax)xmax=absx;
}
return(xmax);
}
/*---------------------------------------------------------------------------*/
