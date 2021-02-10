#include "struct.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

int get_store_ptr(int , name );

/* ---------------------------------------- */
/*     print error message and exit         */
/* ---------------------------------------- */
void nerror(char *error_text, int nn)
/*  standard error handler */
{
printf("*** Error*****\n\t%s :%d\n\n",error_text,nn);
printf("... now exiting ...\n");
exit(0);
}

/* ---------------------------------------- */
/*    returns position within string        */
/* ---------------------------------------- */
int strpos(char* string,char* search)
 {
 char *pos;
 pos=strstr(string,search);
 if (pos != NULL)
    pos=(char*) (pos-string);
 else
    pos=(char*) -1;

 return (int) pos;
 }

/* --------------------------------------------------- */
/*     transfers index to a parameter vector	       */
/* --------------------------------------------------- */
char* set_parvector(unsigned long dv, int out)
{
int d,IsLog;
unsigned long dvl=1,vi;
double f,ff,fm,fac,fr,max;
char *ind,tmp[99];
varvalues[0]='\0';
if(out)
  printf("%ld: ",dv);
for(d=0;d<num_variat;d++)
  {

  vi=dv%(variat_steps[d]*dvl);
  vi=(unsigned long)(vi*1.0/dvl);
  if(out &&0)
     printf("%ld %ld %d\t",dv,vi,variat_steps[d]);
  dv-=vi*dvl;
  dvl*=variat_steps[d];
  f=variat_delt[d]*(variat_steps[d]-1);
  fm=variat_min[d];
  max = fm+f;
  if(fm<0.1*max && fm>1E-8 && 1)/* && variat_steps[d]==3 */
      {
      fr=(double)vi/(variat_steps[d]-1);
      ff=log10(max/(fm+1E-9));
      fac=pow(10.0,fr*ff);
      *par_val[d]=fm*fac;
      IsLog=1;
      }
  else
     *par_val[d]=fm+vi*variat_delt[d],IsLog=0;
  if(out)
    printf("%s:%1.4e\t",variat_names[d],*par_val[d]);
  ind=index((char *)varvalues,'\0');
  if((max<0.1) || IsLog)
     sprintf(ind,"%1.3e\t",*par_val[d]);
  else if(max<1.1)
     sprintf(ind,"%1.3f\t",*par_val[d]);
  else if(max<10.1)
     sprintf(ind,"%1.2f\t",*par_val[d]);
  else
     sprintf(ind,"%1.2f\t",*par_val[d]);

  //sprintf(pres,"%s%1.2e",variat_names[d],*par_val[d]);
  sprintf(tmp,"%s%1.2e",variat_names[d],*par_val[d]);
  if (d==0) strcpy(pres,tmp);
  else strcat(pres,tmp);
  }
if(out)
  printf("\n");

return pres;
}

/* ------------------------------------ */
/*    Re-setting changes to zero        */
/* ------------------------------------ */
void reset_derivs(void) {
int i,j;
long int number;
number=num_states*n_comp*sizeof(double);
for (j=0;j<num_states;j++)
  for (i=0;i<n_comp;i++)
    SCCC [j][i]=0.0;
}

/* --------------------------------------------- */
/*                                               */
/*    Re-setting initial values from matrix      */
/*                                               */
/* ------------------------------- -------------- */
void init_states(int dir)
{
int i,ii,j,j0,d,di,da,d1,experim;

if(Experiment<=-1 && ind_exp<LengthOfIndExp) experim=IndExp[ind_exp];
else   experim=Experiment;

/* -------------------------------------------- */
/*    set initial values from fixed stores      */
/* -------------------------------------------- */
if(dir==1)
  for(d=0;d<num_states;d++)
    for (j=0; j<n_comp; j++) CCC[d][j]=CCC_i[d][j];
/* -------------------------------------------------- */
/*  set fixed stores from initial value in init.sta   */
/* -------------------------------------------------- */
else if(dir==0)
  for(d=0;d<num_states;d++)
    for (j=0; j<n_comp; j++) CCC_i[d][j]=CCC[d][j];
}

/* ---------------------------------------------------- */
/*                                                      */
/*   writes last state variables in XSISI format ..     */
/*                                                     */
/* ---------------------------------------------------- */
void write_last_state()
{
int d,dd,di;
FILE* sp;

printf("writing last state to last.sta ...\t");
printf("%d %d\t\t\t%s\n",num_states,n_comp,MovieResFile);
sp=fopen("last.sta","w");

for(d=0;d<num_states;d++)
  {
  write_header(sp,state_names[d],n_comp,0);
  for(dd=0;dd<n_comp;dd++)
    fprintf(sp,"%1.4f\t",CCC[d][dd]);
  fprintf(sp,"\nend\n");
  }
fclose(sp);
}

/* --------------------------------------------- */
/*                                              */
/*    writes data header in XSISI format ..     */
/*                                              */
/* -------------------------------------------- */
void write_header(FILE *sp,char *nam,int dim1, int dim2)
{
fprintf(sp,"array\t %s\n",nam);
fprintf(sp,"\t\\u umol/ccm\n");
fprintf(sp,"\t typeOfArray\t float\n");
if(dim2==0)
  fprintf(sp,"\t dimension\t %d\n",dim1);
else
  fprintf(sp,"\t dimension\t %d %d\n",dim1,dim2);
fprintf(sp,"data\n\t");
}

void make_store(double **store_ptr,float *vector,int number)
{
int i;
for (i=0; i<number; i++)
   vector[i]=(float) (*store_ptr[i]);
}

/*-------------------------------------------------*/
/*                                                 */
/*    vorbereitung der ausgabe in ergebnisdatei    */
/*                                                 */
/*-------------------------------------------------*/
int store_prep_open(char *output_file)
{
int ok,error,d,dd,i,ne,nem,nsf;
char helpstr[20],*cptr,*sdot;
char outfn[200];
FILE *tmp;

struct header
  {     unsigned int n_stor; float startim;  float endtim;  float maxdelt;
     float outdelt; float year; float cycle; char sdir[12];
  } oh;

/* ---------------------------------- */
/*       Allocate Buffer memory       */
/* ---------------------------------- */
nsf=(num_stores+1)*sizeof(float);
store_vector=(float *) malloc((n_comp+2)*nsf);
store_ptr=(double **) malloc((num_stores+1)*(n_comp+1)*
                                        sizeof(double *));
error=0;
/* ---------------------------------- */
/*     Get double address from name   */
/* ---------------------------------- */
for(d=0;d<num_stores;d++)
   {
   ok=1;
   if(get_store_ptr(0,store_names[d])==0)
     {
     printf("%s not in varlist:\n\n",store_names[d]);
     ok=0;
     for (i=0;i<num_states;i++)
        printf("%s\t",state_names[i]);
     printf("\n");
     }
   if(ok==1)
      for(dd=0;dd<n_comp;dd++)
	 store_ptr[d*n_comp+dd]=&dptr[dd];
   else
     error=1;
/*   printf("%s %d/%d: %1.3f %1.3f\t%1.3f %1.3f\n",store_names[d],d,n_comp,
             dptr[0],dptr[2],*(store_ptr[d*n_comp]),*(store_ptr[d*n_comp+2]));*/
   }

/*------------------------------------------*/
/* directory filtern aus modpath und in
   den result-header                        */
/*------------------------------------------*/

strncpy(oh.sdir,"Multi",12);

/*----------------------------------------------*/
/*  result-header information about simulation  */
/*----------------------------------------------*/
oh.startim = (float) storetim;
oh.endtim = (float) TimeEnd;
oh.maxdelt = (float) TimeStep;
oh.outdelt = (float) OutputStep;
oh.year = (float) 00;
oh.cycle = (float) 365.0;

/*----------------------------------------------*/
/*   open result file and writes header infos   */
/*----------------------------------------------*/
if(Experiment<0)  /* open array of files in case of loop
      			over experiments */
    nem=-Experiment;
else
  nem=1;

for(ne=0;ne<nem;ne++)
  {
  strcpy(outfn, output_file);
  if(Experiment<0) {
  sdot=strchr(outfn,'.');
  //sprintf(sdot-1,"%s.outb\0",ExpName[ne]);
  sprintf(sdot-1,"%d.outb\0",ne);
}
  //printf("opening (%s) %s ...\n",output_file,outfn);
  if( (outfile=fopen(outfn,"wb"))==NULL) error=1;

  outfile_exp[ne]=outfile;
  oh.n_stor= (int)((num_stores)*(n_comp+0));

  fwrite(&oh,sizeof(oh),1,outfile);
  for (i=0;i<num_stores;i++)
    for (dd=0;dd<n_comp+0;dd++)
     {
     sprintf(helpstr,"%s(%d)",store_names[i],dd);
     for (d=(int)strlen(helpstr);d<12;d++)
       helpstr[d]= (char) 32;
     helpstr[12]= (char) NULL;
     fwrite(helpstr,12,1,outfile);
     }
  }
return error;
}

/*----------------------------------------------*/
/*   find address of variale
                      according name string     */
/*----------------------------------------------*/
int get_store_ptr(int box, name search) {
double *ptr;
int i,j,sp,ok=0;
for (i=0;i<num_states;i++)
  {
  sp=strpos(search,state_names[i]);
  if (sp>0 && sp<2)
   {
   if (strpos(search,"S") == 0) dptr=&SCCC[i][box],ok=1;
   }
  else
    if (sp==0) dptr=&CCC[i][box],ok=1;
   }
if(ok==0)
  for (i=0;i<num_others;i++)
    if (strcmp(search,other_names[i]) == 0) dptr=other_ptr[i],ok=1,j=i;
return ok;
}

/* ------------------------------------------ */
/*        simple random number generator      */
/*          depends on bit-architecture       */
/* ------------------------------------------ */
double random2()
{
if((l0=l0*65539)<0)
        l0+=l1;
if(sizeof(long)==4)
   return (double)(l0)*(double)2.3283064369E-10;
else
   return (double)(l0)*2.3283064369E-10*2.3283064369E-10;
}


/* ------------------------------------------------------ */
/*   writes random number to file for statistical tests   */
/* ------------------------------------------------------ */
void random_out(int r_init, unsigned long num)
{
FILE *sp;
char RandomFile[21]="RandomNumbers_0.res";
unsigned long  d,dd;
double r0,r1;
RandomFile[14]=48+(r_init%10);
l0=(r_init*4712)%(513+r_init*r_init);
printf("writing %ld random numbers in %s ...\n",num,RandomFile);

sp=fopen(RandomFile,"w");

r0=random2();
for(d=0;d<num;d++)
   {
   r1=random2();
   fprintf(sp,"%1.5f %1.5f\n",r0,r1);
   r0=r1;
   }
printf("\t\t\tready\n");
fclose(sp);
}

/* ------------------------------------------------------ */
/*      astronomical formula for daylength (TODO source)  */
/* ------------------------------------------------------ */
double daylength(double day, double lat)
{
//   lat = latitude
//   day of the year
double f,fa,fb,pi=3.1415,dl;

f  = asin(0.39795*cos(0.2163108 + 2*atan(0.9671396*tan(0.00860*(day-186)))));
fa = sin(0.8333*pi/180) + sin(lat*pi/180)*sin(f);
fb = cos(lat*pi/180)*cos(f);
dl = 1. - (1./pi)*acos(fa/fb);
//printf(" \n%1.3f %1.3f   %1.3f\n",fa/fb,acos(fa/fb),dl);
return dl;
}
