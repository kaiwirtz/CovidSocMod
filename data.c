/*-------------------------------------------------------*/
/*  Routines for data handling and model validation      */
/*                                                       */
/* Author	 Kai Wirtz (kai.wirtz@hzg.de)                  */
/*-------------------------------------------------------*/
#include "struct.h"
#define CMP_ALL 0
double var[NUM_FIT_VARS],scalfac[2]={1E-6,1.};
char NameData[NUM_FIT_VARS][12]={"fatal","Mobility"};
void dev_clean();
int NumData,mask,d_logval[NUM_FIT_VARS],act_dat,ds_count,di_start[NUM_FIT_VARS];
double sim_max[NUM_FIT_VARS],data_max[NUM_FIT_VARS],mean[NUM_FIT_VARS],InitFatalfac[NUM_FIT_VARS],dat_avg,sim_avg,dat_var,sim_var;

/* ----------------------------------------------------- */
/*             Preparation of data array                 */
/* ----------------------------------------------------- */
long data_prep(int mode)
{
int d,dd,ti,da,da0,di,de,partyp,act_tim,dl[NUM_FIT_VARS],dlv[NUM_FIT_VARS],i0,dii,ds=0;
double f,ff; //data_max[NUM_FIT_VARS]
double x,dx,data_min[NUM_FIT_VARS],mean2[NUM_FIT_VARS],meanv[NUM_FIT_VARS],norm_dat_i[NUM_FIT_VARS];

// set country specific correction factor for data
act_dat=ColumnsOfdata-1;
if (mode==0) {
/* ----------------------------------------------------- */
/*    calculates renormalisation factors from            */
/*                       1) min & max                    */
/*                       2) mean standard deviation      */
/* ----------------------------------------------------- */
da0=DataCountryIndex-1;

// after correction, Iran already started before late Jan 2020
if(DataCountryIndex==4) dat_off=20;
else dat_off=0;

/* --------------------------------------------- */
/*    loop over data entries:
      scan all values for range and statistics   */
/* --------------------------------------------- */
for(da=da0,de=ti=0;da<act_dat;de++,da+=ncountry)
  {
  data_max[da]=-1,sim_max[da]=-1;data_min[da]=9E9;
  mean[da]=meanv[da]=var[da]=0;
  dl[da]=0; di_start[da]=-9;
  /* bulk variables are compared as log -> d_logval[da]=1 */
  d_logval[da]=0; // 0
/* --------------------------------------------------------- */
/*   loop over all rows of data matrix, with 2 formats ...   */
/* --------------------------------------------------------- */
  for(di=NumData=0;di <(1-CMP_ALL)*180+CMP_ALL*RowsOfdata; di++)
    {
  /* --------------------------------------------------- */
  /*  index format w 2nd entry: experiment/region index  */
  /* --------------------------------------------------- */
    if(di>0)
      if( data[di][0]< data[di-1][0])
        printf("wrong order:%d \t %1.3f %1.3f\n",di,data[di][0],data[di-1][0]);

    f=data[di][1+da];
    if (data[di][1+da0]>=critfatal*1E6)
      {
      if (di_start[da]<0) di_start[da]=di;
      if (d_logval[da]==1)  f=log10(f*scalfac[de]+1E-7);
      if(f>data_max[da]) data_max[da]=f; // min and max
      if(f<data_min[da]) data_min[da]=f;
      mean[da]+=f; // mean
      dl[da]++;
      var[da]+=f*f;// variance
      }
    }
  if(sim_out>-1 & 0) printf("\nfound %d/%d in data-set\t act_dat=%d\n",NumData+1,dl[da],act_dat);
  /* -------------------------------------------------- */
  /*    calculating min&max-> normalization factor
  				                       for each variable      */
  /* -------------------------------------------------- */
  if(dl[da]>0)
    {
    mean[da]/=dl[da];
    var[da]/=dl[da];
    var[da]-=mean[da]*mean[da];
    if (d_logval[da]==1)
      {
      norm_dat[da] = sqrt(var[da]);
      if(norm_dat[da]<0.25) norm_dat[da]=0.25;
      }
    else // add base variance to homogenize intercomparison between regions
      norm_dat[da] = sqrt(var[da])*scalfac[de]+(de==1)*30;
    //norm_dat[da] = mean[da];
    if(norm_dat[da]<EPS*scalfac[de]) norm_dat[da]=EPS*scalfac[de];
    InitFatal=data[di_start[da0]][1+da0];
    // norm_dat[da] = 1.0;
  if(sim_out>=10)
    printf("\033[1m\033[34m %d %s: %d avg=%1.3f\tnorm=%1.2e\t%1.2f-%1.2f\tstart at %d  with %1.2f\033[0m\n", da,HeadOfdata[da+1],dl[da],mean[da],norm_dat[da],  data_min[da],data_max[da],di_start[da],InitFatal);
    }
  else
  if(sim_out&& 1) printf("%d %s: No %s data found!\n",da,HeadOfdata[da+1],NameData[ds]);
  }
/* ---------------------------------------------------------- */
/*   create pointer vector to model variables for comparison  */
/* ---------------------------------------------------------- */
valid_prep(0);
/* correction factor for small individual differences in starting fatality */
InitFatalfac[da0]=1;
InitFatal0=InitFatal*1E-6; //TODO;
}

di=di_start[DataCountryIndex-1];
data_step= (long)((data[di][0]+dat_off)/TimeStep);
/* clean stat variables */
dat_avg=dat_var=sim_avg=sim_var=ds_count=0;

return di;
}

/* ---------------------------------------------- */
/*    calculating model skill: RMS from data      */
/* ---------------------------------------------- */
double calc_err()
{
int d,di,dd,de,da,da0,dp,dss;
double ff,f,fv,fg,fn,err_tot,norm,err_fmin,errf[4];

norm=err_tot=0;
da0=DataCountryIndex-1;
/* loop over all entries = regions' mortality and mobility */
for(da=da0;da<act_dat;da+=ncountry)
  {
  if(err_in[da]>EPS)
    {
  /* ------------------------------------------------------- */
  /*   calculate total error with country-specific weights   */
  /* ------------------------------------------------------- */
//    printf("calc_err %d   %d %1.1f norm=%1.1f\t%1.2f\n",da,err_in[da],err_f[da],norm_dat[da],sqrt(err_f[da]/err_in[da])/norm_dat[da]);
    fn = 1.0/norm_dat[da];
    ff = sqrt(err_f[da]/err_in[da])*fn; // Root from summed squared deviation
    if (da<ncountry)
        err_f[da]=ff;
    else
        err_f[da0]+=ff;

    err_tot +=ff;
  /* ------------------------ */
  /*   ...  normalization     */
  /* ------------------------ */
    norm+=0.5;  // equal weight of all mortality and mobility data
  if(da0==-4)  printf("calc_err %d %1.3f  %1.3f\t total=%1.3f norm %1.2e %1.1f\n",da,err_f[da],err_f[da0],err_tot,norm_dat[da],norm);
    }
  }
if(norm>EPS)
  err_tot/=norm;

/* ----------------------------------------- */
/*     evaluate match of lockdown timing     */
/* ----------------------------------------- */
ff = add_lockdown_err(da0);
err_tot += ff;
err_f[da0]+= ff;

return err_tot;
}

/* ---------------------------------------------------- */
/*   calc deviation of data series from simulation
      results   using a time screening algorithm       */
/* ---------------------------------------------------- */
int calc_dev(long di)
{
int d,de=0,da,da0,dd;
double sum,ff,f,tf,f0,f1,ff0,datat,fg,err_min,simt,err_a;
simt = Time;
/* ---------------------------------------------------- */
/*    calculating min square mean deviation from data   */
/* ---------------------------------------------------- */
da0=DataCountryIndex-1;
//printf("di=%d %d\tdata=%1.2e sim=%1.2e\t%1.2e\n",di,da,data[di][1+da],*valid_ptr[de],InitFatal);
datat=data[di][0];

/* ------------------------------------------------ */
/*     using either first half or full data set     */
/* ------------------------------------------------ */
if(simt<190+CMP_ALL*180){

/* ------------------------ */
/*      unit conversion     */
/* ------------------------ */
if (mobil0>0) InitFatalfac[da0+ncountry]=100*mobn/mobil0;
else InitFatalfac[da0+ncountry]=1E6;

/* -------------------------------------------------- */
/*    loop over data entries: mortality + mobility    */
/* -------------------------------------------------- */
for(da=da0,de=0;da<act_dat ;de++,da+=ncountry) //-(simt>42.1)*ncountry
  { // valid data & after first three days for mobility?
   if((f=data[di][1+da]*scalfac[de])>=0 && (simt>3.1 || de==0))
    {
    if (d_logval[da]==1)
      {
      f=log10(f+1E-7);
      fg=log10((*valid_ptr[de]*InitFatalfac[da]+1E-7));
      // no cut-off for outlier-data/errors TODO
      }
    else
      {
      fg=*valid_ptr[de];
      if(de==0) fg*=InitFatalfac[da];
      }
/* ------------------------------------------- */
/*      calibration specific weighing          */
/* ------------------------------------------- */
    if(CMP_ALL==0) tf=1;
    else // lower weight of late mobility data for vaccination calibration
      tf=1-((de==1)*1+1)*0.4*(simt>200)-0.5*(simt>90&&da==3);
    if(tf<0) tf=0;
    // higher weight of final mortality data for vaccination calibration
    if(de==0 && simt>300) tf=10;
/* ------------------------------------------- */
/*             simple RMS deviation            */
/* ------------------------------------------- */
    ff=f-fg;
/* ------------------------------------- */
/*   add error to cumulative mismatch    */
/* ------------------------------------- */
    err_f[da]  += tf*ff*ff;
    lim[1+de][0]=tf*ff*ff;
    err_in[da]++;
/* -------------------------- */
/*      optional output       */
/* -------------------------- */
    if(di>=di_start[da] && (simt>0)*(simt<335)*(da0==-4) ==1 &&de==0)
      {
      f1=1./scalfac[de];
      printf("%d \033[34m %d/%d tim=%1.0f:%1.0f %d\tD %1.2f ~ %1.2f \terr=%1.2f*\033[1m %1.2f\033[0m\t",da,di,di_start[da],simt,data[di][0],de,f*f1,fg*f1,tf,fabs(ff*f1));
      printf("err_f=%1.2e tot=%1.2f %d\033[0m\n",err_f[da]*f1,sqrt(err_f[da]/err_in[da])/norm_dat[da],err_in[da]);
      }
    }
  }
}
/* ----------------------------------------------------- */
/*   increment to next data time, if reached by SimTime  */
/* ----------------------------------------------------- */
if(di+1<RowsOfdata)
  return di+1;
else
  return -1;
}

/* ---------------------------------------------------- */
/*       writes header-information in var-resfile &     */
/*       sets pointer of results to be validated        */
/* ---------------------------------------------------- */
int valid_prep(int box)
{
int d,dd,du,ok,di;
FILE *sp;
valid_ptr=(double **) malloc(2*sizeof(double *));
for(d=0;d<2;d++)
/* if(err_data_weights[d]>EPS)*/
  {
  ok=1;
  if(get_store_ptr(0,NameData[d])==0) // error: variable not found in list
	   {
     printf("data %s not in varlist!\n",NameData[d]),
     ok=0;
     }
  if(ok==1) valid_ptr[d]=&dptr[box];
//  printf("valid_prep: %d %s\t%ld\t%1.4f\n",d,NameData[d],valid_ptr[d],*valid_ptr[d]);
  }
return ok;
}

/* ----------------------------------------- */
/*     evaluate match of lockdown timing     */
/* ----------------------------------------- */
double add_lockdown_err(int da0)
{
double lterr,ff;
lterr=0;
if ((da0 != 3)&&(da0 != -9)) // neither Iran nor Sweden due to lacking lockdown date
  { // high tolerance of 3 days, low one for <4.5d
  ff=fabs(locktarget-lockday[da0][1])-1;
  if(ff>3.01)  lterr=1;
  if(ff>4.51)  lterr=2,err_f[da0]+=2;
  if(ff>5.01)  lterr=100,err_f[da0]+=100;
  //if(sim_out)
  //printf("%d lockday %1.1f %1.1f vs %1.1f (sim)\terr %1.2f\n",da0,lockday[da0][0],locktarget,lockday[da0][1],ff);
  }
return lterr;
}
