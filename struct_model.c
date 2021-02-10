/*--------------------------------------------------*/
/*						   					   					   					  */
/*  definition and initialisation of
                          global variables          */
/*					                                        */
/*--------------------------------------------------*/
typedef char *String;
typedef char name[16];
#include "variables.h"
#include "defines.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
/*----------------------------------------------------*/
/* Declaration of global C variables		              */
/*                          (lists, fields, names)    */
/*----------------------------------------------------*/
char resfile[92],*cptr;

/* -------------------------------------------------- */
/*    simulation control parameters and diagnostics   */
/* -------------------------------------------------- */
double startim,endtim,start_of_int,end_of_int,cc,delt_of_int,sim_time,delt,mob_scaled,fatal_min,fac_syn,fatal_max,locktarget,norm_dat[NUM_FIT_VARS],var2[NUM_FIT_VARS],**valid_ptr,*dptr,err_f[NUM_FIT_VARS],err[N_VARIAT][NUM_DPV],*other_ptr[NUM_PD],*par_val[N_VARIAT],preruntime,weekf,InitFatal,InitFatal0,mrate_max,fac_exp,mobil0,expo,attackfac,mobil,mfatal,fsoc,Immune[10],season,fvacc,ftemp,sunf,dayl_max,fmode,flex,sub_0,SD_min,fatalavg_sim,fatalavg_dat,dat_off,Susc_old;
double Temp,simavg[NUM_FIT_VARS],simvar[NUM_FIT_VARS],errv[4],par[8],tweight;
unsigned int int_time,day_of_year,PeriodNum,SteadyState,IsTime,IsAware,mobn,ncountry,mrim,mrim0,accs;
long step_of_int,step_event,stepsum_sim,pdf_out[NUM_EXP];
double CTemp[MAX_VAR][MAX_BOX],CCC[MAX_VAR][MAX_BOX],STemp[MAX_VAR][MAX_BOX],
  *SCCC[MAX_VAR],CCC_i[MAX_VAR][MAX_BOX],CCC_1[MAX_VAR][MAX_BOX];
unsigned long l0,l1,step_total,store_step_next,target_i[4];
double SumC,Rtrans[NUM_IA],time_Hcap,fconfid0,mtran;
int  ind_old[MAX_FORCE+3],ind_new[MAX_FORCE+3],ind_end[MAX_FORCE+3],postNew,ForcedLock;
int sortind[NUM_FIT_VARS][N_OPTSIM],sortlen[NUM_FIT_VARS],nclass,data_off;
double *forcing[MAX_FORCE+3],*force_time[MAX_FORCE+3],sorterr[NUM_FIT_VARS][N_OPTSIM],*mmrates,mmrate,mrate,wave_min[4][2],wave_max[4][2];
float bet[NUM_IA][NUM_IA][20],cmob[NUM_IA][NUM_IA],adist[NUM_IA][20],lockday[20][2];

/* ------------------------------------------------- */
/*     Total rates of change in state variables      */
/* ------------------------------------------------- */
double SConc[NUM_IA][MAX_BOX],Strans[NUM_IA][NUM_IA][MAX_BOX],SaRtrans[NUM_IA][MAX_BOX],Sexposition[MAX_BOX],SCostAcc[MAX_BOX],SocialDist0[MAX_BOX],fatal[NUM_IA];
double *transM[NUM_IA][NUM_IA],transM0[NUM_IA][NUM_IA][MAX_BOX],contact[NUM_IA][NUM_IA][MAX_BOX],trans_base[NUM_IA][NUM_IA][MAX_BOX],*transO[NUM_IA][NUM_IA],mixr[NUM_IA][NUM_IA];
//,transM[NUM_IA][NUM_IA][MAX_BOX];

/* 1: Kids (-19) 2:Young adults (20-29) 3: Mid-age (30-49) 4: Elderly workers (50-59) 5: Old workers (60-69)  6: Young Seniors (70-79)  7: Old Seniors (80-)*/

/* -------------------------------------- */
/*           Auxiliary variables          */
/* -------------------------------------- */
double TotConc[MAX_BOX],dyncoup,MeanTrans[MAX_BOX],TotMort[MAX_BOX],Mtrans[NUM_IA][MAX_BOX],lim[5][MAX_BOX],fatalrate[MAX_BOX],SocialCosts[MAX_BOX],TotLoss[MAX_BOX],Contacts[MAX_BOX],Mobility[MAX_BOX],FConc[NUM_IA][MAX_BOX],InfectAge[MAX_BOX],DeathAge[MAX_BOX],rate[NUM_IA],Susc[NUM_IA],age_rat[NUM_IA][NUM_IA][MAX_BOX],TotSusc[MAX_BOX],TotRecov[MAX_BOX],SocialCostAcc[MAX_BOX],SocialCostAvg[MAX_BOX],*SocialCostMeM[MAX_BOX],*MortalityMeM[MAX_BOX],IndCosts[MAX_BOX],MortAcc[MAX_BOX][4],ExpoEff[MAX_BOX],fac[MAX_BOX],Exposure[MAX_BOX],SocialDist[MAX_BOX];
double ages[NUM_IA]={10,25,40,55,65,75,89},fage,season0,catchupf;

/* -------------------------------------------------- */
/*    Pointer to everything to be shown as result     */
/* -------------------------------------------------- */
name other_names[NUM_PD]={"TotConc","TotSusc","TotRecov","MeanTrans","fatal","SocialCosts","Contact","Mobility","ExpoEff","fac","lim0","lim1","lim2","InfectAge","DeathAge","Mtrans1","Mtrans2","Mtrans3","Mtrans4","Mtrans5","Mtrans6","Mtrans7","Exposure","SocialCostAcc","-"};
double *other_ptr[NUM_PD]={TotConc,TotSusc,TotRecov,MeanTrans,fatalrate,SocialCosts,Contacts,Mobility,ExpoEff,fac,&lim[0][0],&lim[1][0],&lim[2][0],InfectAge,DeathAge,&Mtrans[0][0], &Mtrans[1][0],&Mtrans[2][0],&Mtrans[3][0],&Mtrans[4][0],&Mtrans[5][0],&Mtrans[6][0],Exposure,SocialCostAcc};
name ExpName[NUM_EXP]={"Ref","CatchUp","SocAge","EpiAge","BothAge","SCV","ref","-"};//"fastfuture","fast","slow",

/* --------------------------------------- */
/*      Simulation control  variables      */
/* --------------------------------------- */
unsigned long store_step,num_total_variat,nvar,d_var,step_var,step_out,step_st,data_step;
int num_stores, num_others,num_variat,n_comp,variat_steps[N_VARIAT],ind_exp,
    sim_out,num_states,store_prec[NUM_PD],svar_in[NUM_FIT_VARS],
    err_in[NUM_FIT_VARS],ie,bulk[MAX_VAR],lim_i,num_tr;
float variat_min[N_VARIAT],variat_max[N_VARIAT],variat_err[N_VARIAT],variat_delt[N_VARIAT],
   variat_delt0[N_VARIAT],*store_vector;
double var_out[N_OUTVAR],var_out0[N_OUTVAR],var_out1[N_OUTVAR],**store_ptr;
name variat_names[N_VARIAT],store_names[NUM_PD],*state_names;
char pres[199],*varvalues[199], aggfile[99],aggfile2[99];
int setup; /* 0: chemo 1:mesocosm 2:reede.spring 3: rgr */
char VarFile[99];
FILE *outfile,*outfile_exp[NUM_EXP];
double	week_lag=6; // 6:start on a Monday
/* ------------------------------------------------- */
/*                                                   */
/*   creates fields for handling state variables     */
/*                                                   */
/* ------------------------------------------------- */
void struct_model(void)
{
int num,off,d,dd,d0,d1,d2,b,j,di;
double ft;
void init_states(int);

/* ------------------------------------ */
/*    check dimensions of age groups    */
/* ------------------------------------ */
if(NUM_IA > RowsOfConc )
   {
   printf("Class array not consistent\t NUM=%d/%d \n",NUM_IA,RowsOfConc);
   exit(0);
   }

n_comp=MAX_BOX;
state_names=(name *) malloc((MAX_VAR+1)*sizeof(name));
num_others=0;
while(strpos(other_names[num_others++],"-")<0);
  num_others--;
if(num_others<=0)
   printf("Error: No other names detected!\n");
/* store_names=(name *) malloc(num_stores*n_comp*sizeof(name));*/
num=(n_comp)*sizeof(double);

/*-------------------------------------------------------*/
/*      sets time offset for model data comparison       */
/*-------------------------------------------------------*/
preruntime = data_off_tim;

/* ------------------------------------------------- */
/*    state vector assignement for each variables    */
/* ------------------------------------------------- */
d0=0;
for(d=0;d<NUM_IA;d++)
  {
  (void)memcpy(&CCC[d0+d][0],Conc[d],num);
  Conc[d] = (double*) &CCC[d0+d][0];
  sprintf(state_names[d0+d],"%s%d","Conc",d+1);
  SCCC[d0+d]=&SConc[d][0];
  }
d0+=NUM_IA;
for(d=0;d<NUM_IA;d++)
  {
  (void)memcpy(&CCC[d0+d][0],aRtrans[d],num);
  aRtrans[d] = (double*) &CCC[d0+d][0];
  sprintf(state_names[d0+d],"%s%d","aRtrans",d+1);
  SCCC[d0+d]=&SaRtrans[d][0];
  }
d0+=NUM_IA;

for(d=num_tr=0;d<NUM_IA;d++)
 {
 for(dd=d;dd<NUM_IA;dd++)
   {
   CCC[d0+num_tr][0]=trans[dd][d];
   transM[d][dd] = (double*) &CCC[d0+num_tr][0];
   sprintf(state_names[d0+num_tr],"%s%d%d","trans",d+1,dd+1);
   SCCC[d0+num_tr]=&Strans[d][dd][0];
   num_tr++;
   }
 }
d0+=num_tr;
(void)memcpy(&CCC[d0][0],exposition,num);
exposition = (double*) &CCC[d0][0];
sprintf(state_names[d0],"%s","exposition");
SCCC[d0]=&Sexposition[0];
d0+=1;

(void)memcpy(&CCC[d0][0],CostAcc,num);
CostAcc = (double*) &CCC[d0][0];
sprintf(state_names[d0],"%s","CostAcc");
SCCC[d0]=&SCostAcc[0];
d0+=1;

/* ----------------------------------------------- */
/*       first check of correct  dimensions        */
/* ----------------------------------------------- */
num_states= (int)(2*NUM_IA + num_tr +2 ); /*  */
printf("num_states=%d  num_tr=%d n_comp=%d\n",num_states,num_tr,n_comp);

if(num_states>MAX_VAR)
   nerror("Number of variables exceeds MAX_VAR",num_states);

/* --------------------------------------------- */
/*   Stores initial values in a second matrix
                             for variation runs  */
/* --------------------------------------------- */
init_states(0);

/* --------------------------------------------------- */
/*     prepares Forcing data input, here NOT USED      */
/* --------------------------------------------------- */
if(ColumnsOfForcData>MAX_FORCE) printf("ColumnsOfForcData %d > MAX_FORCE\n",ColumnsOfForcData);
num = (int)(1*(RowsOfForcData+9)*sizeof(double));
/* TODO: store array size for later memory check */
for(d=0;d<ColumnsOfForcData;d++)
  {
  forcing[d] = (double *)malloc(num);
  force_time[d] = (double *)malloc(num);
  }
season0=0;

return;
}

/*--------------------------------------------------*/
/*   prepares simulation run:
      initialize boundary input, output flags, etc  */
/*--------------------------------------------------*/
void sim_prep(int ie)
{
int d,n,b,j,i,iz,fi,di,dd,da;
double ff,fg,aw,bm;
FILE *sp;

/* tartget variable init */
//for(d=0;d<4;d++)  target[d]=target_i[d]=0;

/*    setting storetim as offset	   */
storetim+=TimeStart;
/* ------------------------------------------------- */
/*   	discretization of simulation event steps		   */
/* ------------------------------------------------- */
/* adjust outdelt, if nceessary */
if (TimeStep > OutputStep )
   {
   OutputStep=TimeStep;
   if ( sim_out) {
     printf("Warning: outdelt adapted to calculation step!");
     printf("\n New outdelt=%f\n",OutputStep); }
   }
/* -------------------------------------- */
/*    calculate output/event/etc times	  */
/* -------------------------------------- */
store_step=(long) (OutputStep/TimeStep);
step_var=(long)(20.0/TimeStep);
stepsum_sim=(long)((TimeEnd-TimeStart)/TimeStep);
sim_time=0;
pdf_out[0] =0;
store_step=(long) (OutputStep/TimeStep);
step_event=(long)(95./TimeStep); // event at May,15

/*-------------------------------------------------------*/
/*      array for 21-day accumulated infection rate      */
/*-------------------------------------------------------*/
mrim=(int)(21./OutputStep);
mrim0=(int)(7./OutputStep)+1; //initial vector has to grow
//if (mrim0>mrim) mrim0=mrim;
mmrates=(double *)malloc(mrim*sizeof(double));
fatalavg_sim=0;
fatal_max=fatal_min=-1;
Susc_old=-1;
/* ------------------------------------------- */
/*    preparation of model-data comparison     */
/* ------------------------------------------- */
if(DataActive)
  {
  //sunf0 = 1 + (daylength(41,latitude[DataCountryIndex])-0.5)*sun_ampl;
  da=DataCountryIndex-1;
  err_f[da]=err_in[da]=0;
  err_f[da+ncountry]=err_in[da+ncountry]=0;
  lockday[da][1]=-9E9;
  if(DataCountryIndex>0) dyncoup=0;  //
/* ------------------------------------------------------- */
/*   continentally synchronous date of first lockdowns     */
/* ------------------------------------------------------- */
  locktarget=36;  // March 18-19
  if(da>10) locktarget = 38; //US+UK
  if(da==5) locktarget = 26; // Italy 9 March 2020
  if(da==3) locktarget = 22; // Iran
//  if(da==8) locktarget=30; // Spain earlier than reported
  }
else
  {
  locktarget=0;
  init_spread(InfectDens);
  for (i=0;i<4;i++)
    wave_max[i][0]=wave_max[i][1]=wave_min[i][1]=0,wave_min[i][0]=9E9;
  }

/*-----------------------------------------------------------*/
/*   clears/set array for time-accumulated infection rate    */
/*-----------------------------------------------------------*/
for(i=0;i<mrim;i++) mmrates[i]=0.2;  // average value for initial spread

/*---------------------------------------------------*/
/*  ratios between contact rates and transmission    */
/*---------------------------------------------------*/
for(i=0;i<NUM_IA;i++)
  for(j=0;j<NUM_IA;j++)
    mixr[i][j]=age_dist[i]/(attack[i]*attack[j]);

/* ----------------------------- */
/*      minimal SD measure       */
/* ----------------------------- */
b=0;
for(i=0,SocialDist[b]=SD_min=0;i<NUM_IA;i++)
  for(j=i;j<NUM_IA;j++) // loop over interacting age classes
    {
    aw = transM0[i][j][b] * mixr[i][j] *age_dist[i]*age_dist[j];
  /* well-being decreases with increasing social distancing */
    fg = fmax(sub_0,cmob[i][j]);
    SD_min += aw*pow(1-fg,2);
    }
//printf("\nSD_min %1.3f %1.3f\n\n",SocialDist[b],SD_min);
}

/*---------------------------------------------------------*/
/*  calculates daily mobility fluctuation based on
                              weekly and seasonal cycles   */
/*---------------------------------------------------------*/
double calc_mobfluc(double t)
{
double wphase,fluc,dayl,ff,fg,temp,tempf,temp_neutral;
double daylength(double,double );

t-=dat_off;  // only needed for Iran due to preceding epidemic start
/*----------------------------------------------------*/
/*     change in mobility demand during the week      */
/*----------------------------------------------------*/
wphase=(t-week_lag)/7;  // weekly phase==4
wphase-=floor(wphase);
if(week_ampl<0.01)  wphase=0.5;

fluc = (1+(2*week_ampl*wphase-week_ampl));
/*-------------------------------------------*/
/*    add seasonal component:  daylength     */
/*-------------------------------------------*/
dayl = daylength(t+41,latitude[DataCountryIndex]);
ff = (dayl-0.5)/(1*(0+dayl_max)-0.5);
sunf = sun_ampl*atan(2*ff);

/*--------------------------------------------------*/
/*   outdoor activity dependent on thermal comfort  */
/*--------------------------------------------------*/
fg = sin((t+0+1.5*dat_off)/365.25*3.1415);
temp = temp_min[DataCountryIndex]+(temp_max[DataCountryIndex]-temp_min[DataCountryIndex])*fg*fg;
/*-----------------------------------------------*/
/*   partial adaptation of neutral temperatures  */
/*-----------------------------------------------*/
temp_neutral=20+0.*(temp_max[DataCountryIndex]-20);

/*------------------------------------------*/
/*   Gaussian function for thermal comfort  */
/*------------------------------------------*/
ftemp = exp(-0.5*pow((temp-temp_neutral)/temp_width,2))-exp(-1);
//ftemp *= temp_ampl*(1+0.3*(dat_off>0));
ftemp *= temp_ampl;

// total seasonality from both proxies
season = fmode*sunf + ftemp;

/*--------------------------------*/
/*   check for consistent value   */
/*--------------------------------*/
if(1+season<0) season=-1;

fluc *= 1+season;
return fluc;
}

/*-------------------------------------------------------------*/
/*   set initial fields for contact rates and age structure    */
/*-------------------------------------------------------------*/
void reset_contactage(int da, int da2)
{
int i,j,b=0,nc;
double fg,ff;

if(da2<0) da2=da;
/* ---------------------------------------------------- */
/*     country specific initial transmission rates
                            based on initial slope      */
/* ---------------------------------------------------- */
fg=1;
for(i=0;i<NUM_IA;i++)
  for(j=i;j<NUM_IA;j++)
     transM[i][j][b] =trans[j][i] = fg*bet[j][i][da];

for(i=0;i<NUM_IA;i++) age_dist[i]=adist[i][da];

// TODO return check of reciprocity
/*----------------------------------------------------*/
/*   household contribution estimated from Prem2017   */
/*----------------------------------------------------*/
ff=1*1.5/houseold[da+1];
//printf("house %1.2f %d %1.1f\n",ff,da,houseold[da+1]);
for(i=0;i<NUM_IA;i++) // 3rd loop over target age classes: RHS of ODEs
  for(j=i;j<NUM_IA;j++) // loop over interacting age classes
    {
    cmob[i][j]=0;
    if(i==j)
      {
      if(i==0) cmob[i][j]=0.35*ff+houseold[da+1]*0.0;
      else
        if(i<5)
          cmob[i][j]=0.5*ff+houseold[da+1]*0.0;
        else
          cmob[i][j]=0.7*ff+houseold[da+1]*0.0;
      }
    if((i==0)&&(j==2)) cmob[i][j]=0.75*ff;
    if(((i==0)&&(j==5))||((i==2)&&(j==5))||((i==2)&&(j==6))||((i==3)&&(j==6)))
       cmob[i][j]=houseold[da+1]*0.02;
    }
}

/*-----------------------------------------------------*/
/*  set initial fields for infection absent/present    */
/*-----------------------------------------------------*/
void init_spread(double InitialInfect)
{
int d,n,b=0,j,i,iz,fi,di,dd;
double ff,fg,dage=20,betall,betall2,avgage,bm,bm0;
double calc_mobfluc(double );

// weekly fluctuation; start with average
weekf = 0.333*(calc_mobfluc(0)+calc_mobfluc(2)+calc_mobfluc(4));
//printf("init_spread weekf=%1.3f\n",weekf);
for(i=avgage=0;i<NUM_IA;i++)
  avgage += ages[i]*age_dist[i];

if (InitialInfect<1E-6)
  {
  for(i=betall=betall2=0;i<NUM_IA;i++)
     {
     FConc[i][b] = age_dist[i];  // susceptibles from age distribution
     if(InitialInfect>-EPS) Conc[i][b]=0; // initialisation of infected age distribution
  /* ---------------------------------------------- */
  /*     IFR depending on age     */
  /* ---------------------------------------------- */
     ff = pow(q10_fatal,(ages[i]-ages[NUM_IA-1])/10);
     fatal[i] = mort*ff;   // relate Q10 to oldest age-class
     // increase mortality in Iran due to less medical infrastructure
     // if(DataCountryIndex==4) fatal[i]*=1.;//excessm[DataCountryIndex];

  /* ---------------------------------------------- */
  /*     renormalize initial transmission rates     */
  /* ---------------------------------------------- */
     fg = avgage-ages[i];
     ff = exp(-fg*fg/(2*dage*dage));

     bm0 = beta_max*attack[i];
     bm  = weekf*bm0;
     for(j=i;j<NUM_IA;j++)
       {
       transM[i][j][b] = trans[j][i]*bm*attack[j];
       transM0[i][j][b] = trans[j][i]*bm0*attack[j];
       betall += transM[i][j][b];
       betall2+= transM[i][j][b]*ff;//   *age_dist[i];
       age_rat[i][j][b] = age_dist[i]/age_dist[j];
       }
     }
   if(sim_out && 0) printf("weekf=%1.3f\n bet all=%1.3f %1.3f\t%1.4f %1.4f %1.4f\n",weekf,betall*1E3,betall2*1E3,beta_max*attack[0],trans[0][0],trans[6][6]);
   }
else
  /* ------------------------------------------------------------- */
  /* initialize infected age distribution as a normal distribution */
  /* ------------------------------------------------------------- */
  {
  if (init_age>0)
    for(i=0,SumC=0;i<NUM_IA;i++)
      {
     //   InfectAge[b]+=Conc[i][b]*ages[i];
      fg=init_age-ages[i];
      Conc[i][b] = exp(-fg*fg/(2*dage*dage));
      SumC+=Conc[i][b];
    //   printf("%1.3f %1.3f\t",fg,exp(-fg*fg/(2*dage*dage)));age_dist[i]
      }
  else
    for(i=0,SumC=1;i<NUM_IA;i++)
      Conc[i][b] = age_dist[i];
  for(i=0;i<NUM_IA;i++)
     Conc[i][b]*=InitialInfect;
  }
}
