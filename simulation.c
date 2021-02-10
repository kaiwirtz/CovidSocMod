/*------------------------------------------------------*/
/*  Integrated Covid19 Social and Epidemiological Model */
/*  Main Simulation Routine    			    			    			*/
/*  Author	 Kai Wirtz                			    				*/
/*  2020/09/10	                	    			    				*/
/*------------------------------------------------------*/
#include "struct.h"

double simulation(int da1)
{
long step,di;
int num,j,i,iz,rgr_iz,b=0,mri,mtf,ie,ddt,da,da0,wi,esi,memn,memf,iVacc=NUM_IA-1;
double ff,fl,fg,ft,critime,time_close;
double NoVacc,VaccInit=323; //vaccination starts Dec, 30 2020
FILE *sp;
const char *const normal = "\033[0m",bold = "\033[1m";

/* --------------------------------------- */
/*      flags for simulation control       */
/* --------------------------------------- */
data_off=1;PeriodNum=0; esi=0;postNew=0;IsAware=0;
sun_ampl = temp_ampl;
/* ------------------------------------------ */
/*       clear/set diagnostic variables       */
/* ------------------------------------------ */
mobil0=Time=SocialDist[b]=0;mobn=0;mobil=1;Mobility[b]=100;sunf=1;
mtf=0;fvacc=1.;
SocialCostAcc[b]=MortAcc[b][0]=MortAcc[b][1]=catchupf=TotLoss[b]=TotMort[b]=0;
da=da0=DataCountryIndex-1;/*   reset diagnosis variables */
CostAcc[b]=1;

/*---------------------------------------------------------------*/
/*  set regional characteristics in seasonality, traveling,
                                                and subsistence  */
/*---------------------------------------------------------------*/
dayl_max=daylength(182,latitude[DataCountryIndex]); //daylength at summer solstice
fmode=1;
if (DataCountryIndex>12 || DataCountryIndex==4) //smaller holiday seasonality in US vs. EU
   fmode=travel;

/*   minimal subsistence transmission */
sub_0 = sub_beta/houseold[DataCountryIndex];

init_spread(0);  // set initial infected and transmission
calc_mobfluc(0); // calc seasonality
/* ------------------------------------- */
/*      simulation loop settings         */
/* ------------------------------------- */
sim_prep(0); // prepares simulation run: boundary input, output steps, ...
if(DataActive)  di=data_prep(1); // prepares comparison with data

/* -------------------------------------------------------- */
/*		             Main simulation time loop		            */
/* -------------------------------------------------------- */
for(step=step_st=step_total=mri=0,ddt=0; step<=stepsum_sim+1 ;step++)
  {
  Time=start_of_int=sim_time=TimeStart+step*TimeStep;
  SocialDist0[b]=Mobility[b];
  /* --------------------------------------------- */
  /*    invoke model routine ..                    */
  /* --------------------------------------------- */
  model();
  /* ------------------------------------ */
  /*    vaccination ..                    */
  /* ------------------------------------ */
  // store immunity level at vaccination start
  if(Time>VaccInit && Susc_old<0) Susc_old=Susc[NUM_IA-1],mrate_max=fatalrate[b];
  if(period_vacc>0) // on or off?
    if(Time>VaccInit && iVacc>=0) //start Dec 2020
     {
     FConc[iVacc][b] -= TimeStep/period_vacc;  // percentage per day, oldest group first
     fvacc-=TimeStep/period_vacc;
     if(fvacc<0) fvacc=0;
     NoVacc=(1-vacc_frac)*age_dist[iVacc];

  // decrement age class if cohort is done
     if(FConc[iVacc][b]<NoVacc)
      {
      if((--iVacc)>=0)
         FConc[iVacc][b]+=(FConc[iVacc+1][b]-NoVacc);  // shift leftover to next class
      FConc[iVacc+1][b]=NoVacc;
      }
     }

  /* --------------------------------------------- */
  /*    total mortality & social losses ..          */
  /* --------------------------------------------- */
  TotLoss[b]+=(fatalrate[b]+SocialCosts[b])*TimeStep;
  TotMort[b]+=(fatalrate[b])*TimeStep;

  /* --------------------------------------------------- */
  /*     stores susceptible density at discrete times    */
  /* --------------------------------------------------- */
  if(step==step_event)
     {
     Immune[esi] =(1.-TotSusc[b])*1E2;
     step_event += (long)(213./TimeStep) + (esi>3)*9E9; //7-month step
     esi++;
     }

  /* ------------------------------------------ */
  /*     detects 1st and 2nd fatality peak      */
  /* ------------------------------------------ */
  if(Time>30 && Time<110)
      if(fatalrate[b]>fatal_min) fatal_min=fatalrate[b];
  if(Time>120 )
     if(fatalrate[b]>fatal_max) fatal_max=fatalrate[b];

  /* ------------------------------------------ */
  /*     storing results for later analysis     */
  /* ------------------------------------------ */
  if((step_st++)%store_step==0)
    {
  /* ------------------------------------------------- */
  /*   calculates average ages of infected and deaths  */
  /* ------------------------------------------------- */
    InfectAge[b]=DeathAge[b]=0;
    for(i=0,ff=ft=fl=0;i<NUM_IA;i++)
       {
       InfectAge[b]+=Conc[i][b]*ages[i];
       ff+=Conc[i][b];
       DeathAge[b] +=fatal[i]*Rtrans[i]*Susc[i]*ages[i];
       fl+=fatal[i]*Rtrans[i]*Susc[i];
       }
    InfectAge[b]/=ff;DeathAge[b]/=fl;

/* -------------------------------------------------- */
/*        calculates contact and mobility rates       */
/* --------------------------------------------------- */
  /* --------------------------------------- */
  /*        initial baseline mobility        */
  /* --------------------------------------- */
    if(Time<2.01)
      {
      mobil0+=Contacts[b];
      mobn++;
      }
  if(fabs(Time-2.01)<0.01 && mobn>0)
       mobil0/=mobn,mobn=0;

    /* --------------------------------------------- */
    /*      detects and stores lockdown time         */
    /* --------------------------------------------- */
    if(Mobility[b]<70&& Time>7) IsAware=1;
    if(DataActive && da>=0)  //lockdown date is needed for calibration of tlag
      {
      if(lockday[da][1]<-9E4)  // first day of transition
        if(Mobility[b]<60.+13.*(da==9 || da==3)+6*(da<13 || da==19))
        {
          lockday[da][1]=-Time;
          fac_syn=fac_exp;
          if(da==-3)
             printf("1 lock %1.2f\tfac_exp=%1.2f\n",lockday[da][1],fac_syn);
        }
      if(lockday[da][1]<0 && lockday[da][1]>-9E4 && Time<-lockday[da][1]+6*(1+(da>11)))
        if(Mobility[b]<50.+13*(da==9 || da==3)+6*(da<13|| da==19)  && Time<-lockday[da][1]*7)
        { // 2nd day of transition -> intermediate defines lockday
          lockday[da][1]=(Time-lockday[da][1])/2;
          if(da==-3)
             printf("2 lock %1.2f\n",lockday[da][1]);
        }
      }

  /* --------------------------------------------- */
  /*      running 7-day avg of infection rates     */
  /* --------------------------------------------- */
    if(DataActive==0 || data_off==-1)
      {
      mmrates[mri%mrim]=mrate;
      // average over entire stored field
      //for(i=0;i<j;i++) printf("\t%1.2f",mmrates[(mri-i+mrim)%mrim]);
      //printf("\n");
      for(i=mmrate=0;i<mrim0;i++) mmrate+=mmrates[(mri-i+mrim)%mrim];
      mmrate/=mrim0;
      mri++;
      }
    // density of immune/recovered (incl deaths)
    TotRecov[b] = 1.-TotSusc[b]-TotConc[b];
  /* --------------------------------------- */
  /*     output results to binary file       */
  /* --------------------------------------- */
    if((sim_out>0)&&(Time>=storetim-0.5))
      {
      num=num_stores*n_comp;
      make_store(store_ptr,store_vector,num); /* pass actual state to output vector */
      fwrite(&store_vector[0],sizeof(float),num,outfile);
      }
    }

/* ------------------------------------------------------- */
/*    initial setting of infected density
                            to match threshold fatality    */
/* ------------------------------------------------------- */
  if(step==data_step-1 && DataActive && data_off==1)
    {
    InitFatal=InitFatal0; PeriodNum=1;
    init_spread(0.5);
    model();
/* re-set infected density to match the observed average initial mortality */
    InfectDens=0.5*InitFatal/(fatalrate[0]+1E-9);
    if (sim_out>0 && 0)
       printf("transini InfectDens=%1.2e fatalrate=%1.2e \tInitFatal=%1.2e\n",InfectDens,fatalrate[0],InitFatal);
    init_spread(InfectDens);
/* ----------------------------------------------------------- */
/*      calc threat awareness from preset lockdown timing      */
/* ----------------------------------------------------------- */
    model();

    dyncoup= Time+tlag+(DataCountryIndex==4||DataCountryIndex==6||DataCountryIndex==20)*10;
//    if(MonteCarlo==-1) dyncoup=0;
    if(sim_out>0&&0) printf("%1.3f I=%1.2e S=%1.3f SCV=%1.2e tlag=%1.2f  dyncoup=%1.3f\n",Time,TotConc[b],Susc[b],social,tlag,dyncoup);
    data_off=-1;PeriodNum=0;
    }

/* ------------------------------------------------------ */
/*     calculate deviation from data, here every day      */
/* ------------------------------------------------------ */
  if(step==data_step && DataActive &&di>=0 && DataCountryIndex>0)
   {
   di=calc_dev(di);
   if( di>=0 ) data_step+=(long)((data[di][0]-Time+dat_off)/TimeStep);//+data_off_tim)

/* ----------------------------------------------- */
/*      post-lockdown expectation behavior         */
/* ----------------------------------------------- */
  if(mmrate<0 && postNew==0 && (Time>17 || (DataCountryIndex==4 & Time>10)) )  {
    postNew=1; time_close=Time;
    }
   }
// preset degradation date ?
if(DecStart>0)
   {if(Time>DecStart) postNew=3;}
else  // check for resurrecting infection wave
  if(((mrate+fvacc*relaxs*mtran*(postNew>=1)>0) ) && postNew==1 && Time>time_close+21 )
   if(fatalrate[b]<1.5E-6 || Time>time_close+120) postNew=3;

/* -------------------------------------------- */
/*        numerical integration ....            */
/* -------------------------------------------- */
  integrate_model(method,0,num_states); /* 1:Runge-Kutta  0:Euler-Cauchy  */
  } // end integration loop

if(DataActive)
  return (double)calc_err();
else return 0.0;
}

/*----------------------------------------------------*/
/*  simple integration routine                     		*/
/*   (runge-kutta-2-3 without time step adaptation)   */
/*----------------------------------------------------*/
void integrate_model(int type,int j1,int j2)/* 1:Runge-Kutta  0:Euler-Cauchy  */
{
int i,j,i0,j0,step_of_int,ds;
double del6,del2,f,tim_old=0,relch;

relch=relchange;
/* ------------------------------------------- */
/*    Calculating maximal relative change      */
/* ------------------------------------------- */
sim_time=Time;
delt=TimeStep;
for(i=0,i0=j0=0;i<n_comp;i++)
  for (j=j1;j<j2;j++)
    if( CCC[j][i]+delt*SCCC[j][i]>0)
      if( ( f=fabs(SCCC[j][i])*delt/(fabs(CCC[j][i])+1E-3) ) >relch )
        {
        delt*=relch/(f+EPS);
        i0=i;
        j0=j;
        }
if(delt<mindelt)
     delt=mindelt;

/* ------------------------------------------------ */
/*   Increase resolution
         if time-step exceededs critical threshold  */
/* ------------------------------------------------ */
if(delt<TimeStep-2*EPS )
  step_of_int=(long)(sqrt(5*TimeStep/delt+0.5));
else
  step_of_int=1;
if(step_of_int>99 ) // warning at stiff situations
  printf("%1.3f steps=%ld\t%s at %d reg=%d\t%1.3e\t%1.3e\n",\
    Time,step_of_int,state_names[j0],i0,DataCountryIndex,SCCC[j0][i0],CCC[j0][i0]);

if((sim_out)&&step_of_int>999)
  {
  printf("Warning: int_steps reset from %ld to 999 !\n",step_of_int);
  step_of_int=999;
  }
step_total+=step_of_int;

/* ---------------------------------------------------- */
/*    variable interval length for Runge-Kutta scheme   */
/* ---------------------------------------------------- */
delt=TimeStep/step_of_int;
del2=delt*0.5;
del6=delt/6;
/* ------------------------------------------------------- */
/*    sub-integration routine (loop within one TimeStep)   */
/* ------------------------------------------------------- */
for(ds=0;ds<step_of_int;ds++)
/* --------------------------------- */
/*       simple Euler-Cauchy         */
/* --------------------------------- */
  if(type==0)
    {
    for(i=0;i<n_comp;i++)
      for (j=j1;j<j2;j++)
      	{
      	CCC[j][i]+=SCCC[j][i]*delt;
      	if(CCC[j][i]<0)
      	  CCC[j][i]=0;
      	}
    sim_time+=delt;
    /* ----------------------------------- */
    /*   re-evaluation of change rates     */
    /* ----------------------------------- */
    if(ds<step_of_int-1)
      model();
    }
  else
/* --------------------------------- */
/*        Runge-Kutta 3-2            */
/* --------------------------------- */
    {
  /* ------------------------------------------- */
  /*    storing actual state and calc 1. base    */
  /* ------------------------------------------- */
    for(i=0;i<n_comp;i++)
      for (j=j1;j<j2;j++)
      	{
      	CTemp[j][i]=CCC[j][i];
      	CCC[j][i]+=(STemp[j][i]=SCCC[j][i])*del2;
      	 if(CCC[j][i]<0)
      	  CCC[j][i]=0;
      	}
    sim_time+=del2;
    model();

  /* ----------------------------------------------------- */
  /*   calc 2nd base while preventing shifts below zero    */
  /* ----------------------------------------------------- */
    for(i=0;i<n_comp;i++)
      for (j=j1;j<j2;j++)
	      {
      	f=2*SCCC[j][i]-STemp[j][i];
      	if(f>0 && STemp[j][i]<0 && CCC[j][i]<mindelt)
      	  f=0;
      	CCC[j][i]+=f*delt;
      	if(CCC[j][i]<0)
      	  CCC[j][i]=0;
      	STemp[j][i]+=4*SCCC[j][i];
      	}
      sim_time+=del2;

      model();
  /* ----------------------------------------- */
  /*    add 3rd part SCCC to overall change    */
  /* ----------------------------------------- */
    for(i=0;i<n_comp;i++)
      for (j=j1;j<j2;j++)
      	{
      	CCC[j][i]=CTemp[j][i]+del6*(STemp[j][i]+SCCC[j][i]);
      	if(CCC[j][i]<0)
      	  CCC[j][i]=0;
      	}
    }
}
