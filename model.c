/*-------------------------------------------------------*/
/*  Integrated Covid19 Social and Epidemiological Model  */
/*	Age Based Virus Transmission and Regulation 	    	 */
/*  Model equations   		  	                      	 	 */
/*  Author	 Kai Wirtz (kai.wirtz@hzg.de)		 		 		 		 */
/*  last edit 2021/01/22               					 		 		 */
/*-------------------------------------------------------*/
#include "struct.h"
#include "colour.h"
void model(void)
{
int nt,i,j,k,b,in,ni,m,IsCheck;
double ff,fg,f,dayl,socv,wphase,sweight[NUM_IA][NUM_IA],contact_base[NUM_IA][NUM_IA],transAll,fatalAll,beta_min,beta_min0,expo_target,cb,avgfatal,fmobil;
double f0,f1,dC_dbeta,dM_dbeta,dI_dbeta,aweight,cweight,ampl0,import,exfat;
/*------------------------------*/
/*   flag for output on screen  */
/*------------------------------*/
IsCheck = (fmod(sim_time,10*OutputStep)<1E-3*TimeStep)*(fabs(sim_time-Time)<2E-5)*1;
IsTime  = IsCheck*(Time>=0 && Time<-51)*((DataCountryIndex==4)||(DataCountryIndex==-6));

/*-----------------------------------------------------------------------------*/
/*  calculates daily mobility fluctuation based on weekly and seasonal cycles  */
/*-----------------------------------------------------------------------------*/
weekf = calc_mobfluc(sim_time);

/*----------------------------------------------------*/
/*         loop over spatial boxes, if any            */
/*----------------------------------------------------*/
for(b=0;b<MAX_BOX;b++) {

import = Mobility[b]/100*fvacc*relaxs*fmode*(data_off<0);

/* total density of currently infected individuals */
for(i=0,TotConc[b]=0;i<NUM_IA;i++)
    TotConc[b]+=Conc[i][b];
/*-------------------------------------------------*/
/*      Rebound of socio-economic activities       */
/*-------------------------------------------------*/
catchupf =  CostAcc[b];//*pow(age_dist[0]/0.2,2);

// cut-off for numerical reasons
ff=3E4;
if(catchupf>ff) catchupf=ff;
//if(catchupf>1E-2) catchupf=1E-2;
//printf("%1.1f %1.1e > %1.1e\n",Time,catchupf,social*136);

/*---------------------------------------------------*/
/*      density of susceptibles per age class        */
/*---------------------------------------------------*/
for(i=0,TotSusc[b]=0,SocialCosts[b]=SocialDist[b]=0;i<NUM_IA;i++)
  {
  Susc[i] = FConc[i][b]*exp(-aRtrans[i][b]);
  TotSusc[b]+=Susc[i];
    for(j=i;j<NUM_IA;j++) // loop over interacting age classes
      {
  /*-----------------------------------------------------------*/
  /*  initial (optimal) value incorporates a weekly cyclicity  */
  /*-----------------------------------------------------------*/
      contact_base[i][j] = weekf * transM0[i][j][b] * mixr[i][j];

  /*------------------------------------------*/
  /*     contact rates from transmission      */
  /*------------------------------------------*/
      contact[i][j][b] = transM[i][j][b] * mixr[i][j];
  /*-----------------------------------------------------------*/
  /*     social value of contacts accounts for relative
      distribution of interacting age classes and basal rate   */
  /*-----------------------------------------------------------*/
      aweight = contact_base[i][j]*age_dist[i]*age_dist[j];
      sweight[i][j] = aweight * social * (1+ catchupf);

  /* well-being decreases with increasing social distancing */
      ff = pow(1-contact[i][j][b]/contact_base[i][j],2);
      SocialCosts[b] += sweight[i][j]*ff;
      SocialDist[b]  += aweight*ff;
 /*if(fmod(sim_time,1*OutputStep)<1E-3*TimeStep && fabs(sim_time-Time)<2E-5 && Time<5)
printf("%d %d %1.3f %1.3f\t%1.3f*%1.3f ->%1.3f\t%1.0f\n",i,j,contact[i][j][b],contact_base[i][j],sweight[i][j]*1E6,pow(1-fg,2),SocialCosts[b]/(social + catchupf),contact[i][j][b]*age_dist[j]*age_dist[i]*1E4);*/
      }
  }
//SCostAcc[b] = (postNew>=1)*(1+acc_catch*fmax(exposition[b]-0.5,0)*CostAcc[b]);

exfat = 1;//+0.5/(1+exp(10*(1-fatalrate[b]/15E-6)));

/*----------------------------------------------------------------*/
/* Exposure dynamics
    willingness to wear face masks and to keep distancing rules   */
/*----------------------------------------------------------------*/
mob_scaled =  sqrt(SocialDist[b]);// + sqrt(CostAcc[b]/MemTime);
//mob_scaled =  SocialDist[b];

lim[1][b]=mob_scaled;


expo_target= behave*mob_scaled;//*pow(age_dist[NUM_IA-1]/0.1,2);//*(1.-TotSusc[b]);
//expo_target = expo_tar0*(SocialCostAcc[b]/1E-4);//*(1.-TotSusc[b]);
if(expo_target>0.99) expo_target=0.99;

//if(Time>59.0 && Time<296.4 &&IsCheck)  printf("%03.0f C%5.0f L%5.0f M%4.0f SC%4.0f\t%4.2f Acc=%04d  \te=%6.2f %6.2f catr=%6.2f\n",Time,TotConc[0]*1E5,TotLoss[b]*1E5,fatalrate[b]*1E7,SocialCosts[b]*1E7,SocialDist[b]/SD_min,(int)(CostAcc[b]),exposition[b],1-expo_target,catchupf/social);


/* ------------------------------------------------------------ */
/*    reduction in effective contacts (masks) after shutdown    */
/* ------------------------------------------------------------ */
if(1){
 //ff=social/4E-6;
 //ff=soc_bet/(social+soc_bet);//*age_dist[NUM_IA-1]/0.06;
 Sexposition[b] = behave_rate*( 1 - expo_target - exposition[b]);
 //if(DataCountryIndex==4) Sexposition[b] *=2.5;
 //expo = exposition[b]/sqrt(\n);

 //ff = 1- 2*fmode * sunf - 2*ftemp;
 ff = 1 - 2*((0+fmode) * sunf + ftemp);
 //ff = 1- 2*sunf - ftemp;
 if (ff<5E-2) ff=5E-2;
 expo = exposition[b]*ff;
 Exposure[b] =exposition[b];
 }
else
 {
 Sexposition[b] = 0;
 expo = (1-expo_target)*exp(-1*(season));
 Exposure[b] = (1-expo_target);
 }
ExpoEff[b]  = expo;
// printf(" exposition[b]=%1.2f\n",exposition[b]);

/*-----------------------------------------------------*/
/*     calculate auxiliary variables over age classes  */
/*-----------------------------------------------------*/
transAll=fatalAll=MeanTrans[b]=mrate=mfatal=avgfatal=0;
tweight = IndCosts[b]=cweight=0;
for(i=0,Contacts[b]=mtran=0;i<NUM_IA;i++) // loop over target age classes
  {
  /* average fatality */
  avgfatal+=exfat*fatal[i]*Conc[i][b];
  /* average specific and total transmission of each age class */
  for(j=i,Mtrans[i][b] = Rtrans[i] = 0;j<NUM_IA;j++) // loop over source age classes
     {
     Mtrans[i][b]+= transM[i][j][b]*age_rat[i][j][b];  //expo* average transmission, only for output
     Rtrans[i]   += expo*transM[i][j][b]*Conc[j][b]*age_rat[i][j][b]; // effective integral transmission
//if(IsTime&&sim_out==1 &&i==5) printf("%d %d\t%1.1e tr=%1.1e C=%1.1e ar=%1.1e\t%1.1e\n",i,j, expo,transM[i][j][b],Conc[j][b],age_rat[i][j][b],expo*transM[i][j][b]*Conc[j][b]*age_rat[i][j][b]);
     }
  for(j=0;j<i;j++) // do not double count reciprocal cross-age interactions
     {
  /*   average transmission, only for output   */
     Mtrans[i][b]+= transM[j][i][b]; //expo*
  /*  effective integral transmission for target age class i */
     Rtrans[i]   += expo*transM[j][i][b]*Conc[j][b];
//     if(IsTime&&sim_out==1&&i==5) printf("%d %d\t%1.1e tr=%1.1e C=%1.1e\t%1.1e\n",i,j, expo,transM[j][i][b],Conc[j][b],expo*transM[j][i][b]*Conc[j][b]);
     }

  /* gross infection rate of target group i */
  rate[i] = Rtrans[i]*Susc[i]; // -recov*Conc[i][b];
  mrate += rate[i];

  /* average transmission rate of age class i, only for output */
  Mtrans[i][b]/=NUM_IA;
  mtran += Mtrans[i][b]*age_dist[i];
  /*---------------------------------------------------------*/
  /*  case fatality contribution to marginal community gain  */
  /*---------------------------------------------------------*/
  mfatal += exfat*fatal[i]*age_dist[i];

/*---------------------------------------*/
/*   calculate averaged infection rate   */
/*---------------------------------------*/
  for(k=0; k<NUM_IA;k++) // loop over interacting age classes
    if (i<=k)
      {
    //  mrate += ff*expo*transM[i][k][b];
      Contacts[b]+=fmax(contact[i][k][b]-cmob[i][k]*contact_base[i][k],0)*age_dist[k]*age_dist[i];
      }
    else
      {
    //  mrate += ff*expo*transM[k][i][b]*age_rat[k][i][b];
      Contacts[b]+=fmax(contact[k][i][b]-cmob[k][i]*contact_base[k][i],0)*age_dist[i]*age_dist[k];
      }

  } // end loop i

/* --------------------------------------- */
/*         baseline mobility ratio         */
/* --------------------------------------- */
if(mobn==0 && mobil0>0)
    mobil= Contacts[b]/mobil0;
  Mobility[b]=100*mobil;

if(fmod(sim_time,20*OutputStep)<1E-3*TimeStep && fabs(sim_time-Time)<2E-5 && Time>=0 && Time<-450) printf("%1.0f\tweek %s %1.2f %1.2f%s %1.2e\t %1.2f %s%3.0f%s\n",Time,bold,weekf,season,normal,SocialCosts[b]/(social + catchupf),1E2*Contacts[b],bold,1E2*mobil,normal);

/* average relative growth rate of infected */
mrate /= TotConc[b]+1E-8;
avgfatal/= TotConc[b]+1E-8;
mtran  /= TotConc[b]+1E-8;

/*-----------------------------------------------------*/
/*   subtract recovery rate -> relative growth rate    */
/*-----------------------------------------------------*/
mrate -= recov;

/*-------------------------------------------------------------*/
/*   extrapolated change during the epidemic spread
                          to correct for time-lagged recovery  */
/*-------------------------------------------------------------*/
fac_exp = 1;
/*if(Time<locktarget-7)
   fac_exp = 0.01;
else*/
//fg=1./(1+exp(-0.5*(Time-locktarget+7)));
if(postNew==0)
   fac_exp = exp((mmrate+import*mtran)*dyncoup);  // mmrate is time-averaged mrate updated every day (in simulation.c)

//if((fmod(sim_time,10*OutputStep)<1E-3*TimeStep) && (fabs(sim_time-Time)<2E-5) && (Time>9 && Time<=41.5) && DataCountryIndex==4)  printf("%1.1f %d mu=%1.3f %s%1.2e%s %1.3f\n",Time,postNew,mrate,blue,exp(mrate*dyncoup),normal,fac_exp);
if(fac_exp<0.02)fac_exp=0.02;
//fac[b] = fac_exp;
fac[b] = 1./(social * (1 + catchupf));
for(i=0;i<NUM_IA;i++) // 2nd loop over target age classes
  {
  /*----------------------------------------*/
  /*   recovery, here without time delay    */
  /*----------------------------------------*/
  rate[i] -= recov*Conc[i][b]; //*exp(-mrate/recov);

  /*----------------------------------*/
  /*       total fatality rate        */
  /*----------------------------------*/
  fatalAll += exfat*fatal[i]*Rtrans[i]*Susc[i];
  /*---------------------------------------------------------*/
  /*    relative infection growth rate, includes recovery    */
  /*---------------------------------------------------------*/
  transAll += rate[i];
  /*---------------------------------------------------------------*/
  /*   weighted average of overall transmission rate (no recovery) */
  /*---------------------------------------------------------------*/
  MeanTrans[b] += Rtrans[i]*age_dist[i];
  }
/*  average relative increase of   */
transAll     /= TotConc[b]+1E-8;
MeanTrans[b] /= TotConc[b]+1E-8;

/*------------------------------------------------------------------*/
/*    total fatality rate including excess hospitalization effect   */
/*------------------------------------------------------------------*/
fatalrate[b] = fatalAll;
/*  marginal increase in mortality due to excess hospitalization  */
/* catchup depending on M/C ratio */
//ff=fatalrate[b]/SocialCosts[b]*(mrate>0);
//fg=2./(1+exp(parvar1*3*mrate*21));

SCostAcc[b] = (postNew>=3)*(acc_catch*CostAcc[b]*(2.+0.E-5/social)*(1-CostAcc[b]/catchup));
// IsAware*(1+0*SocialDist[b]/SD_min+CostAcc[b]/parvar1);catchup
//if(SCostAcc[b]>50) SCostAcc[b]=50;

//fg = social/1E-6;
// sub_0 = sub_beta;
//ampl0 = (1-weekampl_sub +2*weekampl_sub*wphase);
/*-----------------------------------------------------*/
/*     output of some aggregated variables to screen   */
/*-----------------------------------------------------*/
if(IsTime&&sim_out==1){
printf("%1.0f %d %s%1.3f%s Tot=%1.0f M=%1.0f %1.3f %1.3f week=%1.2f\tmob= %1.1f\n",Time,postNew,bold,fatalrate[b]/SocialCosts[b],normal,TotConc[b]*1E6,fatalrate[b]*1E8,mrate,mrate+import*mtran,weekf,Mobility[b]);

printf("C=%1.0f Dist=%1.0f%s AC=%1.0f%s e=%1.2f/%1.2f S=%1.0f\t%3.0f %sc=%3.0f%s\n",SocialCosts[b]*1E8,SocialDist[b]*1E5,bold,CostAcc[b],normal,exposition[b],expo,TotSusc[b]*1E3,bold,social*1E6,catchupf*1E6,normal);
//for(i=0;i<NUM_IA;i++)printf("%3.0f ",Mtrans[i][b]*1E4*age_dist[i]);
//printf("\n");
//for(i=0;i<NUM_IA;i++)printf("%3.0f ",Conc[i][b]*1E6); printf("\n");
}

/*------------------------------------------------------------------*/
/*  retrospection/expectation:                                      */
/*  expected contribution of case fatality accounting for time-lag  */
/*------------------------------------------------------------------*/
mfatal *= fac_exp;
avgfatal*= fac_exp;
//if((fmod(sim_time,40*OutputStep)<1E-3*TimeStep) && (fabs(sim_time-Time)<2E-5) && (Time>0 && Time<=180))  printf("%1.1f %6.4f %6.4f\n",Time,mfatal,avgfatal);
  //  60*exp(x/5-1)/7
//    if (ff>contact_base) ff=contact_base;sqrt
//lim[0][b]=lim[1][b]=lim[2][b]=0; // diagnostic variables, here RCI
/*-----------------------------------------------------*/
/*    social interactions: loop over all age groups   */
/*-----------------------------------------------------*/
for(i=0;i<NUM_IA;i++) // 3rd loop over target age classes: RHS of ODEs
  {
/*--------------------------------------------*/
/*    total change in infected  population    */
/*--------------------------------------------*/
  SConc[i][b]= rate[i] + import*(Mtrans[i][b]*age_dist[i]);//*(postNew>=1);// fmobil;
  // no sub-individual account; assuming total pop size of 2E7
  if(SConc[i][b]<0 && Conc[i][b]<1E-7) Conc[i][b]=0;

//  if()
/*-------------------------------------------------------------------*/
/*     accumulated infection rate defines density of susceptibles    */
/*-------------------------------------------------------------------*/
  SaRtrans[i][b] = Rtrans[i];

  for(j=i;j<NUM_IA;j++) // loop over interacting age classes
    {
    //cweight += pow(1-fg,2);
    //if(fg<0.5)  IndCosts[b] += aweight*(pow(0.5-fg,2));
    //if(fg<1)  IndCosts[b] += aweight*(pow(1-fg,2));
  /* pressure increases at growing distance from initial (optimal) value */
    cb = contact_base[i][j];
    dC_dbeta = sweight[i][j] * 2*(cb-contact[i][j][b])/(cb*cb);
    dC_dbeta *= mixr[i][j];

  /* -------------------------------------------------------- */
  /*     marginal dependency of infected density on rates     */
  /* -------------------------------------------------------- */
/*-----------------------------------------------------------*/
/* derivative of mortality induced by indirect transmission
                                    to sensitive age classes */
/*-----------------------------------------------------------*/
//    dI_dbeta = age_dist[i]*Susc[j]/age_dist[j];  //TODO: check
//    dI_dbeta += Susc[i]*age_rat[i][j][b];   /* ... obeying reciprocity */
    ff = Conc[i][b] * Susc[j];
    dI_dbeta = ff;   //TODO: check
/*--------------------------------------------------------*/
/*   direct derivative of mortality rate on transmission  */
/*--------------------------------------------------------*/
    dM_dbeta = exfat*fatal[j] * ff;
  /* ... obeying reciprocity */
    ff = Conc[j][b] * Susc[i] / age_rat[i][j][b];
    dI_dbeta += ff;
    dM_dbeta += exfat*fatal[i] * ff;

    dI_dbeta *= 0.5*expo;
    dM_dbeta *= 0.5*expo*pow(fac_exp,0.2);

    //dM_dbeta *= pow(fac_exp,0.5);
/*---------------------------------------*/
/*     minimal maintenance contact       */
/*---------------------------------------*/
    // how many individuals does each person need to meet?
    //  higher in aged and
//     if (sub_0 > ) ff=weekf * transM0[i][j][b];
//     else ff=sub_0;
//     ff = sub_0*weekf * transM0[i][j][b];
    beta_min  =  fmax(sub_0,cmob[i][j]) * weekf * transM0[i][j][b];

/*-------------------------------------------------------------*/
/*   RHS of specific transmission rates                        */
/*     adaptive dynamics: societal regulation in transmission  */
/*-------------------------------------------------------------*/
    Strans[i][j][b]= flexibility*(dC_dbeta - (dM_dbeta + avgfatal*dI_dbeta ) );

    beta_min0 =  beta_min/weekf; // threshold

//    Strans[i][j][b]= flexibility*(dC_dbeta - (dM_dbeta + 0.018*dI_dbeta ) );
    //if ((Time>=110)*(Time<150)) Strans[i][j][b]=0;

    /* dependency of excess hosp death neglected  */

    if((fmod(sim_time,40*OutputStep)<1E-3*TimeStep) && (fabs(sim_time-Time)<2E-5) && (Time>61 && Time<=-192.5) &&((i==5)||(j==5)||(i==j)))
      printf("%d %d\tc=%1.3e\tbase=%1.2e %s%1.2f%s %1.2f  %5.1f\n",i,j,contact[i][j][b],cb,bold,contact[i][j][b]/cb,normal,pow(1-contact[i][j][b]/cb,2),transM[i][j][b]*1E3);//Susc[i]

  /*--------------------------------------------*/
  /*  contact rates near subsistence levels ?   */
  /*--------------------------------------------*/
    f0=transM[i][j][b]+Strans[i][j][b]*TimeStep; //
  f1=Strans[i][j][b];
    if (f0<beta_min0*1.3 && Strans[i][j][b]<0)
      {
     /*--------------------------------------------------*/
     /*  numerical relaxation of negative over-shooting  */
     /*--------------------------------------------------*/
      ff = (transM[i][j][b]-beta_min);
      fg = 0.5;
      if (transM[i][j][b]-fg*ff<0) fg *= exp((transM[i][j][b]-fg*ff)*1E3);
      Strans[i][j][b]=-fg*ff/TimeStep;
      if(Strans[i][j][b]>1E-2 &0)printf("%d %d \t\t\t%1.3e S=%1.3e\tb=%1.2e/%1.2e %1.2f %1.2e\tmin=%1.2e %1.2e\n",i,j,ff,Strans[i][j][b],contact[i][j][b],cb, weekf,transM0[i][j][b],beta_min0,beta_min);//Susc[i]
      }

 /*----------------------------------------------*/
 /*       store selective/adaptive forces        */
 /*----------------------------------------------*/
  //  lim[0][b] += age_dist[i]/(attack[i]*attack[j])*dC_dbeta;
    lim[0][b] = 1+season;

    if(IsTime&&sim_out==-1 &&((i==-j)||(i==5)||(j==5)))
      printf("%d %d\t%s%7.1e/%7.1e %7.1e/%7.1e %5.2f%s%s %3.0f%s%s\t%3.1f-%3.1f-%3.1f\t%s=%1.1e %1.1e %1.1e\t%7.1e\n",i,j,blue,transM[i][j][b],transM0[i][j][b],beta_min0*1.3,beta_min,fmax(sub_0,cmob[i][j]),normal,bold,100*contact[i][j][b]/contact_base[i][j],normal,red,1E8*dC_dbeta,1E8*dM_dbeta,1E8*avgfatal*dI_dbeta,normal,flexibility*(dC_dbeta - (dM_dbeta + avgfatal*dI_dbeta ) ),Strans[i][j][b],f1,transM[i][j][b]+Strans[i][j][b]*20);//\t%1.1e<%1.1e %1.2f %1.2f,f0,beta_min0,age_dist[i]/(attack[i]*attack[j]),f1 beta_min,sweight[i][j] ,contact[i][j][b]
    }
  }

/* --------------------------------------------------------------------*/
/*-------------------------------------------------*/
/*                output to screen                 */
/*-------------------------------------------------*/
fg=100; ff=1E3;
/* if(sim_time>=2*TimeStep && sim_time<2.1*TimeStep)
   printf("weight=%1.2e (%1.5f)\n",tweight,sim_time);*/

if(IsTime&&sim_out==11)
 printf("%1.1f e=%1.3f+%1.3f p:%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\t%s %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f%s %1.3f\n",
   Time,exposition[b],Sexposition[b],fg*Conc[0][b]/age_dist[0],fg*Conc[1][b]/age_dist[1],fg*Conc[2][b]/age_dist[2],fg*Conc[3][b]/age_dist[3],fg*Conc[4][b]/age_dist[4],fg*Conc[5][b]/age_dist[5],fg*Conc[6][b]/age_dist[6],red,ff*transM[0][0][b],ff*transM[1][1][b],ff*transM[5][5][b],ff*transM[6][6][b],ff*transM[0][6][b],ff*transM[2][4][b],normal,SocialCosts[b]);

} // end loop over spatial compartments
/*-----------------------------------------*/
/*       exchange between model boxes       */
/*------------------------------------------*/
//   for(b=0;b<MAX_BOX-1;b++)
//      advection(b, b+1);
} // end main routine
