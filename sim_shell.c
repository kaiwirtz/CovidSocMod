/*-------------------------------------------------------*/
/*  Integrated Covid19 Social and Epidemiological Model  */
/*    Simulation control incl. variation runs    	       */
/*					defining distinct numerical experiments		   */
/*  Author	 Kai Wirtz (wirtz@hzg.de)	                 	 */
/*  2020/08/11	                                         */
/*-------------------------------------------------------*/
#include "struct.h"

double InSz,VaSz,InFu,MiSz,peak[22],minerr[NUM_FIT_VARS],rerr[NUM_FIT_VARS];
long mindvar[NUM_FIT_VARS];
/* -------------------------------------------------- */
/*    	 Loops for parameter variation		      */
/* -------------------------------------------------- */
int simulation_shell()
{
FILE *sp,*sp2,*sps,*spm;
int d,dd,di,ds,da,da1,da2,nex,stored_sim,in,em_n,j,i,i0,vi,vn,dv,ve,mtargi,np,nv,v1,v2,v3,simrun,mc0;
unsigned int pind[NUM_FIT_VARS][2*N_OPTSIM],pindl[NUM_FIT_VARS];
unsigned long indmin,ssn,stepvarlog=1,ivarlog,dss,nv1=6,nv2=6,nv3=1,nvp,dvi,pind_IE;
double em,em_s,sumd,emi,relch,val,f,ff,fg,err0,par0[7],parv[7],s0,m0,*parp,parsv0;
char indfile[99],sufn[4],*pres,outfn[100],infn[100],cmd[200],*sdot;
char par_names[13][44]={"initial transmission $\\beta_0'$","social value $1/H$","attack rate $\\alpha_1$","attack rate $\\alpha_7$","learning rate $r_b$","responsiveness $\\delta$","age-fatality increase $Q_{10}$","infection fatality $\\omega_7$","recovery $r$"," sensitivity $\\epsilon$","seasonal ampl. $a_y$","external input $\\gamma'$","human value $H$"};
char abbrev[20][3]={"BE","FR","DE","IR","IE","IT","NL","PT","ES","SE","CH","UK","GA","IN","LA","MA","MI","NJ","NY","WA"};

init_states(0);  /* store state variables at first time */

/* -------------------------------------------------- */
/*   number of regions from size of data input file   */
/* -------------------------------------------------- */
if(DataActive)
   {
   ncountry = (ColumnsOfdata-1)/2;
   if(DataCountryIndex >0) di=data_prep(0);
   }
else
   ncountry=0;
da=DataCountryIndex-1;

/* ---------------------------------------------------------- */
/*    retrieving index for for calibrated parameter values    */
/* ---------------------------------------------------------- */
if(DataActive && (AnalysisType>=5))
  {
  printf("opening index file %s\n",VarIndexFile);
  sp = fopen(VarIndexFile,"r");
  for(d=0;d<ncountry;d++)
    {
    fread(&di,sizeof(unsigned int),1,sp);
    if (di>=2*N_OPTSIM) di=2*N_OPTSIM-1;
    pindl[d]=di;
    fread(&pind[d][0],sizeof(unsigned int),di,sp);
    //printf("%d %d\t%ld %ld  ... %d\n",d+1,di,pind[d][0],pind[d][1],pind[d][111]);
    }
  fclose(sp);
  }

/*------------------------------------------------------------*/
/*  load initial fields of contact rates and age structure
     for 20 regions; CAUTION if country list changes  */
/*------------------------------------------------------------*/
if(DataActive || DataCountryIndex !=0 || AnalysisType==8)
  {
  printf("load contact matrices\n");
  sp= fopen("contact_matrices.bin","r");
  fread(&bet,sizeof(bet),1,sp);
  fclose(sp);
  sp= fopen("agedist.bin","r");
  fread(&adist,sizeof(adist),1,sp);
  fclose(sp);

  // clear cumulative mismatch errors
  reset_contactage(da,-1);
  err_f[da]=err_f[da+ncountry]=err_in[da]=err_in[da+ncountry]=0;
  /* if(NumDice>-1) for(d=0;d<NUM_DPV;d++) err[da][d]=0; */ // NOT USED
  }

if(NumDice==1111) /* writes result files during variation */
  sim_out=1;
for(d=0;d<ncountry;d++)   // index of best fitting parameter set; NOT USED
   mindvar[d]=-1;

/* --------------------------------------------- */
/*	      Variation is active ...            */
/* --------------------------------------------- */
if(VarActive && (num_variat>0||MonteCarlo>0) )
  {
/* -------------------------------------------------- */
/*   major control switch for
   	   type of num eperiment/parameter variation		  */
/* -------------------------------------------------- */
  switch(AnalysisType)
   {
   case 1:
 /* ----------------------------------------------- */
 /*                                                 */
 /*         systematic parameter variation          */
 /*                                                 */
 /* ----------------------------------------------- */
   d_var = write_index(1);  /* clears old var files */
 /* --------------------------------------------------- */
 /*    controlling output step of variation results     */
 /* --------------------------------------------------- */
  stepvarlog=VarOutputStep + (VarOutputStep<=0) ;

  if (DataActive)  // not really used ; TODO remove
     for(d=0;d<NUM_FIT_VARS;d++)
        for(minerr[d]=9E9,sortlen[d]=0,dd=0;dd<N_OPTSIM;dd++)
            sorterr[d][dd]=9E9;

   nvp=num_total_variat;
/* ----------------------------------------------------------------------- */
/*    partitions loop into for multiple CPUs, needs extra args at call     */
/* ----------------------------------------------------------------------- */
//nvp = nv1*nv2; /* number of total variations over RandomInit */
   if(nvar>1)
      {
      nv1   = (int)ceil(num_total_variat/nvar);
      d_var = RandomInit*nv1;  // lower bound
      nvp   = (RandomInit+1)*nv1; // upper bound
      if (nvp>num_total_variat)nvp=num_total_variat;
    //  printf("split nvar=%d nv1=%d nvp=%d\n",nvar,nv1,nvp);
      }
   for(ivarlog=0,em=9E9;d_var<nvp;d_var++)
     {
/* --------------------------------------------------- */
/*     sets current values of parameters	       */
/* --------------------------------------------------- */
     pres=set_parvector(d_var,0);  // set calibrated parameter values
     sim_out=0; // no output during variation

 /* -------------------------------------------------------- */
 /*   loop over experiments as defined by IndExp & remin_s.  */
 /* -------------------------------------------------------- */
     nex=-Experiment; if(nex<1) nex=1;
//      printf("dvar=%d\t nex=%d %d\n",d_var,nex,simrun);
     for(d=0;d<N_OUTVAR;d++) var_out[d]=0; // diagnosis; NOT USED

 /* ----------------------------------------- */
 /*     open result files during variation    */
 /* ----------------------------------------- */
   if (sim_out>0)
       {
       fclose(outfile);
       strcpy(outfn,MovieResFile);
       if (DataActive==0)
         {
         sdot=strchr(outfn,'.');
         sprintf(sdot-1,"%s.outb\0",pres);
         printf("opening in var %s ...\n",outfn);
         }
       store_prep_open(outfn);
       }

 /* ---------------------------------------- */
 /*           loop over all regions          */
 /* ---------------------------------------- */
  if(ncountry>0)
    for(DataCountryIndex=1;DataCountryIndex<=ncountry;DataCountryIndex++)
      {
      da=DataCountryIndex-1;
    //  if(NumDice>-1) for(d=0;d<NUM_DPV;d++) err[da][d]=0;
  /* ---------------------------------------- */
  /*       sets region specific environment
               such as mixing rates or data   */
  /* ---------------------------------------- */
      reset_contactage(da,-1);
      data_prep(0);   // prepares data arrays for comparison
      init_states(1); // initialize state variables
      em=(double)simulation(da+1);  // runs simulation
      //printf("%s\t\033[1m\033[34m %1.2f\033[0m \n",varvalues,err_f[da]);
      }
  else
    {
    init_states(1); // initialize state variables
    em=(double)simulation(0); // runs simulation
    }

  if(sim_out>0) // closes result file of model dynamics
     fclose(outfile);

 /* ---------------------------------------- */
 /*     writes varied par values in file     */
 /* ---------------------------------------- */
  sp=fopen(VarResFile,"a");
  fprintf(sp,"%d %d ",RandomInit,d_var);
  fprintf(sp,"%s\t",varvalues);

 /* ---------------------------------------- */
 /*    writes varied par values in 2nd file  */
 /* ---------------------------------------- */
 /*  if(NumDice>-1)
    {
    sp2=fopen(VarFile,"a");
    fprintf(sp2,"%d %d ",RandomInit,d_var);
    fprintf(sp2,"%s\t",varvalues);
    for(d=0;d<NUM_DPV;d++)
      {
      for(da=0;da<20;da+=1)
        fprintf(sp2,"%1.4f ",err[da][d]);
      fprintf(sp2,"\t");
      }
    fprintf(sp2,"\n");
    fclose(sp2);
  }
  */
 /* ------------------------------------------ */
 /*    writes model skill in variation file    */
 /* ------------------------------------------ */
  if (DataActive)
    {
    for(d=0;d<ncountry;d++)
      {
      fprintf(sp,"%1.4f ",err_f[d]);
      // index of best fitting parameter set; NOT USED
      if (err_f[d]<minerr[d]) minerr[d]=err_f[d], mindvar[d]=d_var;
      } //d<ncountry
    fprintf(sp,"%1.3f\n",em);
    }
  else
    {
    fprintf(sp,"%1.3f\n",SocialCostAcc[0]*1E6);
    }
   fclose(sp);

 /* -------------------------------------- */
 /*    progess of variation in log file    */
 /* -------------------------------------- */
  if((d_var)%stepvarlog==0 && 1)
      {
      progress_var_log(d_var+1);
      printf(".. %ld/%ld %ld\n",d_var,nvp,num_total_variat);
      }
   } /* for d_var */
    break;
  case 6:
 /* ------------------------------------------------------ */
 /*     reruns optimal parameter-set for each region       */
 /* ------------------------------------------------------ */
   sim_out=1;TimeEnd=450;TimeEnd=690;
   // construct names of output files
   strcpy(indfile,"extra_mrate_");
   sdot=strstr(MovieResFile,"Res");
   printf("indfile %s %s %s\n",indfile,sdot,sdot+3);
   if(sdot!=NULL) strncat(indfile,sdot+3,2);
   strcat(indfile,".dat");
   printf("%s\n",indfile);
   // open output files
   sp=fopen(indfile,"w"); /* fprintf(sp,"Indices\t");*/
   sp2=fopen("parvaluesC.txt","w"); // for extended calibration
 /* ------------------------------------------------------- */
 /*     evaluate optimal parameter-set for each region      */
 /* ------------------------------------------------------- */
   for(d=0;d<ncountry;d++)
      {
      for(dd=0;dd<(MonteCarlo==0)+(MonteCarlo!=0)*pindl[d];dd++)
        {
        dvi=pind[d][dd]+1; // index for optimal parameter set
        pres=set_parvector(dvi,1); // set optimal parameter values
        // extended calibration results
        fprintf(sp2,"%s:%1.1f/%1.0f/%1.0f, ",abbrev[d],acc_catch*1E2,relaxs*1E4,DecStart+41);
        par0[0]=relaxs;
    /* ---------------------------------------- */
    /*        loop over import scenarios        */
    /* ---------------------------------------- */
        for(di=0;di<1+2*(MonteCarlo!=0);di++)
          {
          relaxs = (di+1-3*(di==2))*par0[0]; // import higher/off
          printf("%d/%d\033[1m %s\033[0m\t%d\t%s\trelaxs=%1.2e\n",dd,pindl[d],HeadOfdata[d+1],dvi,pres,relaxs);
          // construct names of output files
          strcpy(outfn,MovieResFile);
          sdot=strchr(outfn,'.');
          sprintf(sdot-1,"opt%d_%d.outb\0",d,dd*1+di+1);
          //printf("opening %s ...\n",outfn);
          store_prep_open(outfn);
          DataCountryIndex=d+1;
          reset_contactage(d,-1);
          data_prep(0); // prepares data arrays for comparison
          init_states(1); // initialize state variables
          em=(double)simulation(d+1); // runs simulation
          printf("simulation(s)\033[1m\033[34m %d %s finished (w error %1.4f) \033[0m ..\n",d,HeadOfdata[d+1],em);
          fclose(outfile);
          if(dd==0&& di==0) rerr[d]=em;
          if(dd==0 && di==0) fprintf(sp,"%1.2e %1.2f %1.3f %1.1f %1.1f %1.1f %1.3e %1.3e %1.3f %1.3f %1.3f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f\n",1./social,tlag,beta_max,med_ages[1+d],temp_min[1+d],temp_max[1+d],fatal_min,fatal_max,dyncoup,fac_syn,houseold[1+d],latitude[1+d],Immune[0],Immune[1],Immune[2],Susc_old*100,mrate_max*1E6);
          }
        relaxs=par0[0];
        }
      }
    fclose(sp);
    fclose(sp2);
    break;

   case 8:
    /* ---------------------------------------- */
    /*      systematic parameter variation      */
    /* ---------------------------------------- */
     sim_out=1;
     par0[0]=beta_max;
     for(di=0;di<Experiment;di++)// loop over regions
        {
        d=IndExp[di]-1;    // region index from list in control file
        dvi=pind[d][0]+1;  // index for optimal parameter set
        pres=set_parvector(dvi,1); // set calibrated parameter values
        printf("%d/%d \033[1m %s\033[0m\t%d\t%s\n",di,Experiment,HeadOfdata[d+1],dvi,pres);
        par0[1]=sub_beta;
        nv=11; // total number of variations
      /* ---------------------------------- */
      /*        loop over parameters        */
      /* ---------------------------------- */
        for(dd=0;dd<=nv;dd++)
          {
          switch(dd) {
            case 0:
               parp=&beta_max; break;
            case 1:
                parp=&social; break;
            case 2:
                parp=&attack[0]; break;
            case 3:
               parp=&attack[6]; break;
            case 4:
               parp=&behave_rate; break;
            case 5:
               parp=&flexibility; break;
            case 6:
               parp=&q10_fatal; break;
            case 7:
               parp=&mort; break;
            case 8:
               parp=&recov; break;
            case 9:
               parp=&behave; break;
            case 10:
               parp=&temp_ampl; break;
            case 11:
               parp=&relaxs; break;
             }
          //printf("%d \t %d %d %d %d %d\n",dd,(dd&(1<<0))>0,(dd&(1<<1))>0,(dd&(1<<2))>0,(dd&(1<<3))>0,(dd&(1<<4))>0);
          // store reference parameter value
          parsv0=*parp;
          ff=floor(log10(parsv0));
      /* ---------------------------------- */
      /*       formatted latex output       */
      /* ---------------------------------- */
          if(di==0)
           if(fabs(ff)<2)
            printf("\\put(\\the\\numexpr \\xs+%d*\\dx,\\the\\numexpr \\ys+%d*\\dy)\{\\footnotesize %s = %1.2f$\\pm$%1.2f\}\n",dd%2,(int)floor(dd*0.5),par_names[dd],parsv0,parsv0*(0.5-(dd==0)*0.3));
           else
            printf("\\put(\\the\\numexpr \\xs+%d*\\dx,\\the\\numexpr \\ys+%d*\\dy)\{\\footnotesize %s = %1.1f$\\pm$%1.1f 10$^{%1.0f}$\}\n",dd%2,(int)floor(dd*0.5),par_names[dd],parsv0/pow(10,ff),parsv0/pow(10,ff)*(0.5-(dd==0)*0.3),ff);
      /* ----------------------------------------- */
      /*        loop over 3 parameter values       */
      /* ----------------------------------------- */
          for(i=0;i<3;i++)
            {
            da2=-1;
            // ref, 50% up, 50% down
            *parp=parsv0*(1+(i-1)*(0.5-(dd==0)*0.3));
            //printf("%d %d\t%ld\t%1.2f %1.2f\n",d,dd,parp,parsv0,*parp);

            //  name & open simulation output file for each variation
            strcpy(outfn,MovieResFile);
            sdot=strchr(outfn,'.');
            sprintf(sdot-1,"sens%d_%d_%d.outb\0",d,dd+1,i);
            //printf("opening %s ...\n",outfn);
            store_prep_open(outfn);

            DataCountryIndex=d+1; // sets region index
            reset_contactage(d,-1);
            data_prep(0); // prepares data arrays for comparison
            init_states(1); // initialize state variables
            em=(double)simulation(d+1); // runs simulation
            if(di>0) printf("simulation(s)\033[1m\033[34m %d %d finished (w error %1.4f) \033[0m ..\n",d,i,em);
            fclose(outfile);
            }
          *parp=parsv0;  // resets original par value
          d=IndExp[di]-1;
          } // loop over parameters
        } // loop over regions
   break;
  } /* switch 0: */
 }
else
/* --------------------------------------------------- */
/*                                                     */
/*  no variation -> 1 single simulation run.....       */
/*                                                     */
/* --------------------------------------------------- */
  if(AnalysisType==0)
   {
 /* ------------------------------------------------------ */
 /*    systematic par variation, but no calibration        */
 /* ------------------------------------------------------ */
 /* --------------------------------------------------- */
 /*    controlling output step of variation results     */
 /* --------------------------------------------------- */
   for(d_var=0,ivarlog=0,em=9E9;d_var<num_total_variat;d_var++)
     {
/* --------------------------------------------- */
/*     sets current values of parameters	       */
/* --------------------------------------------- */
     pres=set_parvector(d_var,1);  // set calibrated parameter values
/*    set_logparvector(d_var,sp);*/
 /* -------------------------------------------------------- */
 /*   loop over experiments as defined by IndExp & remin_s.  */
 /* -------------------------------------------------------- */
     nex=-Experiment; if(nex<1) nex=1;
     for(ind_exp=0,em_n=em_s=0;ind_exp<nex;ind_exp++,em_n++)
       {
       init_states(1); // initialize state variables
       em=(double)simulation(ind_exp); // runs simulation
       }
     }
   }
 else  // AnalysisType>0
  {
/* ------------------------------------------------- */
/*     variation of single scenario parameters	     */
/* ------------------------------------------------- */
  printf("starting %d single simulation(s) .. \n",1-(Experiment<-1)*(Experiment+1));
  da=DataCountryIndex-1;
  if(Experiment<0)
    // loop over regions
    for(DataCountryIndex=da+1+(da<0);DataCountryIndex<=ncountry;DataCountryIndex+=1+(da>=0)*ncountry)
      {
/* ----------------------------------------------- */
/*      open specific variation result files	     */
/* ----------------------------------------------- */
      di=strlen(VarResFile);
      VarResFile[di-6]=48+(int)((DataCountryIndex)*0.1);
      VarResFile[di-5]=48+DataCountryIndex%10;
      //printf("opening %s \n",VarResFile);
      sp=fopen(VarResFile,"w");
      // set calibrated parameter values
      pres=set_parvector(pind[DataCountryIndex-1][0]+1,1);

  /* ----------------------------------- */
  /*     which scenario parameter ?	     */
  /* ----------------------------------- */
     switch(MonteCarlo) {
        case -1:
          parp=&flexibility; break;// varies flexibility
        case -2:
//          parp=&social;  break;// varies mitigation readiness H
          parp=&social;
          par0[1]=social;par0[2]=relaxs;par0[3]=acc_catch; break;
        case -6:
          parp=&period_vacc;  break;// varies vaccination speed
        }
      par0[0]=*(parp);
      di=data_prep(0);
      strcpy(outfn,MovieResFile);
      sdot=strchr(outfn,'.');
      sprintf(sdot-1,"_%d__.outb\0",DataCountryIndex);

      store_prep_open(outfn);
      reset_contactage(DataCountryIndex-1,-1);
/* ---------------------------------------------- */
/*   loop over experiments as defined by IndExp .  */
/* ----------------------------------------------- */
      for(ind_exp=0,em_n=em_s=0;ind_exp<-Experiment;ind_exp++,em_n++)
        {
  /* ---------------------------------------------- */
  /*     prepares & runs simulation .....	          */
  /* ---------------------------------------------- */
        init_states(1); // initialize state variables
  /* ------------------------------------ */
  /*   set values of scenario parameter   */
  /* ------------------------------------ */
        switch(MonteCarlo) {
          case -1:
            *(parp)=exp(ind_exp*0.4-2),relaxs=acc_catch=0; break;// varies flexibility
          case -2:
             *(parp)=par0[1]*exp(-ind_exp*0.13);
             relaxs=acc_catch=0;
             break;// varies SCV catchup=0
          case -6:
             *(parp)=(double)IndExp[ind_exp];  break;// varies vaccination speed
          }

        outfile=outfile_exp[ind_exp];
        // runs simulation
        em=(double)simulation(DataCountryIndex);

        printf("\033[1m\033[34m experiment %d %d\033[0m catc=%1.2f %5.1f \t%7.1f %7.1f\n", ind_exp,DataCountryIndex,acc_catch,*(parp),TotMort[0]*1E6,TotLoss[0]*1E6);
        fclose(outfile);//,ExpName[ind_exp]
    /* ------------------------------ */
    /*   output into variation file   */
    /* ------------------------------ */
        fprintf(sp,"%1.2e %1.2f %1.2f\n",*(parp),TotMort[0]*1E6,TotLoss[0]*1E6);
        *(parp)=par0[0];  // resets original par value
        }
      relaxs=par0[2]; social=par0[1];acc_catch=par0[3]; // resets original par value

      fclose(sp);
      }
  else // Experiment>=0
/* ----------------------------------------- */
/*       single simple / no experiments      */
/* ----------------------------------------- */
      {
      init_states(1);
      em_n=1;
      if(DataCountryIndex>0)
        {
        // runs single simulation
        em_s=(double)simulation(da1);
        }
      else // loop over regions
        for(DataCountryIndex=1;DataCountryIndex<=ncountry;DataCountryIndex++)
         {
          strcpy(outfn,MovieResFile);
          sdot=strchr(outfn,'.');
          da=DataCountryIndex-1;
          sprintf(sdot-1,"%d.outb\0",da);
          //printf("opening %s ...\n",outfn);
          store_prep_open(outfn);
         reset_contactage(da,-1);
        data_prep(0);init_states(1);
        em=(double)simulation(da+1);
      }
      if(sim_out)
         fclose(outfile);
      }
  write_last_state();
  if(DataCountryIndex>0) ff=err_f[DataCountryIndex-1];
  else ff=0;
  printf("simulation(s) %d finished (w error %1.4f %1.4f %d) L=%1.2f..\n",DataCountryIndex,em_s/em_n,ff,TotLoss[0]*1E6);
/*  printf("simerr\t\t%1.3f\n",em);*/
  }
}

int write_index(int new)
{
int di,dd,d,dn;
char indfile[99],lin[399];
FILE *sp,*spi;
/* ---------------------------------------------- */
/*    denoting and open variation result files    */
/* ---------------------------------------------- */
di=strlen(VarResFile);
VarResFile[di-7]=48+(int)((RandomInit)*0.01);
VarResFile[di-6]=48+(int)((RandomInit%100)*0.1);
VarResFile[di-5]=48+RandomInit%10;
strcpy(indfile,VarResFile);
strcat(indfile,".ind");
//printf("%s %s\n",VarResFile,indfile);
strcpy(aggfile,VarResFile);
strncpy(&aggfile[di-3],"res",3);
//printf("%s %d\n",aggfile,new);

/* --------------------------------------------- */
/*     write formatted header of index file      */
/* --------------------------------------------- */
switch (new)
  {
 case 1:
  sp=fopen(VarResFile,"w"); fclose(sp);
  printf("%s %d\n",VarResFile,NumDice);
  if(NumDice>-1)
   //for(d=0;d<NUM_DPV;d++)
     {
     strcpy(VarFile,VarResFile);
     printf("%s %d\n",VarFile,strlen(VarFile));

     strcat(VarFile,".dyn");
     printf("cleaning %s ... \n",VarFile);
  //sprintf(outfn,'%s\0',output_file);
     sp=fopen(VarFile,"w"); fclose(sp);
     }
  sp=fopen(aggfile,"w"); fclose(sp);
  printf("cleaning %s %s ... \t",VarResFile,aggfile);
  printf("writing %s ... \n",indfile);
  sp=fopen(indfile,"w"); /* fprintf(sp,"Indices\t");*/
  for(dd=0;dd<num_variat;dd++)
    fprintf(sp,"%s\t\t",variat_names[dd]);
  fprintf(sp,"\n\t");
  for(d=0;d<num_variat;d++)
     fprintf(sp,"%1.3e %1.3e\t",variat_min[d],variat_min[d]+variat_delt[d]*(variat_steps[d]-1));
  if(variat_min[d]<0.2*variat_delt[d]*(variat_steps[d]-1) && variat_steps[d]==3 )
    fprintf(sp,"1\t\t");
  for(d=0;d<num_variat;d++)
    if(variat_min[d]<0.2*variat_delt[d]*(variat_steps[d]-1) && variat_steps[d]==3 )
      fprintf(sp,"1\t\t");
    else
      fprintf(sp,"0\t\t");
  fprintf(sp,"\nNtotVariat = %d\n",num_total_variat);
  fclose(sp);
  return 0;
case 0:
  sp=fopen(VarResFile,"r");
  dn=0;
  while(!feof(sp))
    {
    fgets(lin,399,sp);
//   printf("%s\n%d\n",lin,dn);
    dn++;
    }
  fclose(sp);
  printf("%d variation entries found\n",dn);
  if(dn>0) dn--;
  return dn;
  }
}
