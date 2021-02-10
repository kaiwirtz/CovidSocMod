#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "values.h"
#include <ctype.h>
#include "defines.h"
#ifdef __C_PRE
 void progress_log(void);
 void progress_var_log(long );
 typedef char *String;
#else
  extern "C" { int simulation_shell(void);
               int store_prep_open(char*);
 	        void struct_model(void);
 }
#endif
#include "variables.h"

typedef char name[16];
typedef char string[120];
typedef char long_name[80];

/* ---------------------------------------------- */
/*         Simulation control parameters          */
/* ---------------------------------------------- */
extern double startim,endtim,dat_off,start_of_int,end_of_int,cc,delt_of_int,sim_time,delt,mob_scaled,fac_syn,fatal_min,fatal_max,locktarget,dayl_max,fmode,sub_0,norm_dat[NUM_FIT_VARS],var2[NUM_FIT_VARS],**valid_ptr,*dptr,err[N_VARIAT][NUM_DPV],err_f[NUM_FIT_VARS],*other_ptr[NUM_PD],
  *par_val[N_VARIAT],preruntime,weekf,InitFatal,InitFatal0,mrate_max,fac_exp,mobil0,expo,attackfac,fsoc,mfatal,Immune[10],season,ftemp,catchupf,flex,SD_min,fatalavg_sim,fatalavg_dat,Susc_old;
extern double Temp,simavg[NUM_FIT_VARS],simvar[NUM_FIT_VARS],errv[4],par[8],tweight;
extern unsigned int int_time,day_of_year,PeriodNum,SteadyState,IsTime,IsAware,mobn,ncountry,mrim,mrim0,accs;
extern long step_of_int,step_event,stepsum_sim,pdf_out[NUM_EXP];
extern unsigned long store_step,num_total_variat,nvar,d_var,step_var,step_out,step_st,data_step;
extern double CTemp[MAX_VAR][MAX_BOX],CCC[MAX_VAR][MAX_BOX],STemp[MAX_VAR][MAX_BOX],
  *SCCC[MAX_VAR],CCC_i[MAX_VAR][MAX_BOX],CCC_1[MAX_VAR][MAX_BOX];
extern unsigned long l0,l1,step_total,store_step_next,target_i[4];
extern double SumC,Rtrans[NUM_IA],time_Hcap,fconfid0,mobil;
extern int ind_old[MAX_FORCE+3],ind_new[MAX_FORCE+3],ind_end[MAX_FORCE+3],ForcedLock,postNew;
extern int sortind[NUM_FIT_VARS][N_OPTSIM],sortlen[NUM_FIT_VARS],nclass,data_off;
extern double *forcing[MAX_FORCE+3],*force_time[MAX_FORCE+3],sorterr[NUM_FIT_VARS][N_OPTSIM],*mmrates,mmrate,mrate,mtran;
extern float bet[NUM_IA][NUM_IA][20],cmob[NUM_IA][NUM_IA],adist[NUM_IA][20],lockday[20][2],rgrinit[20];

 /* ------------------------------------------------- */
/*     Total rates of change in state variables      */
/* ------------------------------------------------- */
extern double SConc[NUM_IA][MAX_BOX],Strans[NUM_IA][NUM_IA][MAX_BOX],SaRtrans[NUM_IA][MAX_BOX],Sexposition[MAX_BOX],SCostAcc[MAX_BOX],SocialDist0[MAX_BOX];
extern double *transM[NUM_IA][NUM_IA],*transO[NUM_IA][NUM_IA],transM0[NUM_IA][NUM_IA][MAX_BOX],contact[NUM_IA][NUM_IA][MAX_BOX],trans_base[NUM_IA][NUM_IA][MAX_BOX],fac[MAX_BOX],mixr[NUM_IA][NUM_IA],fatal[NUM_IA];

/* -------------------------------------- */
/*           Auxiliary variables          */
/* -------------------------------------- */
extern double TotConc[MAX_BOX],dyncoup,MeanTrans[MAX_BOX],TotMort[MAX_BOX],Mtrans[NUM_IA][MAX_BOX],lim[5][MAX_BOX],fatalrate[MAX_BOX],SocialCosts[MAX_BOX],TotLoss[MAX_BOX],Mobility[MAX_BOX],FConc[NUM_IA][MAX_BOX],InfectAge[MAX_BOX],DeathAge[MAX_BOX],ages[NUM_IA],rate[NUM_IA],Susc[NUM_IA],age_rat[NUM_IA][NUM_IA][MAX_BOX],Contacts[MAX_BOX],TotSusc[MAX_BOX],TotRecov[MAX_BOX],SocialCostAcc[MAX_BOX],IndCosts[MAX_BOX],MortAcc[MAX_BOX][4],ExpoEff[MAX_BOX],Exposure[MAX_BOX],SocialCostAvg[MAX_BOX],*SocialCostMeM[MAX_BOX],*MortalityMeM[MAX_BOX],wave_min[4][2],wave_max[4][2],season0,fvacc,sunf,SocialDist[MAX_BOX];
extern int num_stores, num_others,num_variat,n_comp,variat_steps[N_VARIAT],ind_exp,
    sim_out,num_states,store_prec[NUM_PD],svar_in[NUM_FIT_VARS],err_in[NUM_FIT_VARS],ie,bulk[MAX_VAR],lim_i,num_tr;// act_dat
extern float variat_min[N_VARIAT],variat_max[N_VARIAT],variat_err[N_VARIAT],variat_delt[N_VARIAT],variat_delt0[N_VARIAT],*store_vector;
extern double var_out[N_OUTVAR],var_out0[N_OUTVAR],var_out1[N_OUTVAR],**store_ptr;
extern name variat_names[N_VARIAT],other_names[NUM_PD],store_names[NUM_PD],*state_names,ExpName[NUM_EXP];
extern FILE *outfile,*outfile_exp[NUM_EXP];
extern char pres[199],*varvalues[199], aggfile[99],aggfile2[99],VarFile[99];
extern int setup;

/* -------------------------------------- */
/*          function definitions          */
/* -------------------------------------- */

//int strpos(char *,char *)
int up_init(char *,int),write_index(int),index_nut(int), get_store_ptr(int, name);
void struct_model(void),model(void),integrate_model(int,int,int);
double simulation(int),calc_err(void), exp_simulation(int , int);
/* int simulation_shell(void); */
void set_logparvector(long , FILE*),fix_lastparvar(unsigned long ,int );
char *set_parvector(unsigned long, int);
double random2(void);
void sim_prep(int), init_spread(double),write_rgr(int),init_states(int ), preydist_update(int , int );
void write_last_state(void),write_header(FILE *,char*,int , int ), write_allometry(void);
int valid_prep(int ), calc_dev(long );
double daylength(double,double ),add_lockdown_err(int),radiation(double ,double ,double ,double );
long data_prep(int);
unsigned long up_init_vel(int),up_init_tke(int);
double update_env(), nq_min(double , double ), logsize(int), calc_mobfluc(double );
void reset_contactage(int, int ),  nerror(char *, int );
int strpos(char* ,char* );
char* set_parvector(unsigned long, int);
