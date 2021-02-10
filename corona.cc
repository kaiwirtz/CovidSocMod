/////////////////////////////////////////////////////////////////////////////
//
//  Age Based Virus Transmission and Regulation Model
//
//  Author       Kai Wirtz
//
//  Initial Revision 0.0  2020/04/10
//
////////////////////////////////////////////////////////////////////////////

#include "callSiSi.hh"
#include <stdio.h>
#include "struct.h"

void generate_varlist();
void generate_reslist();
double LeftFromSep (String );
double RightFromSep (String );
double *par_val[N_VARIAT];

int main (int argc, char* argv[])
{
int d;
double *tmpp,err;
char resfile[92],*cptr,input[120];


if( !SiSi::parseSimulation(argc, argv) ) // Read all parameters. <------- !!!
    return 0;                            // Error!!!

nvar=1;
if(argc>3)
  {RandomInit=atoi(argv[2]); nvar=atoi(argv[3]);}
//  printf("args: %s %s\n",argv[2]);
printf("RI: %d %d\n",RandomInit,nvar);

/*--------------------------------------------*/
/* speichern der C++ liste "OutputVariables"  */
/*--------------------------------------------*/
num_stores=num_variat=0;
TwoWayList* el = NULL;
  //  TwoWayListElement* ele = NULL;
el = SiSi::OutputVariables;
ResultParameter* el2 = (ResultParameter*) el->resetIterator();
while( el2 )
  {
  if(  el2->isActive() )
    { 	//MessageHandler::information(el2->getName());
    strcpy(input,el2->getName());
    strncpy(store_names[num_stores++],input,18);
    store_prec[num_stores-1] = (int) el2->getPrecision();
  //    printf("%s stored w %d\n",store_names[num_stores-1],strlen(input));
    }
  el2 = (ResultParameter*) el->nextElement();
  }
/*----------------------------------------------*/
/* speichern der C++ liste "VariationVariables" */
/*----------------------------------------------*/
num_total_variat=1;
el = SiSi::VariationVariables;
el2 = (ResultParameter*) el->resetIterator();
while( el2 )
  {
  if(  el2->isActive() )
    { 	//MessageHandler::information(el2->getName());
//       strcpy(store_names[num_stores++],el2->getName());
    strcpy(variat_names[num_variat],el2->getName());
/*--------------------------------------------*/
/* Bestimmung der randwerte fuer variation    */
/*--------------------------------------------*/

/*---------------------------------------------------*/
/*  Identifikation aus der Parserliste anhand Name   */
/*---------------------------------------------------*/
  for(d=0;d<num_variat_parser;d++)
    if(strcmp(VAR_NAMES[d],variat_names[num_variat])==0)
      {
      par_val[num_variat]=VAR_VAL[d];
      //      cout << num_variat<<" -> "<<d<<"\n";
      d=num_variat_parser+3;
      }
  if(d<num_variat_parser+2)
      cout << variat_names[num_variat]<<"/"<< VAR_NAMES[d]<<" not found in callSiSi.cc !!\n";

   /*--------------------------------------------*/
   /*    calc min, number of steps and
                  intervallength of variation    */
   /*--------------------------------------------*/
      variat_min[num_variat]=LeftFromSep(el2->getRange());
      double max = RightFromSep(el2->getRange());
      variat_max[num_variat] = max;
      variat_steps[num_variat] = (int) el2->getPrecision();
//        cout << variat_names[num_variat] << " :" << variat_steps[num_variat];
//        cout << " max=" << max << " min=" <<  variat_min[num_variat];
//        cout << " addr_1=" << VAR_VAL[num_variat]<<" ";
//        cout<<END_OF_LINE;

      if(max >variat_min[num_variat] && variat_steps[num_variat] > 1)
        {
        variat_delt[num_variat]=(max-variat_min[num_variat])/
                                  (variat_steps[num_variat]-1);
//      cout<<   variat_delt[num_variat]<<" "<<variat_steps[num_variat]-1;
//      cout<<END_OF_LINE;
   /*---------------------------------------------------------*/
   /*    calc total number of runs in systematic variation    */
   /*---------------------------------------------------------*/
        num_total_variat*=variat_steps[num_variat++];
        }
      else
        {
   //      variat_min[num_variat]=*VAR_VAL[num_variat];
        variat_delt[num_variat++]=0;
        }

      }
   el2 = (ResultParameter*) el->nextElement();
   }

/*-------------------------------------------------*/
/*  setup definition according to name  0: chemo 1:mesocosm 2:reede.spring */
/*-------------------------------------------------*/
setup=-1;
if(strncmp(SiSi::SimulationName,"UK",9)==0)
  setup=3;

if(num_variat!=num_variat_parser)
  {
  cout << num_variat << " ACTIVE variation_pars found, ";
  cout << num_variat_parser;
  cout << " at maximum expected from parsing files...\n";
  }

/*-------------------------------------------------*/
/*  display some output during single simulations  */
/*-------------------------------------------------*/
if((VarActive&&AnalysisType<10) ) sim_out=0; //&& NumDice>0
else  sim_out=1;
if(AnalysisType==7 && (-Experiment>LengthOfIndExp)) sim_out=0;
/*---------------------------------------*/
/*     set global variables and fields   */
/*---------------------------------------*/
 struct_model();
/*---------------------------------------------*/
/* vorbereitung der ausgabe in ergebnisdatei   */
/*---------------------------------------------*/
printf("sim_out=%d\n",sim_out);
if(sim_out)
  if(store_prep_open(MovieResFile)==1)
    {
    printf("Some error occurred in opening file %s !\n",resfile);
    cout << "Simulation stops.\n";
    abort();
    }
/*------------------------------------------*/
/* deklaration vieler variablen und kopie
    in globale arrays                       */
/*------------------------------------------*/
simulation_shell();

//
// Closes all files and finishes ...
//
 // SiSi::finalize();
return 0;
}

/*----------------------------------------------------*/
/*  write progress of single simulation in log file   */
/*----------------------------------------------------*/
extern "C" void progress_log(void)
{
SiSi::logFile.progress(Time-TimeStart,TimeEnd-TimeStart);
}

/*----------------------------------------------------*/
/*   write progress of total variation in log file   */
/*----------------------------------------------------*/
extern "C" void progress_var_log(long vstep)
{
SiSi::logFile.progress(vstep,num_total_variat);
}

/*-------------------------------------------------*/
/*  2 routines which split a string into 2 halfes  */
/*-------------------------------------------------*/
double RightFromSep (String s) {
    bool sepin = false;
    double val;
    String result = "";
    for( int i=0; i<s.length(); i++ )
      {
      if(sepin)
	result += s.charAt(i);
      if( s.charAt(i) ==':')
         sepin = true;
      }
    val = (double) atof(result);
    return (double) val;
}

double LeftFromSep (String s) {
    double val;
    String result = "";
    for( int i=0; i<s.length(); i++ )
      if( s.charAt(i) !=':')
	result += s.charAt(i);
      else
	break;

    val = (double) atof(result);
    return val;
}


/*----------------------------------------------------*/
/* writing list of variation variables in SiSi format
      using default values for min and max            */
/*----------------------------------------------------*/
void generate_varlist()
{
FILE  *sp;
int i;
double val;

cout<<"writing VariationVariables in varlist.ins...\n";
sp=fopen("varlist.ins","w");
for(i=0;i<num_variat_parser;i++)
  {
  fprintf(sp,"\tresult\t%s\n",VAR_NAMES[i]);
  fprintf(sp,"\t\ttype \tfloat\n");
  val=*VAR_VAL[i];
  fprintf(sp,"\t\t\\r %2.2f:%2.2f\n",0.5*val,1.5*val);
  fprintf(sp,"\t\tactive\ttrue\n");
  fprintf(sp,"\t\tprecision\t3\n");
  }
fclose(sp);
}

void generate_reslist()
{
FILE  *sp;
int i;
double val;

cout<<"writing StateVariables in reslist.ins...\n";
sp=fopen("reslist.ins","w");
for(i=0;i<num_states;i++)
  {
  fprintf(sp,"\tresult\t%s\n",state_names[i]);
  fprintf(sp,"\t\ttype \tfloat\n");
  fprintf(sp,"\t\tactive\ttrue\n");
  fprintf(sp,"\t\tprecision\t3\n");
  fprintf(sp,"\tresult\tS%s\n",state_names[i]);
  fprintf(sp,"\t\ttype \tfloat\n");
  fprintf(sp,"\t\tactive\tfalse\n");
  fprintf(sp,"\t\tprecision\t3\n");
  }
fclose(sp);
}
