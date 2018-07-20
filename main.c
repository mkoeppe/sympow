#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

int main(int argc,char **argv)
{char INSTRING[1024]="2w3s1p32,3bp16d1,4p8\0",CSTR[1024]="\0",LSTR[1024]="\0";
 char TYPE[16]; int i=1; llint NT,UB=(((llint) 1)<<45);
 int NO_CM=FALSE,info=0,ROOTNO=FALSE,SLOPPY=0,QD_CHECK;
 NO_QT=FALSE; VERBOSE=TRUE; GLOBAL=TRUE; HECKE=FALSE; TWIST=FALSE; AP_SAVE=0;
 CM_CASE=FALSE; GET=malloc(1024); COND0=1; fp3=0; fp2=0; MAX_TABLE=1<<27;
 MODDEG=FALSE; ANAL_RANK=FALSE; ZEROCHECK=FALSE; RERUN=FALSE; QD_CHECK=TRUE;
#if defined(ISOC99_FENV) || defined(FPUCONTROLH) || defined(x86)
 fpu_53bits();
#endif
 strcpy(TYPE,"RELEASE"); MD_SPEED=2.0;
 while(i<argc)
 {if (!strcmp(argv[i],"-quiet")) {VERBOSE=FALSE; i++;}
  else if ((!strcmp(argv[i],"-sympow")) || (!strcmp(argv[i],"-sp")))
  {strcpy(INSTRING,argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-verbose")) {VERBOSE=2; i++;}
  else if (!strcmp(argv[i],"-help")) help_message();
  else if (!strcmp(argv[i],"-curve")) {strcpy(CSTR,argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-label")) {strcpy(LSTR,argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-info")) {info=i+1; i+=3;}
  else if (!strcmp(argv[i],"-bound")) {UB=atoll(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-maxtable")) {MAX_TABLE=atoll(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-hecke")) {HECKE=TRUE; i++;}
  else if (!strcmp(argv[i],"-mdspeed")) {MD_SPEED=atof(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-nocm")) {NO_CM=TRUE; i++;}
  else if (!strcmp(argv[i],"-noqdcheck")) {QD_CHECK=FALSE; i++;}
  else if (!strcmp(argv[i],"-noqt")) {NO_QT=TRUE; i++;}
  else if (!strcmp(argv[i],"-local")) {GLOBAL=FALSE; i++;}
  else if (!strcmp(argv[i],"-moddeg")) {MODDEG=TRUE; i++;}
  else if (!strcmp(argv[i],"-analrank")) {ANAL_RANK=TRUE; i++;}
  else if (!strcmp(argv[i],"-sloppy")) {SLOPPY=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-rootno")) {ROOTNO=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-rewarp")) {rewarp_params(); return 0;}
#ifdef NEW_DATA
  else if (!strcmp(argv[i],"-shell1")) {new_sympow_s1(argv[i+1]); return 0;}
  else if (!strcmp(argv[i],"-pari")) {new_sympow_pari(argv[i+1]); return 0;}
  else if (!strcmp(argv[i],"-shell2")) {new_sympow_s2(argv[i+1]); return 0;}
  else if (!strcmp(argv[i],"-txt2bin"))
  {txt2bin(atoi(argv[i+1]),argv[i+2],stdin); return 0;}
  else if (!strcmp(argv[i],"-new_data")) {new_data(argv[i+1]); return 0;}
#else
  else if ((!strcmp(argv[i],"-shell1")) || (!strcmp(argv[i],"-pari")) || 
	   (!strcmp(argv[i],"-shell2")) || (!strcmp(argv[i],"-txt2bin")) ||
	   (!strcmp(argv[i],"-new_data")))
    errorit("new_data not possible --- try re-configuring SYMPOW");
#endif
  else errorit("Command not recognised");}

 if (VERBOSE) printf("sympow %s %s  (c) Mark Watkins -",VERSION,TYPE);
 if (QD_CHECK) QD_check();
 if (VERBOSE) {printf("-- see README and COPYING for details\n");}
 init_primes(); if (!CSTR[0]) getline0(stdin,CSTR,256); curve_init(CSTR,LSTR);
 if (HECKE && !CM_CASE) errorit("Curve does not have CM");
 if (NO_CM && CM_CASE) errorit("Curve has CM");
 if (VERBOSE && CM_CASE)
   printf("Curve has complex multiplication %i\n",CM_CASE);
 if (info) {localinfos(argv[info],argv[info+1]); return 0;}
 if (ROOTNO) {printf("Root number for sp:%i is %i\n",
		     ROOTNO,global_rootno(ROOTNO)); return 0;}
 if (MODDEG) {HECKE=FALSE; prepare_moddeg(INSTRING);}
 if (ANAL_RANK) {HECKE=FALSE; prep_analrank(UB,SLOPPY); return 0;}
 if (VERBOSE && HECKE) printf("Using Hecke powers\n");
 NT=preparation(CURR_MAX,INSTRING,UB);
 if (VERBOSE>=2) printf("NUM_SUMS is %i\n",NUM_SUMS);
 if (NUM_SUMS>CURR_MAX) errorit("Too many sums!");
 if (NUM_SUMS==0) printf("Nothing to do!\n");
 fflush(stdout); if (GLOBAL && NUM_SUMS) go(NT,NT);
 if (ZEROCHECK)
 {NUM_SUMS=0; WHICH=0; free_data(); RERUN=TRUE;
  NT=process_string(INSTRING,UB); go(NT,NT);}
 return 0;}
