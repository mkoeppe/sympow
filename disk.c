#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

int CHEAT=FALSE;
static void read_file_mesh_bin(FILE *A,int which)
{int i,j; QD **P,*Q; if (!GLOBAL) return;
 P=TABLE[which]=malloc(MESH_COUNT[which]*sizeof(QD*));
 for(i=0;i<MESH_COUNT[which];i++)
 {Q=P[i]=malloc((1+LPT)*sizeof(QD));
  for(j=0;j<=LPT;j++) fread(Q[j],sizeof(QD),1,A);} }

static void read_file_series(FILE *A,int which)
{int i,j,k,AB; char INPUT[64],BLAH[64]; QD **P,*Q; if (!GLOBAL) return;
 P=POWSER[which]=malloc((1+NUM_LOGS[which])*sizeof(QD*));
 for (i=0;i<=NUM_LOGS[which];i++)
 {AB=A0PT; if ((i==NUM_LOGS[which]) && (HALF_ZERO[which])) AB>>=1;
  Q=P[i]=malloc(AB*sizeof(QD));
  for(j=0;j<AB;j++)
  {fscanf(A,"%s %s",INPUT,BLAH);
   for(k=0;k<4;k++) {fscanf(A,"%s",INPUT); Q[j][k]=atof(INPUT);}}} 
 for (i=0;i<TACKON[which];i++)
 {fscanf(A,"%s %s",INPUT,BLAH); ASSERT(!strcmp(INPUT,"TACKS"));
  for(k=0;k<4;k++) {fscanf(A,"%s",INPUT); TACKS[which][i][k]=atof(INPUT);}}}

static int get_params(int which,int sp,int ep,int dv)
{FILE *F; int ms=0; char TOK[16]=",",LINE[64],S[16],*STR,U[16];
 if (!GLOBAL) return 0;
 if (HECKE)
 {HALF_ZERO[which]=FALSE; NUM_LOGS[which]=0;
  if (ep==1) S[0]='P'; else S[0]='H';}
 else
 {if (sp&1) HALF_ZERO[which]=FALSE;
  else if (ep&1) HALF_ZERO[which]=2; else HALF_ZERO[which]=1;
  NUM_LOGS[which]=sp/2; S[0]='P';}
 if (!HECKE) {S[1]='0'+(sp/10); S[2]='0'+(sp%10);}
 else {S[1]='0'+(ep/10); S[2]='0'+(ep%10);}
 if ((sp&1) || HECKE)
 {if (dv==0) S[3]='E'; if (dv==1) S[3]='O'; if (dv>1) S[3]='0'+dv;}
 else if (((sp&3)==0) && CM_CASE) {if (2*ep==sp) S[3]='l'; else S[3]='h';}
 else {if (2*ep==sp) S[3]='L'; else S[3]='H';}
 if (HECKE && dv) {TACKS[which]=malloc(dv*sizeof(QD)); TACKON[which]=dv;}
 else if (dv<=sp/2) TACKON[which]=0;
 else {TACKS[which]=malloc((dv-sp/2)*sizeof(QD)); TACKON[which]=dv-sp/2;}
 S[4]=0; F=fopen("datafiles/param_data","r"); strcpy(U,S);
 if (ANAL_RANK) {if (dv>0) U[0]='A'; else U[0]='m';}
 if (MODDEG) {if (HECKE) U[0]='m'; else U[0]='M';}
 while (1)
 {if (!getline0(F,LINE,64))
  {printf("**ERROR** %s not found in param_data file\n",S);
   if (!HECKE) printf("It can be added with './sympow -new_data %i",sp);
   else printf("It can be added with './sympow -new_data %i",ep);
   if ((sp&1) || HECKE) printf("d%i",dv); if ((HECKE) && (sp>1)) printf("h");
   if ((CM_CASE) && ((sp&3)==0)) printf("c"); printf("'\n"); exit(-1);}
  if ((S[0]==LINE[0]) && (S[1]==LINE[1]) && (S[2]==LINE[2]) && (S[3]==LINE[3])
      && !assure_line(LINE))
  {printf("Problem with datafiles/param_data for %s\n",S);
   errorit("datafiles/param_data entry corrupted!");}
  STR=strtok(LINE,TOK);
  if (!strcmp(STR,S))
  {if (VERBOSE>=2) printf("S Found: %s\n",S);
   TOO_BIG[which]=atof(strtok(NULL,TOK));
   STEP_SIZE[which]=QD_2pow(atoi(strtok(NULL,TOK)));
   ms=atoi(strtok(NULL,TOK)); CHEAT=FALSE; break;}
  if (!strcmp(STR,U))
  {if (VERBOSE>=2) printf("U Found: %s\n",U);
   TOO_BIG[which]=atof(strtok(NULL,TOK));
   STEP_SIZE[which]=QD_2pow(atoi(strtok(NULL,TOK)));
   ms=atoi(strtok(NULL,TOK)); CHEAT=TRUE; break;}}
 EXPAND0_LIM[which]=31.5*STEP_SIZE[which]; fclose(F); return 32*ms;}

void load_files(int which,int sp,int ep,int dv)
{FILE *A,*N; char NM[32],NAME[32]="datafiles/P"; int dl=9;
 if (DEBUG) printf("load_files %i %i %i %i\n",which,sp,ep,dv);
 SYMPOW[which]=sp; evalpt[which]=ep; derivative[which]=dv;
 if (!GLOBAL) return; MESH_COUNT[which]=get_params(which,sp,ep,dv);
 if (sp<10) {NAME[2+dl]='0'; NAME[3+dl]='0'+sp;}
 else {NAME[2+dl]='0'+sp/10; NAME[3+dl]='0'+(sp%10);}
 NAME[5+dl]='M';
 if (CHEAT) {if (ANAL_RANK) {if (dv>0) NAME[1+dl]='A'; else NAME[1+dl]='m';}
             else if (HECKE) NAME[1+dl]='m'; else NAME[1+dl]='M';}
 if (sp&1) {if (dv==0) NAME[4+dl]='E'; if (dv==1) NAME[4+dl]='O';
 if (dv>1) NAME[4+dl]='0'+dv;}
 else {if (2*ep==sp) NAME[4+dl]='L'; else NAME[4+dl]='H';}
 NAME[6+dl]='.'; NAME[7+dl]='b'; NAME[8+dl]='i'; NAME[9+dl]='n'; NAME[10+dl]=0;
 if (((sp&3)==0) && CM_CASE) NAME[4+dl]+=('a'-'A');
 if (VERBOSE>=2) printf("%i Reading %s\n",which,NAME);
 A=fopen(NAME,"r");
 if (!A)
 {strcpy(NM,NAME); NM[7+dl]='t'; NM[8+dl]='x'; NM[9+dl]='t'; N=fopen(NM,"r");
  if (VERBOSE) printf("Creating %s from %s\n",NAME,NM);
  txt2bin(MESH_COUNT[which],NAME,N); A=fopen(NAME,"r"); fclose(N);}
 read_file_mesh_bin(A,which); fclose(A);
 NAME[5+dl]='S'; NAME[7+dl]='t'; NAME[8+dl]='x'; NAME[9+dl]='t';
 if (VERBOSE>=2) printf("%i Reading %s\n",which,NAME);
 A=fopen(NAME,"r"); read_file_series(A,which); fclose(A);}

void load_files_hecke(int which,int sp,int ep,int dv)
{FILE *A,*N; char NM[32],NAME[32]="datafiles/H"; int dl=9;
 if (DEBUG) printf("load_files_hecke %i %i %i %i\n",which,sp,ep,dv);
 if (ep==1) {load_files(which,1,1,dv); return;}
 SYMPOW[which]=sp; evalpt[which]=ep; derivative[which]=dv;
 if (!GLOBAL) return; MESH_COUNT[which]=get_params(which,sp,ep,dv);
 if (ep<10) {NAME[2+dl]='0'; NAME[3+dl]='0'+ep;}
 else {NAME[2+dl]='0'+ep/10; NAME[3+dl]='0'+(ep%10);}
 NAME[5+dl]='M'; if (CHEAT) NAME[1+dl]='m';
 if (dv==0) NAME[4+dl]='E'; if (dv==1) NAME[4+dl]='O';
 if (dv>1) NAME[4+dl]='0'+dv;
 NAME[6+dl]='.'; NAME[7+dl]='b'; NAME[8+dl]='i';
 NAME[9+dl]='n'; NAME[10+dl]=0;
 if (VERBOSE>=2) printf("%i Reading %s\n",which,NAME);
 A=fopen(NAME,"r");
 if (!A)
 {strcpy(NM,NAME); NM[7+dl]='t'; NM[8+dl]='x'; NM[9+dl]='t'; N=fopen(NM,"r");
  if (VERBOSE) printf("Creating %s from %s\n",NAME,NM);
  txt2bin(MESH_COUNT[which],NAME,N); A=fopen(NAME,"r"); fclose(N);}
 read_file_mesh_bin(A,which); fclose(A);
 NAME[5+dl]='S'; NAME[7+dl]='t'; NAME[8+dl]='x'; NAME[9+dl]='t';
 if (VERBOSE>=2) printf("%i Reading %s\n",which,NAME);
 A=fopen(NAME,"r"); read_file_series(A,which); fclose(A);}

int getline0(FILE *F,char *v,int l)
{GET=fgets(v,l,F); if (NULL==GET) return(0); return(1);}

