#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static llint num_terms(int s,int sl)
{QD x={0.0,0.0,0.0,0.0},R={1.0,0.0,0.0,0.0}; int i,which; double STEP,TAR,u;
 which=whi[NUM_SUMS]; x[0]=DECAY[s][0];
 TAR=Pow(10,(double) (-wprec[NUM_SUMS]+sl));
 while (1)
 {if (x[0]>=EXPAND0_LIM[which]) get_wt_large(which,x,R,1,NUM_SUMS);
  if (R[0]>TAR) x[0]=x[0]+x[0]; else break;}
 STEP=x[0]*0.25; R[0]=1.0; x[0]-=STEP; STEP*=0.5;
 for (i=0;i<10;i++)
 {if (x[0]>=EXPAND0_LIM[which]) get_wt_large(which,x,R,1,NUM_SUMS);
  if (R[0]<TAR) x[0]-=STEP; else x[0]+=STEP; STEP*=0.5;}
 u=x[0]/DECAY[s][0]; if (u>1e15) return (((llint) 1)<<50);
 return (llint) (x[0]/DECAY[s][0]);}

static llint prepare_decay_hecke(int s,int W,int sl)
{QD o; if (DEBUG) printf("pdh %i %i %i\n",s,W,sl);
 if (!GLOBAL) return 0; QD_sqr(W,QD_twopi,o);
 QD_div(W,o,COND[s],o); QD_sqrt(W,o,DECAY[s]);
 if (VERBOSE>=2) {printf("DCY %i ",s); QD_output(W,wprec[NUM_SUMS],DECAY[s]);}
 return num_terms(s,sl);}

llint prepare_decay(int s,int W,int sl)
{QD o; if (DEBUG) printf("prepare_decay %i %i %i\n",s,W,sl);
 if (!GLOBAL) return 0; QD_powi(W,QD_twopi,s,o);
 if (s&1) QD_mul(W,o,QD_twopi,o); else QD_mul(W,o,QD_pi,o);
 QD_div(W,o,COND[s],o); QD_sqrt(W,o,DECAY[s]);
 if (VERBOSE>=2) {printf("DCY %i ",s); QD_output(W,wprec[NUM_SUMS],DECAY[s]);}
 return num_terms(s,sl);}

static llint prepare_sympow_hecke(int sp,int dv,int sl)
{llint NT=-1; ASSERT(CM_CASE!=0);
 if (DEBUG) printf("psh %i %i %i\n",sp,dv,sl);
 if (sp&1)
 {whi[NUM_SUMS]=wlo[NUM_SUMS]=WHICH;
 load_files_hecke(WHICH,sp,1+sp/2,dv); WHICH++;}
 else
 {wlo[NUM_SUMS]=WHICH; WHICH++; whi[NUM_SUMS]=WHICH;
  load_files_hecke(wlo[NUM_SUMS],sp,sp/2,dv);
  load_files_hecke(whi[NUM_SUMS],sp,1+sp/2,dv); WHICH++;}
 if (COND[sp][0]==0.0) compute_conductor_hecke(sp);
 NT=prepare_decay_hecke(sp,w[NUM_SUMS],sl);
 if (VERBOSE) printf("NT %id%i: %lli\n",sp,dv,NT);
 NUM_SUMS++; return NT;}

static llint prepare_sympow_check_hecke
(int sp,int dv,int prec,int nw,int sl,int bk)
{int i,W; llint NT; QD wig;
 if (DEBUG) printf("psch %i %i %i %i %i %i\n",sp,dv,prec,nw,sl,bk);
 wprec[NUM_SUMS]=prec; W=w[NUM_SUMS]=1+((prec-1)>>4); NUM_WIGS[NUM_SUMS]=nw;
 QD_copy(W,QD_one,WIGGLE[NUM_SUMS]); QD_copy(W,QD_one,WIGSQI[NUM_SUMS]);
 BLOCH_KATO[NUM_SUMS]=bk; NT=prepare_sympow_hecke(sp,dv,sl);
 if (!GLOBAL) return NT;
 for (i=0;i<nw;i++)
 {w[NUM_SUMS]=W; wprec[NUM_SUMS]=prec; QD_copy(W,QD_one,wig); wig[0]=0.5;
  QD_powi(1,wig,3+i,wig); QD_add(W,QD_one,wig,wig); SLOPPY[NUM_SUMS]=sl;
  if (VERBOSE>=2) {printf("wiggle %i ",i); QD_output(W,16*W,wig);}
  whi[NUM_SUMS]=whi[NUM_SUMS-1]; wlo[NUM_SUMS]=wlo[NUM_SUMS-1];
  BLOCH_KATO[NUM_SUMS]=BLOCH_KATO[NUM_SUMS-1];
  QD_copy(W,wig,WIGGLE[NUM_SUMS]); QD_div(W,QD_one,wig,WIGSQI[NUM_SUMS]);
  QD_sqr(W,WIGSQI[NUM_SUMS],WIGSQI[NUM_SUMS]); NUM_WIGS[NUM_SUMS]=nw;
  if (VERBOSE>=2) {printf("wigsqi %i ",i); QD_output(W,16*W,WIGSQI[NUM_SUMS]);}
  NUM_SUMS++;} return NT;}

static llint prepare_sympow(int sp,int dv,int sl)
{llint NT=-1; if (DEBUG) printf("prepare_sympow %i %i %i\n",sp,dv,sl);
 if (sp&1)
 {whi[NUM_SUMS]=wlo[NUM_SUMS]=WHICH;
  load_files(WHICH,sp,(sp+1)>>1,dv); WHICH++;}
 else
 {wlo[NUM_SUMS]=WHICH; WHICH++; whi[NUM_SUMS]=WHICH;
  load_files(wlo[NUM_SUMS],sp,sp>>1,dv);
  load_files(whi[NUM_SUMS],sp,1+(sp>>1),dv); WHICH++;}
 if (COND[sp][0]==0.0) compute_conductor(sp);
 NT=prepare_decay(sp,w[NUM_SUMS],sl);
 if (VERBOSE) printf("NT %id%i: %lli\n",sp,dv,NT);
 NUM_SUMS++; return NT;}

static llint prepare_sympow_check(int sp,int dv,int prec,int nw,int sl,int bk)
{int i,W; llint NT; QD wig;
 if (DEBUG) printf("psc %i %i %i %i %i %i\n",sp,dv,prec,nw,sl,bk);
 if (HECKE) return prepare_sympow_check_hecke(sp,dv,prec,nw,sl,bk);
 wprec[NUM_SUMS]=prec; W=w[NUM_SUMS]=1+((prec-1)>>4); NUM_WIGS[NUM_SUMS]=nw;
 QD_copy(W,QD_one,WIGGLE[NUM_SUMS]); QD_copy(W,QD_one,WIGSQI[NUM_SUMS]);
 BLOCH_KATO[NUM_SUMS]=bk;
 NT=prepare_sympow(sp,dv,sl); if (!GLOBAL) return NT;
 for (i=0;i<nw;i++)
 {w[NUM_SUMS]=W; wprec[NUM_SUMS]=prec; QD_copy(W,QD_one,wig); wig[0]=0.5;
  QD_powi(1,wig,3+i,wig); QD_add(W,QD_one,wig,wig); SLOPPY[NUM_SUMS]=sl;
  if (VERBOSE>=2) {printf("wiggle %i ",i); QD_output(W,16*W,wig);}
  whi[NUM_SUMS]=whi[NUM_SUMS-1]; wlo[NUM_SUMS]=wlo[NUM_SUMS-1];
  BLOCH_KATO[NUM_SUMS]=BLOCH_KATO[NUM_SUMS-1];
  QD_copy(W,wig,WIGGLE[NUM_SUMS]); QD_div(W,QD_one,wig,WIGSQI[NUM_SUMS]);
  QD_sqr(W,WIGSQI[NUM_SUMS],WIGSQI[NUM_SUMS]); NUM_WIGS[NUM_SUMS]=nw;
  if (VERBOSE>=2) {printf("wigsqi %i ",i); QD_output(W,16*W,WIGSQI[NUM_SUMS]);}
  NUM_SUMS++;} return NT;}

llint process_string(char *IN,llint UB)
{llint NT=0,nt,CP; char LAST[256],*NEXT,*CURR,*c,*n,COMMA[2]=",\0";
 int i,prec,sp,lim=0,nwigs,dvs=0,sl,bk,DO;
 c=CURR=malloc(256); n=NEXT=malloc(256); strcpy(CURR,IN);
 do
 {strcpy(LAST,CURR); ASSERT(ISA_NUMBER(CURR[0]));
  if (!ISA_NUMBER(CURR[1])) sp=(int) (CURR[0]-'0');
  else {sp=(int) ((CURR[0]-'0')*10+(CURR[1]-'0')); CURR++;}
  nwigs=1; CURR++; CP=NT; sl=0; bk=0; DO=FALSE;
  if (CURR[0]=='b') {bk=1; CURR++;}
  if (CURR[0]=='w') {nwigs=CURR[1]-'0'; CURR+=2;}
  if (CURR[0]=='s') {sl=CURR[1]-'0'; CURR+=2;}
  if (CURR[0]=='P') {DO=TRUE;} else ASSERT(CURR[0]=='p');
  if (CURR[1]=='z')
  {if (!RERUN) {prec=6; ZEROCHECK=1;} else {prec=12; ZEROCHECK=0;} CURR+=2;}
  else
  {if (ISA_NUMBER(CURR[2])) {prec=10*(CURR[1]-'0')+(CURR[2]-'0'); CURR+=3;}
   else {prec=(CURR[1]-'0'); CURR+=2;}}
  if ((sp&1)==0)
  {if ((CURR[0]=='d') || (CURR[0]=='D')) errorit("Derivative with even power");
   dvs=1; nt=prepare_sympow_check(sp,0,prec,nwigs,sl,bk); if (nt>NT) NT=nt;}
  else
  {if (CURR[0]=='D')
   {lim=CURR[1]-'0'; CURR+=2; dvs=1;
    nt=prepare_sympow_check(sp,lim,prec,nwigs,sl,bk); if (nt>NT) NT=nt;}
   else
   {ASSERT(CURR[0]=='d'); lim=CURR[1]-'0'; CURR+=2; dvs=(lim+1);
    for (i=0;i<=lim;i++)
    {nt=prepare_sympow_check(sp,i,prec,nwigs,sl,bk); if (nt>NT) NT=nt;}}}
  if ((NT>UB) && (DO)) UB=NT;
  if ((NT>UB) || (((sp&1)==0) && NO_QT))
  {for (i=0;((LAST[i]!=0) && (LAST[i]!=','));i++); LAST[i]=0;
   if (VERBOSE)
   {if (NT>UB) printf("Ignoring %s due to bound\n",LAST);
    else printf("Ignoring %s since curve is non qtwist-minimal\n",LAST);}
   NT=CP; NUM_SUMS-=dvs*(1+nwigs);}
  else CP=NT;
  NEXT=strstr(CURR,COMMA); CURR=NEXT+1; } while(NEXT);
 free(c); free(n); return(NT);}

llint preparation(int K,char *IN,llint UB)
{int i; llint NT;
 POWSER=malloc(K*sizeof(QD**)); TABLE=malloc(K*sizeof(QD**));
 SYMPOW=malloc(K*sizeof(int)); evalpt=malloc(K*sizeof(int));
 BLOCH_KATO=malloc(K*sizeof(int)); MESH_COUNT=malloc(K*sizeof(int));
 whi=malloc(K*sizeof(int)); wlo=malloc(K*sizeof(int));
 derivative=malloc(K*sizeof(int)); HALF_ZERO=malloc(K*sizeof(QD));
 WIGGLE=malloc(K*sizeof(QD)); WIGSQI=malloc(K*sizeof(QD));
 EXPAND0_LIM=malloc(K*sizeof(QD)); STEP_SIZE=malloc(K*sizeof(QD));
 TOO_BIG=malloc(K*sizeof(QD)); NUM_LOGS=malloc(K*sizeof(QD));
 TACKON=malloc(K*sizeof(int)); TACKS=malloc(K*sizeof(QD*));
 w=malloc(K*sizeof(int)); wprec=malloc(K*sizeof(int));
 DECAY=malloc(K*sizeof(QD)); NUM_WIGS=malloc(K*sizeof(QD));
 for (i=0;i<K;i++) QD_copy(wmax,QD_zero,COND[i]); WHICH=0; NUM_SUMS=0;
 NT=process_string(IN,UB);
 if (MODDEG) {if (MD_SPEED>0.0) NT=(llint) ((double) NT/MD_SPEED);
              else NT=postpare_moddeg();}
 printf("Maximal number of terms is %lli\n",NT);
 if (NT<MIN_TERMS) NT=MIN_TERMS; return(NT);}
