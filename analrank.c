#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static int primes_upper_bound(llint x)
{double X=(double) x; return (int) (X/Log(X)*(1.0+1.2762/Log(X)));}

void prep_analrank(llint UB,int sl)
{char STR[128]="1w0s0p06D0",S[128]; llint NU,NT; double f; int l,w1;
 l=(int) Floor(Log10(REAL_PERIOD[0])); if (l>0) l=0;
 if (l<-10-sl) {l=-10-sl; printf("*WARNING* possible precision truncation\n");}
 STR[7]-=l; while (STR[7]>'9') {STR[7]-=10; STR[6]++;}
 STR[4]+=sl; w1=global_rootno(1);
 if (w1==-1) STR[9]++; NT=preparation(CURR_MAX,STR,UB);
 if (NT==0) errorit("Too many terms needed to compute the analytic rank");
 if (NT>(1<<30))
 {ANAL_RANK=-1; NUM_SUMS=0; WHICH=0;
  sprintf(S,"%s,%s,%s,%s,%s",STR,STR,STR,STR,STR);
  S[20]+=2; S[31]+=4; S[42]+=6; S[53]+=8;
  NT=process_string(S,UB); go(NT,NT); return;}
 else {AP_SAVE=1; apsave=malloc(primes_upper_bound(NT)*sizeof(int));}
 if (NT<MIN_TERMS) NT=MIN_TERMS; f=go(NT,NT); free_data();
 if (Abs(f)/REAL_PERIOD[0]>0.0001)
 {if (w1==1) printf("Analytic Rank is 0 : L-value %.5e\n",f);
  else printf("Analytic Rank is 1 : L'-value %.5e\n",f); return;}
 if (VERBOSE)
 {if (w1==1) printf("L(E,1) appears to be zero %.5e\n",f);
  else printf("L'(E,1) appears to be zero %.5e\n",f);}
 NUM_SUMS=0; WHICH=0; STR[9]+=2; NU=process_string(STR,UB);
 f=go(NU,NT); free_data();
 if (Abs(f)/REAL_PERIOD[0]>0.0001)
 {if (w1==1) printf("Analytic Rank is 2 : leading L-term %.5e\n",f);
  else printf("Analytic Rank is 3 : leading L-term %.5e\n",f); return;}
 if (VERBOSE)
 {if (w1==1) printf("L''(E,1) appears to be zero %.5e\n",f);
  else printf("L'''(E,1) appears to be zero %.5e\n",f);}
 NUM_SUMS=0; WHICH=0; STR[9]+=2; NU=process_string(STR,UB);
 f=go(NU,NT); free_data();
 if (Abs(f)/REAL_PERIOD[0]>0.0001)
 {if (w1==1) printf("Analytic Rank is 4 : leading L-term %.5e\n",f);
  else printf("Analytic Rank is 5 : leading L-term %.5e\n",f); return;}
 if (VERBOSE)
 {if (w1==1) printf("4th deriv of L(E,s) at s=1 appears to be zero %.5e\n",f);
  else printf("5th deriv of L(E,s) at s=1 appears to be zero %.5e\n",f);}
 NUM_SUMS=0; WHICH=0; STR[9]+=2; NU=process_string(STR,UB);
 f=go(NU,NT); free_data();
 if (Abs(f)/REAL_PERIOD[0]>0.0001)
 {if (w1==1) printf("Analytic Rank is 6 : leading L-term %.5e\n",f);
  else printf("Analytic Rank is 7 : leading L-term %.5e\n",f); return;}
 if (VERBOSE)
 {if (w1==1) printf("6th deriv of L(E,s) at s=1 appears to be zero %.5e\n",f);
  else printf("7th deriv of L(E,s) at s=1 appears to be zero %.5e\n",f);}
 NUM_SUMS=0; WHICH=0; STR[9]+=2; NU=process_string(STR,UB);
 f=go(NU,NT); free_data();
 if (Abs(f)/REAL_PERIOD[0]>0.0001)
 {if (w1==1) printf("Analytic Rank is 8 : leading L-term %.5e\n",f);
  else printf("Analytic Rank is 9 : leading L-term %.5e\n",f); return;}
 if (VERBOSE)
 {if (w1==1) printf("8th deriv of L(E,s) at s=1 appears to be zero %.5e\n",f);
  else printf("9th deriv of L(E,s) at s=1 appears to be zero %.5e\n",f);}
 errorit("Analytic Rank is too big!!");}
