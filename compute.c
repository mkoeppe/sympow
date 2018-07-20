#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

llint BOUND; int SQBND,*small_prime;
extern double EXPANSION_AROUND_ZERO_LIMIT;
QD *SUM,**coeff_small,***coeff_primepower;
#define NS CURR_MAX

static void malloc_arrays(int SQBND,int NUMP)
{int i,j; QD **P; if (DEBUG) printf("malloc_arrays %i %i\n",SQBND,NUMP);
 SUM=malloc(NUM_SUMS*sizeof(QD)); coeff_small=malloc(SQBND*sizeof(QD*));
 for (i=0;i<SQBND;i++) coeff_small[i]=malloc(NUM_SUMS*sizeof(QD));
 coeff_primepower=malloc(NUMP*sizeof(QD**));
 for (i=0;i<NUMP;i++)
 {P=coeff_primepower[i]=malloc(40*sizeof(QD*));
  for (j=0;j<40;j++) P[j]=malloc(NUM_SUMS*sizeof(QD));}}

static void free_arrays(int SQBND,int NUMP)
{int i,j; QD **P; if (DEBUG) printf("free_arrays %i %i\n",SQBND,NUMP);
 free(SUM); for (i=0;i<SQBND;i++) free(coeff_small[i]); free(coeff_small);
 for (i=0;i<NUMP;i++)
 {P=coeff_primepower[i]; for (j=0;j<40;j++) free(P[j]); free(P);}
 free(coeff_primepower);}

static void iterate(llint n,int nump,int exponent,QD *bp,QD *bq)
{int j,S; llint pj=0,nn; QD WT,bn[NS];
 for(S=0;S<NUM_SUMS;S++)
 {get_weight(n,WT,S);
  if ((S>0) && (SYMPOW[whi[S]]==SYMPOW[whi[S-1]]))
    QD_copy(w[S],bn[S-1],bn[S]); else QD_mul(w[S],bp[S],bq[S],bn[S]);
  QD_mul(w[S],bn[S],WT,WT); QD_add(w[S],SUM[S],WT,SUM[S]);
  if (n<SQBND) QD_copy(w[S],bn[S],coeff_small[n][S]);}
 for (j=0;j<=nump;j++)
 {pj=PRIMES[j];
  if (BOUND/pj<n) j=nump+1;
  else
  {nn=n*(llint) pj;
   if (j!=nump) iterate(nn,j,1,coeff_primepower[j][1],bn);
   else iterate(nn,nump,exponent+1,coeff_primepower[nump][exponent+1],bq);}}}

static void make_primep_coeffs(int s,llint p,llint B)
{QDpoly v,v2; int j=0,S,n,sp;
 n=(int) Floor(Log((double) B)/Log((double) p));
 for(S=0;S<NUM_SUMS;S++)
 {sp=SYMPOW[whi[S]]; /* don't recompute ap?!  --- slow but negligible?! */
  if ((S==0) || (sp!=SYMPOW[whi[S-1]]))
  {euler_factor(p,sp,&v); QDpoly_inv(v,n,&v2);
   for (j=1;j<=n;j++) QD_copy(w[S],v2.coeff[j],coeff_primepower[s][j][S]);
   delQDpoly(&v); delQDpoly(&v2);}}}

llint ec_do(llint p)
{if ((C4C6LL) && (p<(1<<29)) && (p>500))
   return ec_bsgs_ap_AB((-27*Ec4ll)%p,(-54*Ec6ll)%p,p);
 return ec_ap(Ec4,Ec6,p);}

static void get_large_coeffs(llint p,int pnum)
{int S,ap;
 if ((COND0%p)==0) {make_primep_coeffs(0,p,p+64); return;}
 if (pnum<AP_SAVE) ap=(int) apsave[pnum]; else ap=(int) ec_do(p);
 for(S=0;S<NUM_SUMS;S++)
 {if ((S>0) && (SYMPOW[whi[S]]==SYMPOW[whi[S-1]]))
    QD_copy(w[S],coeff_primepower[0][1][S-1],coeff_primepower[0][1][S]);
  else pth_coeff_from_ap(w[S],p,ap,SYMPOW[whi[S]],coeff_primepower[0][1][S]);}
 if (AP_SAVE) {apsave[pnum]=ap; if (pnum>AP_SAVE) AP_SAVE=pnum;}}

static void addup(llint n,QD *bp,QD *bq)
{int S; QD WT,t;
 for(S=0;S<NUM_SUMS;S++)
 {get_weight(n,WT,S);
  if ((S==0) || (SYMPOW[whi[S]]!=SYMPOW[whi[S-1]])) QD_mul(w[S],bp[S],bq[S],t);
  QD_mul(w[S],t,WT,WT); QD_add(w[S],SUM[S],WT,SUM[S]);}}

static int digits(int prec)
{if (prec<=4) return 3; if (prec<=6) return 4;
 if (prec<=8) return 5; if (prec<=10) return 6;
 if (prec<=16) return prec-4; if (prec<=24) return prec-5;
 if (prec<=32) return prec-6; if (prec<=48) return prec-7;
 if (prec<=64) return prec-8; return 0;}

static void sum_print(int S,int t)
{int wh=whi[S],sp=SYMPOW[wh],dv=derivative[wh]; char wig;
 if (WIGGLE[S][0]==1.0) wig='n'; else wig='w'; if (t==1) wig='d';
 printf("%2i%c%i",sp,wig,dv); if (t==0) printf(": "); else printf(" ");}

static void results()
{int i,S,z; double Y; QD X;int ZERO[NS];
 for (i=0;i<NS;i++) ZERO[i]=0; z=0; printf("Computed: ");
 for (S=0;S<NUM_SUMS;S+=NUM_WIGS[S]+1) sum_print(S,1);
 printf("\n"); printf("Checked out: ");
 for (S=0;S<NUM_SUMS;S+=NUM_WIGS[S]+1)
 {Y=0.0; for (i=0;i<NUM_WIGS[S];i++)
  {QD_sub(w[S],SUM[S],SUM[S+i+1],X); if (Abs(X[0])>Y) Y=Abs(X[0]);}
  if ((NUM_WIGS[S]) && (Log10(Y)<(double) -digits(wprec[S]))) sum_print(S,1);
  if (Log10(Abs(SUM[S][0]))<(double) -digits(wprec[S])) ZERO[z++]=S+1;}
 printf("\n");
 if ((ZERO[0]!=0) && !ZEROCHECK)
 {printf("Near Zero: "); for(i=0;i<z;i++) sum_print(ZERO[i],1); printf("\n");}
 if (ZERO[0]==0) ZEROCHECK=0;
 for (S=0;S<NUM_SUMS;S++)
 {if (BLOCH_KATO[S]) special_value(w[S],SUM[S],SYMPOW[whi[S]]);
  sum_print(S,0); QD_output(w[S],wprec[S],SUM[S]);}}

double go(llint lim,llint BAK)
{int S,s,mult,pnum; QD CHECK,LAST,NOW,DIFF;
 if (DEBUG) printf("%lli %lli\n",lim,BAK);
 llint NEXT,pi,multpi,*auxp,PMAX,PSIZE=1000000,NUM_PRIMES;
 QD *QD_one_array; double f; int ptotal,last=0;
 if (!HECKE) IFACT_INIT(lim+1+Sqrt(4*lim));
 QD_one_array=malloc(NUM_SUMS*sizeof(QD));
 auxp=malloc(80000*sizeof(llint));
 for (S=0;S<NUM_SUMS;S++) QD_copy(w[S],QD_one,QD_one_array[S]);
 BOUND=lim; if (BOUND<MIN_TERMS) BOUND=MIN_TERMS;
 SQBND=(int) Ceil(Sqrt(BAK)); if (SQBND<MIN_SQRT) SQBND=MIN_SQRT;
 if ((SQBND<30000) && (TWIST)) SQBND=30000;
 NUM_PRIMES=(llint) Floor((double) SQBND/Log((double) SQBND)*
			(1.0+1.5/Log((double) SQBND)));

 malloc_arrays(SQBND,NUM_PRIMES); NEXT=BOUND/4; QD_copy(wmax,QD_zero,LAST);
 if (VERBOSE>=2) printf("BOUNDS %lli %i\n",BOUND,SQBND);
 if (VERBOSE) QD_output(w[0],wprec[0],SCALE);
 for(S=0;S<NUM_SUMS;S++)
 {get_weight(1,SUM[S],S); QD_copy(w[S],QD_one,coeff_small[1][S]);}
 for (s=0;s<NUM_PRIMES;s++) make_primep_coeffs(s,PRIMES[s],BOUND);
 if (VERBOSE>=2) printf("made pp\n");
 for (s=0;s<NUM_PRIMES;s++)
   iterate(PRIMES[s],s,1,coeff_primepower[s][1],QD_one_array);
 if (VERBOSE) printf("Done with small primes %i\n",PRIMES[s-1]);
 pi=PRIMES[s];
 if ((int) Sqrt((double) (2*PSIZE))>SQBND)
 {PSIZE=(llint) SQBND; PSIZE=PSIZE*PSIZE/2;}
 get_primes_ll(pi,PSIZE,auxp); PMAX=pi+PSIZE; pnum=0; ptotal=0;
 while (pi<BOUND)
 {pi=auxp[pnum++]; ptotal++;
  if (pi>PMAX)
  {get_primes_ll(pi+2,PSIZE,auxp); pnum=0; PMAX+=PSIZE;
   if (VERBOSE) printf("At prime %lli %.16f\n",pi,SUM[0][0]); fflush(stdout);}
  get_large_coeffs(pi,ptotal); mult=1; multpi=pi;
  while (multpi<BOUND)
  {addup(multpi,coeff_primepower[0][1],coeff_small[mult]); mult++; multpi+=pi;}
  if ((MODDEG) && (pi>NEXT))
  {NEXT+=BOUND/32;
   QD_copy(wmax,QD_zero,CHECK); QD_mul(w[0],SUM[0],SCALE,CHECK);
   if (MANIN_TWIST!=1) QD_mul1(w[0],CHECK,(double) MANIN_TWIST,CHECK);
   QD_mul1(w[0],CHECK,(double) EVEN_TOP,CHECK);
   QD_div1(w[0],CHECK,(double) EVEN_BOTTOM,CHECK);
   if (VERBOSE) {printf("At %lli: ",pi); QD_output(w[0],wprec[0],CHECK);}
   QD_round(wmax,CHECK,NOW); QD_sub(wmax,CHECK,NOW,DIFF);
   if (Abs(DIFF[0])<0.01)
   {QD_sub(wmax,NOW,LAST,DIFF);
    if (DIFF[0]!=0.0) {QD_copy(wmax,NOW,LAST); last=1;}
    else
    {last++; if ((last==3) && (MD_SPEED!=0.0))
     {QD_mul(wmax,TW_EFF,NOW,CHECK);
      if (MANIN_TWIST!=1) QD_div1(wmax,CHECK,(double) MANIN_TWIST,CHECK);
      if (MANIN!=1) QD_mul1(wmax,CHECK,(double) MANIN,CHECK);
      QD_round(wmax,CHECK,CHECK);
      printf("Modular Degree is "); QD_intout_noplus(CHECK); printf("\n");
      goto DONE;}}}
   else {QD_copy(wmax,QD_zero,LAST); last=0;}}}
 if (MODDEG)
 {QD_copy(wmax,QD_zero,CHECK); QD_mul(w[0],SUM[0],SCALE,CHECK);
  if (MANIN_TWIST!=1) QD_mul1(w[0],CHECK,(double) MANIN_TWIST,CHECK);
  QD_mul1(w[0],CHECK,(double) EVEN_TOP,CHECK);
  QD_div1(w[0],CHECK,(double) EVEN_BOTTOM,CHECK);
  if (VERBOSE) {printf("At end "); QD_output(w[0],wprec[0],CHECK);}
  QD_round(wmax,CHECK,NOW); QD_sub(wmax,CHECK,NOW,DIFF);
  if (Abs(DIFF[0])<0.05)
  {QD_mul(wmax,TW_EFF,NOW,CHECK);
   if (MANIN_TWIST!=1) QD_div1(wmax,CHECK,(double) MANIN_TWIST,CHECK);
   if (MANIN!=1) QD_mul1(wmax,CHECK,(double) MANIN,CHECK);
   QD_round(wmax,CHECK,CHECK);
   printf("Modular Degree is "); QD_intout_noplus(CHECK); printf("\n");
   goto DONE;} errorit("Not converging to an integer in modular degree\n");}
 else if (!ANAL_RANK) results();
 if (ANAL_RANK==-1) {results(); errorit("Not programmed to test if 0\n");}
 DONE: f=SUM[0][0]; free(QD_one_array); free(auxp);
 free_arrays(SQBND,NUM_PRIMES); return f;}
