#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)
static int my_primes[179]=
  {2,3,5,7,11,13,17,19,23,29,31,37, 41,43,47,53,59,61,67,71,
   73,79,83,89,97,101,103,107,109,113, 127,131,137,139,149,151,157,
   163,167,173,179,181,191,193,197,199, 211,223,227,229,233,239,241,
   251,257,263,269,271,277,281,283,293, 307,311,313,317,331,337,347,
   349,353,359,367,373,379,383,389,397, 401,409,419,421,431,433,439,
   443,449,457,461,463,467,479,487,491, 499,503,509,521,523,541,547,
   557,563,569,571,577,587,593,599,601, 607,613,617,619,631,641,643,
   647,653,659,661,673,677,683,691,701, 709,719,727,733,739,743,751,
   757,761,769,773,787,797,809,811,821, 823,827,829,839,853,857,859,
   863,877,881,883,887,907,911,919,929, 937,941,947,953,967,971,977,
   983,991,997,1009,1013,1019,1021,1031, 1033,1039,1049,1051,1061,1063};

static char *aux;

static int get_primes_loop(int start,int size,int *A)
{int i,l,j=0,k=0; int s,p=2;
 if (DEBUG>=2) printf("init_primesj %i %i\n",start,size);
 for (i=0;i<size;i++) aux[i]=1;
 while(p<=Sqrt(start+size+100))
 {s=start-p*(start/p); if (s!=0) l=p-s; else l=0;
  for (i=l;i<size;i+=p) aux[i]=0; p=PRIMES[j++];}
 for (i=0;i<size;i++) {if (aux[i]) {A[k++]=start+i;}} return k;}

static void get_primes(int start,int sz,int *A)
{int i,size=(int) sz,k=0; int *B;
 if (DEBUG>=2) printf("init_primesi %i %i\n",start,sz);
 size+=1000; aux=malloc(u8((size/10)+72)); B=A;
 for (i=0;i<10;i++) {B=B+k; k=get_primes_loop(start,sz/10,B); start+=sz/10;}
 free(aux);}

void init_primes()
{int i; PRIMES=malloc(80000*sizeof(int));
 for (i=0;i<80000;i++) PRIMES[i]=0; for (i=0;i<179;i++) PRIMES[i]=my_primes[i];
 get_primes(1064,1000000,PRIMES+179);
 badprimes=malloc(32*sizeof(llint)); badprimetype=malloc(32*sizeof(llint));}

static void add_to_list(LIST *L,llint P,int E)
{int i=0; if (DEBUG>=2) printf("add_to_list %lli %i\n",P,E);
 while ((*L).p[i]!=0) i++; (*L).p[i]=P; (*L).e[i]=E; (*L).p[i+1]=0;}

static int is_power(llint x,int n,llint *y)
{int i; llint t=1; *y=Round(Root(x,n));
 for (i=0;i<n;i++) t*=*y; if (t==x) return 1; return 0;}

static int is_square(int x)
{int y=(int) Round(Sqrt((double) x)); if (y*y==x) return y; return 0;}

static void squfof_internal(llint n,llint *y1,llint *y2)
{int Pnow,Qnow,Qprev,Pprev,Qback2,quot; int nS,den_bound,den,j,root;
 int phase,nsmallden,iter,smalldens[256]; /* From A. Steel and A. Lenstra */
 if (DEBUG>=2) printf("squfof_internal %lli\n",n);
 if (n<=0) errorit("Non-positive number in squfof");
 if (n>(((llint) 1)<<59)) errorit("Too large in squfof");
 nS=Pnow=(int) Floor(Sqrt((double) n));
 Qnow=(int) (n-(llint) Pnow*(llint) Pnow); Qprev=1;
 phase=1; nsmallden=0; den_bound=(int) Sqrt(2*nS+1);
 for (iter=1;;iter++)
 {Pprev=Pnow; Qback2=Qprev; Qprev=Qnow;
  if (nS+Pprev>=2*Qprev)
  {quot=(nS+Pprev)/Qprev; Pnow=quot*Qprev-Pprev;
   Qnow=Qback2+quot*(Pprev-Pnow);}
  else {Pnow=Qprev-Pprev; Qnow=Qback2+Pprev-Pnow;}
  den=Qprev; if ((den&1)==0) den>>=1;
  if (Pnow==Pprev) {(*y1)=den; (*y2)=n/(llint) den; return;}
  if (phase && den<den_bound && den>1) smalldens[++nsmallden]=den;
  if (phase && (iter&1) && ((root=is_square(Qnow))>0))
  {j=nsmallden; smalldens[0]=root; while (smalldens[j]!=root) j--;
   if (j==0)
   {phase=0; Pnow=nS-(nS-Pnow)%root;
    Qnow=(int) ((n-(llint) Pnow*(llint) Pnow)/(llint) root); Qprev=root;}}}}

static void squfof(llint n,llint *y1,llint *y2)
{int i=1; llint m=1; if (DEBUG>=2) printf("squfof %lli\n",n);
 while (1) 
 {squfof_internal(n*m,y1,y2);
  if (((*y1)%m)==0) (*y1)/=m; if (((*y2)%m)==0) (*y2)/=m;
  if (((*y1)!=1) && ((*y2)!=1)) return; else m=(llint) my_primes[i++];}}

static short int *IFACT_TABLE;
static void TABLE_FACTOR(int x,int s,LIST *L)
{int p,l=0,e=1,f=0; if (DEBUG>=2) printf("TABLE_FACTOR %i %i\n",x,s);
 if (!(x&1)) {do {f++; x=x>>1;} while (!(x&1)); add_to_list(L,2,f*s);}
 while (1)
 {p=(int) IFACT_TABLE[(x-1)/2];
  if (p==0)
  {if (x==l) {e++; x=1;} if (l>0) add_to_list(L,l,e*s);
   if (x>1) add_to_list(L,x,s); return;}
  if (p==l) {x=x/p; e++;}
  else {if (l>0) add_to_list(L,l,e*s); e=1; x=x/p; l=p;}}}

static int primes_up_to(int B) {int i; for (i=0;PRIMES[i]<=B;i++); return i;}

void modular_exponentiation(llint B,QD n,llint pow,QD M)
{QD b; llint ll,l; if (DEBUG>=2) printf("mod_exp %lli %f %lli\n",B,n[0],pow);
 QD_copy(wmax,QD_zero,b); b[0]=(double) B; QD_copy(wmax,QD_one,M); ll=l=pow;
 for (l>>=1;l>0;l>>=1)
 {QD_sqr(wmax,b,b); QD_mod(b,n,b);
  if (l&1) {QD_mul(wmax,M,b,M); QD_mod(M,n,M);}}
 if (ll&1) {QD_mul1(wmax,M,(double) B,M); QD_mod(M,n,M);}}

static int mod_test(int B,QD n,llint x)
{QD M; modular_exponentiation(B,n,x-1,M); if (M[0]==1.0) return 1; return 0;}

static int pp_test_QD(QD n,llint x)
{if (n[0]<1.0e12)
 {if (!mod_test(2,n,x)) return 0; if (!mod_test(13,n,x)) return 0;
  if (!mod_test(23,n,x)) return 0; if (!mod_test(1662803,n,x)) return 0;
  return 1;}
 if (!mod_test(2,n,x)) return 0; if (!mod_test(3,n,x)) return 0;
 if (!mod_test(5,n,x)) return 0; if (!mod_test(7,n,x)) return 0;
 if (!mod_test(11,n,x)) return 0; if (!mod_test(13,n,x)) return 0;
 if (!mod_test(17,n,x)) return 0; if (!mod_test(19,n,x)) return 0;
 if (!mod_test(23,n,x)) return 0; if (!mod_test(29,n,x)) return 0; return 1;}

static int pseudo_prime_test(llint x)
{QD n; n[0]=(double) x; n[1]=(double) (x-(llint) n[0]); n[2]=0.0; n[3]=0.0;
 return pp_test_QD(n,x);}

void QD_factor(QD x,LIST *L)
{llint f,y1,y2; int i,e,s=1; QD y,z; QD_copy(wmax,x,y);
 for (i=0;PRIMES[i]!=0;i++)
 {e=0; while ((QD_is_divisible(y,(double) PRIMES[i])) && (y[0]!=1.0))
  {e++; QD_div1(wmax,y,(double) PRIMES[i],y); QD_round(wmax,y,y);}
  if (e>0) add_to_list(L,PRIMES[i],e);}
 if (y[0]==1.0) return;
 while (QD_is_power(y,2,z)) {s*=2; QD_copy(wmax,z,y);}
 while (QD_is_power(y,3,z)) {s*=3; QD_copy(wmax,z,y);}
 while (QD_is_power(y,5,z)) {s*=5; QD_copy(wmax,z,y);}
 while (QD_is_power(y,7,z)) {s*=7; QD_copy(wmax,z,y);}
 if (y[0]>576460752303423488.0) errorit("Disc has a large prime factor");
 f=(llint) y[0]+(llint) y[1];
 if (pseudo_prime_test(f)) {add_to_list(L,f,s); return;}
 squfof(f,&y1,&y2); add_to_list(L,y1,s); add_to_list(L,y2,s);}

static int TABLE_BOUND=0;

void IFACT_INIT(llint B)
{int i,j,num,p; if (DEBUG>=2) printf("IFACT_INIT %lli\n",B);
 if (B>MAX_TABLE) B=MAX_TABLE; if (B<=TABLE_BOUND) return;
 if (TABLE_BOUND!=0) free(IFACT_TABLE); TABLE_BOUND=B;
 IFACT_TABLE=(short int*) malloc((B/2)*sizeof(short int));
 for (i=0;i<B/2;i++) IFACT_TABLE[i]=(short int) 0;
 num=primes_up_to((int) Sqrt((double) B));
 for (i=num-1;i>0;i--)
 {p=PRIMES[i]; for (j=(p*p-1)/2;j<B/2;j+=p) IFACT_TABLE[j]=(short int) p;}}

static llint TRIAL_DIV(llint x,LIST *L,int s,int B)
{int i,e; llint p; if (DEBUG) printf("TRIAL_DIV %lli %i\n",x,B);
 if (x==1) return 0;
 for (i=0;PRIMES[i]<B;i++)
 {p=(llint) PRIMES[i];
  if ((x%p)==0)
  {x=x/p; e=1; while ((x%p)==0) {e++; x=x/p;}
   add_to_list(L,p,e*s); if (x==1) return 0;}}
 return x;}

void ifactor(llint x,LIST *L,int s,int B)
{llint y,y1,y2; if (x<=0) errorit("Nonpositive number in ifactor");
 if (DEBUG>=2) {printf("ifactor %lli\n",x); fflush(stdout);}
 if (x<TABLE_BOUND) {TABLE_FACTOR((int) x,s,L); return;}
 x=TRIAL_DIV(x,L,s,B); if (!x) return;
 if (x<TABLE_BOUND) {TABLE_FACTOR((int) x,s,L); return;}
 if (pseudo_prime_test(x)) {add_to_list(L,x,s); return;}
 if (is_power(x,2,&y)) {ifactor(y,L,2*s,0); return;}
 if (is_power(x,3,&y)) {ifactor(y,L,3*s,0); return;}
 if (is_power(x,5,&y)) {ifactor(y,L,5*s,0); return;}
 squfof(x,&y1,&y2); ifactor(y1,L,s,0); ifactor(y2,L,s,0);}
