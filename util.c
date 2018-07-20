#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

const char *VERBOSE2option[3]={"-quiet","-terse","-verbose"};

void errorit(S) char *S; {printf("**ERROR** %s\n",S); exit(-1);}

int u8(int x) {return(((x+7)>>3)<<3);}

static unsigned char *auxp;

static int get_primesll_loop(llint start,llint size,llint *A)
{llint i,l,j=0,k=0,s,p=2;
 if (DEBUG>=2) printf("get_primesll_loop %lli %lli\n",start,size);
 for (i=0;i<size;i++) auxp[i]=1;
 while(p<=Sqrt(start+size+100))
 {s=start-p*(start/p); if (s!=0) l=p-s; else l=0;
  for (i=l;i<size;i+=p) auxp[i]=0; p=PRIMES[j++];}
 for (i=0;i<size;i++) {if (auxp[i]) {A[k++]=start+i;}} return k;}

void get_primes_ll(llint st,llint sz,llint *A)
{int i,k=0; llint *B; if (DEBUG>=2) printf("get_primes_ll %lli %lli\n",st,sz);
 sz+=1000; auxp=malloc(u8((sz/10)+72)); B=A;
 for (i=0;i<10;i++) {B=B+k; k=get_primesll_loop(st,sz/10,B); st+=sz/10;}
 if (DEBUG>=2) printf("Done with get_primes_ll\n"); free(auxp);}

void free_data()
{QD **P; int i; if (DEBUG) printf("free_data\n");
 P=TABLE[0]; for (i=0;i<MESH_COUNT[0];i++) free(P[i]); free(P);
 P=POWSER[0]; for (i=0;i<=NUM_LOGS[0];i++) free(P[i]); free(P);}

int kron(int a,int b) /* assumes b is odd and positive */
{int e=0; a=a%b; if (a<0) a+=b; if (a<=1) return a; if (b==1) return 1;
 while ((a&1)==0) {e++; a>>=1;}
 if ((e&1) && (((b&7)==3) || ((b&7)==5))) e=-1; else e=1; if (a==1) return e;
 if (((a&3)==3) && ((b&3)==3)) return -e*kron(b,a); return e*kron(b,a);}

int kronll(llint a,llint b) /* assumes b is odd and positive */
{int e=0; a=a%b; if (a<0) a+=b; if (a<=1) return a; if (b==1) return 1;
 while ((a&1)==0) {e++; a>>=1;}
 if ((e&1) && (((b&7)==3) || ((b&7)==5))) e=-1; else e=1; if (a==1) return e;
 if (((a&3)==3) && ((b&3)==3)) return -e*kronll(b,a); return e*kronll(b,a);}

int gcd(int a,int b)
{while (1) {a-=(a/b)*b; if (a==0) return b; b-=(b/a)*a; if (b==0) return a;}}

llint gcdll(llint a,llint b)
{while (1) {a-=(a/b)*b; if (a==0) return b; b-=(b/a)*a; if (b==0) return a;}}
