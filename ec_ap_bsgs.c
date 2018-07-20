#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

typedef struct {int x,y; int i;} AP_pt;

#define MUL(x,y,z,p) z=(((llint) x)*((llint) y))%p
#define INV(x,i,p)\
{int rr=0,ss=1,cc,aa=p,bb=x; while (1)\
 {cc=aa/bb; aa-=bb*cc;\
  if (aa==0) {if (bb>0) i=ss; else i=-ss; break;} rr-=cc*ss;\
  cc=bb/aa; bb-=aa*cc;\
  if (bb==0) {if (aa>0) i=rr; else i=-rr; break;} ss-=cc*rr;}}
#define NORMALISE(x,p) while ((x)>=(p)) (x)-=(p); while ((x)<0) (x)+=(p)

#define FROM_LINE(P,Q,R,s,p)\
{int T,x1; AP_pt V;\
 MUL(s,s,T,p); V.x=T-(P.x+Q.x); NORMALISE(V.x,p);\
 x1=P.x-V.x; MUL(x1,s,T,p); V.y=T-P.y; NORMALISE(V.y,p); *(R)=V;}

static int ec_AP_pt_double(int a4,AP_pt P,AP_pt *R,int p)
{int T,n,i,s,tt; /* NORMALISE(P.x,p); NORMALISE(P.y,p); */
 if (P.y==0) return TRUE; MUL(P.x,P.x,T,p); n=T+T+T+a4; NORMALISE(n,p);
 tt=P.y+P.y; INV(tt,i,p); MUL(n,i,s,p); FROM_LINE(P,P,R,s,p); return FALSE;}

static int ec_AP_pt_add(int a4,AP_pt P,AP_pt Q,AP_pt *R,int p)
{int Ld,Ln,i,s;  /* NORMALISE(P.x,p); NORMALISE(P.y,p); etc. */
 if (P.x==Q.x)
 {if (P.y==Q.y) return ec_AP_pt_double(a4,P,R,p); return TRUE;}
 Ld=P.x-Q.x; Ln=P.y-Q.y; /* allow Ld to be negative */
 INV(Ld,i,p); MUL(Ln,i,s,p); FROM_LINE(P,Q,R,s,p); return FALSE;}

#define ADDI(P,Q,R,p,i) {int s; MUL(P.y-Q.y,i,s,p); FROM_LINE(P,Q,R,s,p);}
#define DOUBLEI(a4,P,R,p,i)\
{int TT,nn,s; MUL(P.x,P.x,TT,p); nn=TT+TT+TT+a4;\
 MUL(nn,i,s,p); FROM_LINE(P,P,R,s,p);}

static int ec_AP_pt_mult(int a4,AP_pt P,int n,AP_pt *R,int p)
{int m=n; AP_pt Q=P,T,U; int z=TRUE; /* P is not O */
 if (n==0) return TRUE; if (n==1) {*R=P; return FALSE;}
 if (m&1) {U=Q; z=FALSE;}
 while (m)
 {m=m>>1; if (ec_AP_pt_double(a4,Q,&T,p)) break; Q=T;
  if (m&1)
  {if (z) {U=Q; z=FALSE;} else {z=ec_AP_pt_add(a4,Q,U,&T,p); if (!z) U=T;}}}
 if (!z) *R=U; return z;}

static int ec_AP_pt_order_from_mult(int a4,AP_pt P,int m,int pr)
{int i,e,k,t; int o=1,n=m,p,pe; AP_pt Q=P,R; LIST f;
 f.p[0]=0; ifactor(n,&f,1,1000); NORMALISE(Q.x,pr); NORMALISE(Q.y,pr);
 for (i=0;f.p[i]!=0;i++)
 {p=f.p[i]; e=f.e[i]; for (t=0;t<e;t++) n=n/p;
  if (ec_AP_pt_mult(a4,Q,n,&R,pr)) continue;
  for (k=1;k<e;k++) {if (ec_AP_pt_mult(a4,R,p,&R,pr)) break;}
  e=k; pe=1; for (t=0;t<e;t++) pe=pe*p; o=o*pe;
  if (ec_AP_pt_mult(a4,Q,pe,&R,pr)) break; Q=R;}
 return o;}

#define TYPE AP_pt
#define NAME ec_AP_qsort
#define CMP(A,B) (A->x>B->x)
#define SWAP_VARS(type,A,B) {type _ttemmp; _ttemmp=A; A=B; B=_ttemmp;}
#define SWAP_DATA(A,B) {SWAP_VARS(int,A->x,B->x);\
 SWAP_VARS(int,A->y,B->y); SWAP_VARS(int,A->i,B->i);}

static void NAME(TYPE *ARR,int n)
{TYPE *left,*right,*ip,*jp,*last,*mid; int l;
 if (DEBUG) printf("ec_AP_qsort %i\n",n);
 if (n<=1) return; left=ARR;
 if (n==2) {ip=left+1; if (CMP(left,ip)) SWAP_DATA(left,ip); return;}
 right=left+(n-1);
 if (n<=5)
 {for (ip=left;ip<right;ip++)
  for (jp=ip+1;jp<=right;jp++) if (!CMP(jp,ip)) SWAP_DATA(ip,jp); return;}
 mid=left+(n/2);
 if (!CMP(left,mid))
 {if (!CMP(mid,right)) SWAP_DATA(left,mid)
  else if (!CMP(left,right)) SWAP_DATA(left,right)}
 else
 {if (CMP(mid,right)) SWAP_DATA(left,mid)
  else if (CMP(left,right)) SWAP_DATA(left,right)}
 l=0; last=left;
 for (ip=left+1;ip<=right;ip++)
   if (!CMP(ip,left)) {l++; last++; SWAP_DATA(last,ip)}
 if (left!=last) SWAP_DATA(left,last)
  NAME(left,l); NAME(last+1,n-l-1);}

#define HEAP_SEARCH(n,A,sz,k)\
{int LL=0,HH=sz,w;\
 while (LL<HH) {w=(LL+HH)>>1; if (A[w].x<n) LL=w+1; else HH=w;}\
 if (DEBUG>=2) printf("heap_search %i %i %i\n",HH,A[HH].x,n);\
 if ((HH<sz) && (A[HH].x==n)) k=HH; else k=-1;}
#define DONE(x) {free(ARR); return x;}
#define DDONE(x) {if (mo==0) mo=(x); else {free(ARR); return gcd(mo,x);}}
#define DONEBS(x) {free(B); free(S); DONE(x);}
#define MONTGOMERY_INVERSE_BABY_STEPS TRUE

static llint ec_AP_bsgs(int a4,AP_pt P,int p,int L)
{int i,j,k,bs,ns,gs,G,n; AP_pt *ARR,BS=P,GS,T,U; int *B,*S,u,t,x,mo=0;
 if (DEBUG) printf("ec_AP_bsgs a4:%i x:%i y:%i p:%i\n",a4,P.x,P.y,p);
 ns=1+2*(int) Sqrt(4.0*(double) p); bs=(int) ((1.0+Sqrt((double) ns))/2.0);
 ARR=malloc(bs*sizeof(AP_pt)); ARR[0]=BS; ARR[0].i=1;
 if (!MONTGOMERY_INVERSE_BABY_STEPS)
 {for (i=1;i<bs;i++)
  {if (ec_AP_pt_add(a4,BS,ARR[i-1],&(ARR[i]),p)) DONE(i+1);
   if (ARR[i].y==0) DONE(2*i+2); ARR[i].i=i+1;}}
 else
 {for (i=1;i<=3;i++)
  {if (ec_AP_pt_add(a4,BS,ARR[i-1],&(ARR[i]),p)) DONE(i+1);}
  if (DEBUG) printf("Done with first baby steps, Total: %i\n",bs);
  B=malloc((bs/2)*sizeof(int)); S=malloc((bs/2)*sizeof(int)); n=3;
  while (n<bs-1)
  {x=ARR[n].x; G=minimum(n,bs-n-1);
   if (DEBUG) printf("Baby steps %i %i %i\n",n,ARR[n].x,ARR[n].y);
   for (i=0;i<G;i++)
   {B[i]=ARR[i].x-x; if (B[i]==0) DONEBS(i+n+2);}
   if (G!=bs-n-1)
   {if (ARR[n].y==0) DONEBS(2*n+2); B[n]=ARR[n].y+ARR[n].y; G++;}
   if (DEBUG) printf("Doing Montgomery %i\n",G);
   S[0]=B[0]; for (i=1;i<G;i++) MUL(S[i-1],B[i],S[i],p); INV(S[G-1],u,p);
   for (i=G-1;i>0;i--) {MUL(u,S[i-1],t,p); MUL(u,B[i],u,p); B[i]=t;} B[0]=u;
   if (DEBUG) printf("Montgomery done %i\n",bs-n-1);
   for (i=0;i<G;i++) ADDI(ARR[i],ARR[n],&(ARR[n+i+1]),p,B[i]);
   if (n+n+1<bs) DOUBLEI(a4,ARR[n],&(ARR[n+n+1]),p,B[n]); n=n+n+1;}
  free(B); free(S);
  for (i=0;i<bs;i++) {if (ARR[i].y==0) DONE(2*i+2); ARR[i].i=i+1;}
 }
 if (DEBUG) printf("%i baby steps\n",bs);
 if (DEBUG>=1) for (i=0;i<bs;i++) printf("bs:%i %i %i\n",i,ARR[i].x,ARR[i].y);
 gs=2*bs; if (ec_AP_pt_double(a4,ARR[bs-1],&GS,p)) DONE(gs);
 if (ec_AP_pt_mult(a4,P,L,&T,p)) DONE(L); ec_AP_qsort(ARR,bs);
 if (DEBUG>=1) printf("T.x:%i T.y:%i\n",T.x,T.y);
 HEAP_SEARCH(T.x,ARR,bs,k);
 if (k!=-1) {i=ARR[k].i; if (T.y==ARR[k].y) k=L-i; else k=L+i; if (k) DONE(k);}
 for (j=1;j*(2*bs-1)<=ns;j++)
 {if (ec_AP_pt_add(a4,T,GS,&U,p)) {DDONE(L+gs*j); j++; T=GS;}
  else T=U; HEAP_SEARCH(T.x,ARR,bs,k);
  if (k!=-1)  /* idea of Geoff Bailey; avoid computing exact order */
  {i=ARR[k].i; if (T.y==ARR[k].y) k=L-i; else k=L+i;
   if (DEBUG) printf("FOUND %i %i %i %i\n",i,k,j,gs); DDONE(k+gs*j);}}
 ASSERT(mo!=0); free(ARR); return -mo;}

llint ec_bsgs_ap(QD C4,QD C6,llint pr)
{int aa,bb; if (DEBUG) printf("ec_bsgs_ap %lli\n",pr);
 if (pr>(1<<29)) return ec_bsgs_ap_large(C4,C6,pr);
 aa=(int) ((-27*QD_modll(C4,pr))%pr); bb=(int) ((-54*QD_modll(C6,pr))%pr);
 return ec_bsgs_ap_AB(aa,bb,pr);}

llint ec_bsgs_ap_AB(int aa,int bb,llint pr)
{int a4,L,U,WIDTH,v,S,k,ko,x,m,o; AP_pt P;
 S=(int) Floor(2.0*Sqrt(pr)); L=(int) pr+1-S; U=(int) pr+1+S; WIDTH=U-L+1;
 if (DEBUG) printf("a4:%i a6:%i\n",aa,bb);
 if (DEBUG) printf("L:%i U:%i W:%i\n",L,U,WIDTH);
 v=bb; NORMALISE(aa,pr); NORMALISE(v,pr); k=kron(v,(int) pr); ko=-k; x=0;
 while (1)
 {while (!k || (k==ko))
  {x++; v=(int) ((((llint) aa+(llint) (x*x))*(llint) x+(llint) bb)%pr);
   k=kron(v,pr);}
  ko=k; if (DEBUG) printf("x:%i k:%i v:%i\n",x,k,v);
  MUL(x,v,P.x,pr); MUL(v,v,P.y,pr); MUL(aa,P.y,a4,pr);
  NORMALISE(P.x,pr); NORMALISE(P.y,pr); NORMALISE(a4,pr);
  if (DEBUG) printf("a4:%i Px:%i Py:%i\n",a4,P.x,P.y);
  m=ec_AP_bsgs(a4,P,pr,L); if (DEBUG) printf("multiple of order is %i\n",m);
  if (m>0) {if (m<WIDTH) continue; o=ec_AP_pt_order_from_mult(a4,P,m,pr);}
  else o=-m; if (DEBUG) printf("ord is %i\n",o);
  if (o<WIDTH) continue; m=o*(L/o); if (m!=L) m+=o; m-=((int) pr+1);
  if (k==1) m=-m; return (llint) m;}}
