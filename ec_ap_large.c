#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

typedef struct {double x,y; int i;} ap_pt;

#define QD_split(a,hi,lo)\
{double xt=134217729.0*a,t1; t1=xt-a; hi=xt-t1; lo=a-hi;}
#define QD_prod(a,b,d,e)\
{double ah,al,bh,bl,u1,u2,u3; d=a*b; QD_split(a,ah,al); QD_split(b,bh,bl);\
 u1=ah*bh; u2=u1-d; u3=ah*bl; u1=u2+u3; u2=al*bh; u3=u1+u2; u2=al*bl; e=u2+u3;}

#define BIG_FLOAT 6755399441055744.0 /* 2^52+2^51 */
#define ROUND(x,r) {double y=(x),w1;  w1=y+BIG_FLOAT; r=w1-BIG_FLOAT;}

#define MUL(x,y,z,p)\
{double zhi,zlo,thi,tlo,fl,v1,v2; QD_prod((x),(y),zhi,zlo); ROUND(zhi/p,fl);\
 QD_prod(p,fl,thi,tlo); v1=zhi-thi; v2=zlo-tlo; z=v1+v2;}
#define INV(x,i,p)\
{double rr=0.0,ss=1.0,cc,aa=p,bb=x; while (1)\
 {ROUND(aa/bb,cc); aa-=bb*cc;\
  if (aa==0.0) {if (bb>0.0) i=ss; else i=-ss; break;} rr-=cc*ss;\
  ROUND(bb/aa,cc); bb-=aa*cc;\
  if (bb==0.0) {if (aa>0.0) i=rr; else i=-rr; break;} ss-=cc*rr;}}
#define NORMALISE(x,p) while ((x)>=(p)) (x)-=(p); while ((x)<0) (x)+=(p)

#define FROM_LINE(P,Q,R,s,p)\
{double T,x1; ap_pt V;\
 MUL(s,s,T,p); V.x=T-(P.x+Q.x); NORMALISE(V.x,p);\
 x1=P.x-V.x; MUL(x1,s,T,p); V.y=T-P.y; NORMALISE(V.y,p); *(R)=V;}

static int ec_ap_large_pt_double(double a4,ap_pt P,ap_pt *R,double p)
{double T,n,i,s,tt;
 if (P.y==0.0) return TRUE; MUL(P.x,P.x,T,p); n=T+T+T+a4; /*NORMALISE(n,p);*/
 tt=P.y+P.y; INV(tt,i,p); MUL(n,i,s,p); FROM_LINE(P,P,R,s,p); return FALSE;}

static int ec_ap_large_pt_add(double a4,ap_pt P,ap_pt Q,ap_pt *R,double p)
{double Ld,Ln,i,s;
 if (P.x==Q.x)
 {if (P.y==Q.y) return ec_ap_large_pt_double(a4,P,R,p); return TRUE;}
 Ld=P.x-Q.x; Ln=P.y-Q.y; /* allow Ld to be negative */
 INV(Ld,i,p); MUL(Ln,i,s,p); FROM_LINE(P,Q,R,s,p); return FALSE;}

#define ADDI(P,Q,R,p,i) {double s; MUL(P.y-Q.y,i,s,p); FROM_LINE(P,Q,R,s,p);}
#define DOUBLEI(a4,P,R,p,i)\
{double TT,nn,s; MUL(P.x,P.x,TT,p); nn=TT+TT+TT+a4;\
 MUL(nn,i,s,p); FROM_LINE(P,P,R,s,p);}

static int ec_ap_large_pt_mult(double a4,ap_pt P,llint n,ap_pt *R,double p)
{llint m=n; ap_pt Q=P,T,U; int z=TRUE; /* P is not O */
 if (n==0) return TRUE; if (n==1) {*R=P; return FALSE;}
 if (m&1) {U=Q; z=FALSE;}
 while (m)
 {m=m>>1;
  if (ec_ap_large_pt_double(a4,Q,&T,p)) break; Q=T;
  if (m&1)
  {if (z) {U=Q; z=FALSE;}
   else {z=ec_ap_large_pt_add(a4,Q,U,&T,p); if (!z) U=T;}}}
 if (!z) *R=U; return z;}

static llint ec_ap_large_order_from_mult(double a4,ap_pt P,llint m,double pr)
{int i,e,k,t; llint o=1,n=m,p,pe; ap_pt Q=P,R; LIST f;
 f.p[0]=0; ifactor(n,&f,1,1000);
 for (i=0;f.p[i]!=0;i++)
 {p=f.p[i]; e=f.e[i]; for (t=0;t<e;t++) n=n/p;
  if (ec_ap_large_pt_mult(a4,Q,n,&R,pr)) continue;
  for (k=1;k<e;k++) {if (ec_ap_large_pt_mult(a4,R,p,&R,pr)) break;}
  e=k; pe=1; for (t=0;t<e;t++) pe=pe*p; o=o*pe;
  if (ec_ap_large_pt_mult(a4,Q,pe,&R,pr)) break; Q=R;}
 return o;}

#define TYPE ap_pt
#define NAME ec_ap_large_qsort
#define CMP(A,B) (A->x>B->x)
#define SWAP_VARS(type,A,B) {type _ttemmp; _ttemmp=A; A=B; B=_ttemmp;}
#define SWAP_DATA(A,B) {SWAP_VARS(double,A->x,B->x);\
 SWAP_VARS(double,A->y,B->y); SWAP_VARS(int,A->i,B->i);}

static void NAME(TYPE *ARR,int n)
{TYPE *left,*right,*ip,*jp,*last,*mid; int l;
 if (DEBUG) printf("ec_ap_large_qsort %i\n",n);
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
 if (DEBUG) printf("heap_search %i %f %f\n",HH,A[HH].x,n);\
 if ((HH<sz) && (A[HH].x==n)) k=HH; else k=-1;}
#define DONE(x) {free(ARR); return ((llint) (x));}
#define DONEBS(x) {free(B); free(S); DONE((x));}
#define MONTGOMERY_INVERSE_BABY_STEPS TRUE

static llint ec_ap_large_bsgs(double a4,ap_pt P,double p,llint L)
{int i,j,k,bs,ns,gs,G,n; ap_pt *ARR,BS=P,GS,T,U; double *B,*S,u,t,x; llint K;
 if (DEBUG) printf("ec_ap_large_bsgs a4:%f x:%f y:%f p:%f\n",a4,P.x,P.y,p);
 ns=1+2*(int) Sqrt(4.0*p); bs=(int) ((1.0+Sqrt((double) ns))/2.0);
 ARR=malloc(bs*sizeof(ap_pt));
 NORMALISE(BS.y,p); NORMALISE(BS.x,p); ARR[0]=BS; ARR[0].i=1;
 if (!MONTGOMERY_INVERSE_BABY_STEPS)
 {for (i=1;i<bs;i++)
  {if (ec_ap_large_pt_add(a4,BS,ARR[i-1],&(ARR[i]),p)) DONE(i+1);
   if (ARR[i].y==0.0) DONE(2*i+2); ARR[i].i=i+1;}}
 else
 {for (i=1;i<=3;i++)
  {if (ec_ap_large_pt_add(a4,BS,ARR[i-1],&(ARR[i]),p)) DONE(i+1);}
  if (DEBUG) printf("Done with first baby steps, Total: %i\n",bs);
  B=malloc((bs/2)*sizeof(double)); S=malloc((bs/2)*sizeof(double)); n=3;
  while (n<bs-1)
  {x=ARR[n].x; G=minimum(n,bs-n-1);
   if (DEBUG) printf("Baby steps %i %f %f\n",n,ARR[n].x,ARR[n].y);
   for (i=0;i<G;i++)
   {B[i]=ARR[i].x-x; if (B[i]==0.0) DONEBS(i+n+2);}
   if (G!=bs-n-1)
   {if (ARR[n].y==0.0) DONEBS(2*n+2); B[n]=ARR[n].y+ARR[n].y; G++;}
   if (DEBUG) printf("Doing Montgomery %i\n",G);
   S[0]=B[0]; for (i=1;i<G;i++) MUL(S[i-1],B[i],S[i],p); INV(S[G-1],u,p);
   for (i=G-1;i>0;i--) {MUL(u,S[i-1],t,p); MUL(u,B[i],u,p); B[i]=t;} B[0]=u;
   if (DEBUG) printf("Montgomery done %i\n",bs-n-1);
   for (i=0;i<G;i++) ADDI(ARR[i],ARR[n],&(ARR[n+i+1]),p,B[i]);
   if (n+n+1<bs) DOUBLEI(a4,ARR[n],&(ARR[n+n+1]),p,B[n]); n=n+n+1;}
  free(B); free(S);
  for (i=0;i<bs;i++) {if (ARR[i].y==0.0) DONE(2*i+2); ARR[i].i=i+1;}}
 if (DEBUG) printf("%i baby steps\n",bs);
 if (DEBUG>=1) for (i=0;i<bs;i++) printf("bs:%i %f %f\n",i,ARR[i].x,ARR[i].y);
 gs=2*bs; if (ec_ap_large_pt_double(a4,ARR[bs-1],&GS,p)) DONE(gs);
 if (ec_ap_large_pt_mult(a4,P,L,&T,p)) DONE(L); ec_ap_large_qsort(ARR,bs);
 if (DEBUG>=1) printf("T.x:%f T.y:%f\n",T.x,T.y);
 HEAP_SEARCH(T.x,ARR,bs,k);
 if (k!=-1) {i=ARR[k].i; if (T.y==ARR[k].y) K=L-i; else K=L+i; if (K) DONE(K);}
 for (j=1;;j++)
 {if (ec_ap_large_pt_add(a4,T,GS,&U,p)) DONE((L+(gs*j)));
  T=U; HEAP_SEARCH(T.x,ARR,bs,k);
  if (k!=-1) /* Don't use the GB trick, as it probably isn't faster */
  {i=ARR[k].i; if (T.y==ARR[k].y) K=L-i; else K=L+i;
   if (DEBUG) printf("FOUND %i %lli %i %i\n",i,K,j,gs); DONE((K+(gs*j)));}}
 errorit("ec_bsgs_bsgs error"); return -1;}

llint ec_bsgs_ap_large(QD C4,QD C6,llint pr)
{llint aa,bb,v,U,L,S,WIDTH,m,o; QD AA,BB; int x; double V;
 int k,ko; double a4,p=(double) pr; ap_pt P;
 if (DEBUG) printf("ec_bsgs_ap_large %lli\n",pr);
 if (pr<500) errorit("p<500 in ec_bsgs_ap\n");
 QD_mul1(wmax,C4,-27.0,AA); aa=QD_modll(AA,p);
 QD_mul1(wmax,C6,-54.0,BB); bb=QD_modll(BB,p); v=bb;
 S=(llint) Floor(Sqrt(4*p)); L=pr+1-S; U=pr+1+S; WIDTH=U-L+1;
 if (DEBUG) printf("a4:%lli a6:%lli\n",aa,bb);
 if (DEBUG) printf("L:%lli U:%lli W:%lli\n",L,U,WIDTH);
 k=kronll(v,pr); ko=-k; x=0;
 while (1)
 {while (!k || (k==ko)) {x++; v=((aa+x*x)*x+bb)%pr; k=kronll(v,pr);}
  ko=k; if (DEBUG) printf("x:%i k:%i v:%lli\n",x,k,v); P.x=(double) ((x*v)%pr);
  V=(double) v; MUL(V,V,P.y,p); MUL(P.y,(double) aa,a4,p); if (a4+a4>p) a4-=p;
  if (DEBUG) printf("a4:%f Px:%f Py:%f\n",a4,P.x,P.y);
  m=ec_ap_large_bsgs(a4,P,p,L); if (DEBUG) printf("mult of order is %lli\n",m);
  o=ec_ap_large_order_from_mult(a4,P,m,p); if (DEBUG) printf("ord %lli\n",o);
  if (o<WIDTH) continue; m=o*(L/o); if (m!=L) m+=o; m-=(pr+1);
  if (k==1) m=-m; return m;}}
