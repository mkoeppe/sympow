#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

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

typedef struct {llint U; llint V;} PERALTA;

static void PMUL(PERALTA i1,PERALTA i2,PERALTA *o,llint d,llint p)
{double A,B,P=(double) p,D=(double) d,E,x,y; llint rU,rV;
 if (DEBUG>=2) printf("PMUL %lli %lli %lli %lli\n",i1.U,i1.V,i2.U,i2.V);
 x=(double) i1.V; y=(double) i2.U; MUL(x,y,A,P);
 x=(double) i2.V; y=(double) i1.U; MUL(x,y,B,P);
 rV=(llint) (A+B); while (rV<0) rV+=p; while (rV>=p) rV-=p;
 x=(double) i2.U; y=(double) i1.U; MUL(x,y,A,P);
 x=(double) i1.V; y=(double) i2.V; MUL(x,y,B,P); MUL(B,D,E,P);
 rU=(llint) (A+E); while (rU<0) rU+=p; while (rU>=p) rU-=p;
 (*o).U=rU; (*o).V=rV; if (DEBUG>=2) printf("RET %lli %lli\n",rU,rV);}

static void mod_exp_peralta(llint r,llint d,llint p,llint e,llint* u,llint *v)
{PERALTA c,m,o; llint le; o.U=r; o.V=1; m.U=1; m.V=0; c=o;
 le=e; for (e>>=1;e>0;e>>=1) {PMUL(c,c,&c,d,p); if (e&1) PMUL(m,c,&m,d,p);}
 if (le&1) PMUL(m,o,&m,d,p); (*u)=m.U; (*v)=m.V;}

static llint mod_sqrt_peralta(int d,llint p)
{llint r=314159,u,v; QD P,R,S; double V,i;
 if (DEBUG>=2) printf("mod_sqrt_peralta %i %lli\n",d,p);
 QD_copy(wmax,QD_zero,P); P[0]=(double) p;
 while (1)
 {r=(r+1)%p; R[0]=(double) r; QD_sqr(wmax,R,S);
  if (QD_modll(S,p)==p-d) return r; mod_exp_peralta(r,d,p,(p-1)/2,&u,&v);
  if (u==0) {V=(double) v; INV(V,i,P[0]); return (llint) i;}}}

static llint mod_sqrt(int d,llint p)
{QD P,M; QD_copy(wmax,QD_zero,P); P[0]=(double) p;
 if ((p&3)==3)
 {modular_exponentiation(d,P,(p+1)/4,M); return (llint) Round(M[0]);}
 else return mod_sqrt_peralta(d,p);}

static void cornaccia2(int d,llint p,llint *x,llint *y) /* x^2+d*y^2=4p */
{llint a,b,L,r; if (DEBUG>=2) printf("cornaccia2 %i %lli\n",d,p);
 a=p+p; b=mod_sqrt(-d,p); if ((b&1)!=(d&1)) b-=p;
 L=(llint) Floor(2.0*Sqrt((double) p));
 while ((b>L) || (b<-L)) {r=a%b; a=b; b=r;}
 *y=(llint) Round(Sqrt((((p<<2)-b*b)/d))); *x=b;
 if (*x<0) *x=-(*x); if (*y<0) *y=-(*y);}

static llint ap_j0(QD C6,llint pr)
{llint x,y,d,ap; QD P,R; if ((pr%3)!=1) return 0; cornaccia2(27,pr,&x,&y);
 if ((x%3)==1) x=-x;
 d=(8*QD_modll(C6,pr))%pr; QD_copy(wmax,QD_zero,P); P[0]=(double) pr;
 modular_exponentiation(d,P,(pr-1)/6,R);
 P[0]=(double) x; QD_mul(wmax,R,P,R); ap=QD_modll(R,pr);
 if (ap>pr/2) ap-=pr; return ap;}

static llint ap_j1728(QD C4,llint pr)
{llint x,y,d,ap; QD P,R; if ((pr&3)!=1) return 0; cornaccia2(4,pr,&x,&y);
 if ((x&3)==0) x=y; if ((x&1)==1) x=(x<<1); if ((x&7)==6) x=-x;
 d=(-27*QD_modll(C4,pr))%pr; QD_copy(wmax,QD_zero,P); P[0]=(double) pr;
 modular_exponentiation(d,P,(pr-1)/4,R);
 P[0]=(double) x; QD_mul(wmax,R,P,R); ap=QD_modll(R,pr);
 if (ap>pr/2) ap-=pr; return ap;}

static llint ap_j8000(llint pr)
{llint x,y; if (((pr&7)==5) || ((pr&7)==7)) return 0; cornaccia2(8,pr,&x,&y);
 if (((x&15)==2) || ((x&15)==6)) {if ((y&3)!=0) x=-x;}
 if (((x&15)==10) || ((x&15)==14)) {if ((y&3)==0) x=-x;}
 if (kronll(CM_TWIST,pr)==-1) x=-x; return x;}

static llint ap_j287496(llint pr)
{llint x,y; if ((pr&3)!=1) return 0; cornaccia2(4,pr,&x,&y);
 if ((x&3)==0) x=y; if ((x&1)==1) x=(x<<1); if ((x&7)==6) x=-x;
 if (kronll(2,pr)==-1) x=-x; if (kronll(CM_TWIST,pr)==-1) x=-x; return x;}

static llint ap_cm(int CM,llint pr)
{llint x,y; if (kronll((llint) CM,pr)!=1) return 0; cornaccia2(-CM,pr,&x,&y);
 if ((CM&3)==0) CM>>=2; if ((kronll(x,(llint) -CM)>0)^(CM==-7)) x=-x;
 if (kronll(CM_TWIST,pr)==-1) x=-x; return x;}

static llint ec_ap_cm(QD C4,QD C6,llint pr)
{switch (CM_CASE)
 {case -3: return ap_j0(C6,pr); case -4: return ap_j1728(C4,pr);
  case -8: return ap_j8000(pr); case -16: return ap_j287496(pr);
  case -28: return ap_cm(-7,pr); default: return ap_cm(CM_CASE,pr);}}

static llint ec_ap2(QD C4,QD C6)
{int a1,a2,a3,a4,a6,b2,b4,b6,c4,c6,ap; if (DEBUG) printf("ec_ap2\n");
 c4=QD_modi(C4,10077696); c6=QD_modi(C6,10077696);
 b2=-(c6%12); if (b2<-5) b2+=12; b4=(b2*b2-c4)/24;
 b6=(-b2*b2*b2+36*b2*b4-c6)/216;
 a1=b2%2; if (a1<0) a1=-a1; a2=(b2-a1)/4;
 a3=(b6%2); if (a3<0) a3=-a3; a4=(b4-a1*a3)/2; a6=(b6-a3)/4;
 a2=a2%2; if (a2<0) a2=-a2; a4=a4%2; if (a4<0) a4=-a4;
 a6=a6%2; if (a6<0) a6=-a6; ap=a1;
 if (DEBUG) printf("ec_ap2: %i %i %i %i %i\n",a1,a2,a3,a4,a6);
 if (a1==a3) ap=-ap; else if (a2!=a4) ap-=2; if (a6!=0) ap=-ap; return ap;}

static llint ec_ap3(QD C4,QD C6)
{int b2,b4,b6,c4,c6; if (DEBUG) printf("ec_ap3\n");
 c4=QD_modi(C4,10077696); c6=QD_modi(C6,10077696);
 b2=-(c6%12); if (b2<-5) b2+=12; b4=(b2*b2-c4)/24;
 b6=(-b2*b2*b2+36*b2*b4-c6)/216; b2=b2%3; b4=b4%3; b6=b6%3;
 if (DEBUG) printf("ec_ap3: %i %i %i\n",b2,b4,b6);
 return -(kron(b6,3)+kron(1+b2-b4+b6,3)+kron(-1+b2+b4+b6,3));}

static llint ec_ap_kron(QD C4,QD C6,llint pr)
{llint ap=0; int a4,a6,x; if (DEBUG) printf("ec_ap_kron %lli\n",pr);
 if (pr==2) return ec_ap2(C4,C6); if (pr==3) return ec_ap3(C4,C6);
 a4=(-27*QD_modi(C4,(double) pr))%pr; a6=(-54*QD_modi(C6,(double) pr))%pr;
 for (x=0;x<pr;x++) ap+=kron(x*x*x+a4*x+a6,(int) pr); return -ap;}

llint ec_ap(QD C4,QD C6,llint pr)
{QD D,T; if (DEBUG) printf("ec_ap %lli\n",pr);
 if (pr<500) return ec_ap_kron(C4,C6,pr);
 QD_powi(wmax,C4,3,T); QD_sqr(wmax,C6,D);
 QD_sub(wmax,T,D,T); QD_div1(wmax,T,1728.0,D); QD_round(wmax,D,D);
 return ec_ap_with_disc(C4,C6,D,pr);}

llint ec_ap_with_disc(QD C4,QD C6,QD D,llint pr)
{QD T; if (QD_valuation(D,(double) pr)!=0)
 {QD_neg(wmax,C6,T); return kronll(QD_modll(T,(double) pr),pr);}
 if (CM_CASE) return ec_ap_cm(C4,C6,pr); return ec_bsgs_ap(C4,C6,pr);}
