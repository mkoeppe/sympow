#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static void generic_plusminus(llint p,int m,int A,int B,QDpoly *v)
{QD P,t; QDpoly v1,v2,w;
 if (DEBUG) printf("generic_plusminus %lli %i %i %i\n",p,m,A,B);
 QD_copy(wmax,QD_zero,P); P[0]=(double) p;
 if (m&1)
 {initQDpoly(&w,2); ASSERT(A==B);
  QD_copy(wmax,QD_one,w.coeff[0]);
  QD_powi(wmax,P,m,t); QD_copy(wmax,t,w.coeff[2]);
  QDpoly_pow(w,A,v,-1); delQDpoly(&w); return;}
 initQDpoly(&w,1); QD_copy(wmax,QD_one,w.coeff[0]);
 QD_powi(wmax,P,m/2,t); QD_copy(wmax,t,w.coeff[1]);
 QDpoly_pow(w,A,&v1,-1); QD_neg(wmax,w.coeff[1],w.coeff[1]);
 QDpoly_pow(w,B,&v2,-1); delQDpoly(&w); QDpoly_mul(v1,v2,v,-1);
 delQDpoly(&v1); delQDpoly(&v2);}

static void cyclic_nonabelian(llint p,int m,int eps,QDpoly *v)
{int A,B;
 if (DEBUG) printf("cyclic_nonabelian %lli %i %i\n",p,m,eps);
 if (m&1) {ASSERT((eps&1)==0); A=eps/2; B=eps/2;}
 else if (m&2) {A=(eps+1)/2; B=A-1;} else {A=(eps-1)/2; B=A+1;}
 return generic_plusminus(p,m,A,B,v);}

static void cyclic_abelian_ap(llint p,int ap,int m,int d,QDpoly *v)
{QDpoly qp,lp,A; QD P,t,u; int i;
 initQDpoly(&A,0); QD_copy(wmax,QD_one,A.coeff[0]);
 initQDpoly(&qp,2); QD_copy(wmax,QD_zero,P); P[0]=(double) p;
 QD_powi(wmax,P,m,t); QD_copy(wmax,t,qp.coeff[2]);
 QD_copy(wmax,QD_one,qp.coeff[0]);
 for (i=0;i<(m+1)/2;i++)
 {if (((2*i-m)%d)==0)
  {power_trace_from_ap(wmax,p,ap,m-2*i,t);
   QD_neg(wmax,t,t); QD_powi(wmax,P,i,u); QD_mul(wmax,t,u,qp.coeff[1]);
   QDpoly_mul(A,qp,v,-1); delQDpoly(&A); A=*v;}}
 delQDpoly(&qp);
 if ((m&1)==0)
 {initQDpoly(&lp,1); QD_copy(wmax,QD_one,lp.coeff[0]);
  QD_powi(wmax,P,m/2,t); QD_neg(wmax,t,t); QD_copy(wmax,t,lp.coeff[1]);
  QDpoly_mul(A,lp,v,-1); delQDpoly(&A); delQDpoly(&lp);}}

static void hecke_good(llint p,int ap,int m,QDpoly *v)
{QD T; if (m==0) errorit("m=0 in hecke_good?!");
 initQDpoly(v,2); QD_copy(wmax,QD_one,(*v).coeff[0]);
 if (ap!=0)
 {QD_copy(wmax,QD_zero,T); T[0]=(double) p; QD_powi(wmax,T,m,T);
  power_trace_from_ap(wmax,p,ap,m,(*v).coeff[1]);
  QD_neg(wmax,(*v).coeff[1],(*v).coeff[1]); QD_copy(wmax,T,(*v).coeff[2]);}
 else
 {QD_copy(wmax,QD_zero,T); T[0]=(double) -p; QD_powi(wmax,T,m,T);
  QD_neg(wmax,T,(*v).coeff[2]); QD_copy(wmax,QD_zero,(*v).coeff[1]);}}

static void cyclic_abelian(llint p,int m,int d,QDpoly *v,int bpt)
{double P=(double) p; int v4,v6,ap,i,c4,c6; QD Et4,Et6;
 if (DEBUG) printf("cyclic abelian p:%lli sp:%i Cd:%i bpt:%i\n",p,m,d,bpt);
 if ((d>2) && (p<=3)) {cyclic_abelian_ap(p,p,m,d,v); return;}
 if ((d>1) && (p>3))
 {QD_copy(wmax,QD_zero,Et4); QD_copy(wmax,QD_zero,Et6);
  v4=QD_valuation(Ec4,P); v6=QD_valuation(Ec6,P);
  if (3*v4>=2*v6)
  {QD_copy(wmax,Ec6,Et6); for (i=1;i<=v6;i++) QD_div1(wmax,Et6,P,Et6);}
  if (3*v4<=2*v6)
  {QD_copy(wmax,Ec4,Et4); for (i=1;i<=v4;i++) QD_div1(wmax,Et4,P,Et4);}}
 if ((d==2) && (p==2)) /* need to worry here */
 {if (bpt==16) {QD_div1(wmax,Ec4,16.0,Et4); QD_div1(wmax,Ec6,64.0,Et6);}
  else
  {QD_div1(wmax,Ec4,4.0,Et4); QD_div1(wmax,Ec6,8.0,Et6);
   c4=QD_modi(Et4,1<<29); c6=QD_modi(Et6,1<<29);
   if ((((c4&31)==16) && ((c6&255)==192)) ||
       (((c4&255)==0) && ((c6&2047)==512)))
   {QD_div1(wmax,Et4,16.0,Et4); QD_div1(wmax,Et6,64.0,Et6);}
   else if ((((c4&31)==16) && ((c6&255)==64)) ||
	    (((c4&255)==0) && ((c6&2047)==1536)))
   {QD_div1(wmax,Et4,16.0,Et4); QD_div1(wmax,Et6,-64.0,Et6);}
  }
 }
 if ((d==2) && (p==3)) /* think this is OK */
 {QD_div1(wmax,Ec4,9.0,Et4); QD_div1(wmax,Ec6,-27.0,Et6);}
 if (d==1) ap=(int) ec_do(p); else ap=(int) ec_ap(Et4,Et6,p);
 if ((HECKE) && (d==1)) return hecke_good(p,ap,m,v);
 cyclic_abelian_ap(p,ap,m,d,v);}

static void euler_factor_good(llint p,int m,QDpoly *v)
{cyclic_abelian(p,m,1,v,0);}

void euler_factor_bad(llint p,int bpt,int m,QDpoly *v)
{int ap,eps,lc; if (DEBUG) printf("euler_factor_bad %lli %i %i\n",p,bpt,m);
 lc=tame_local_conductor(bpt,m); eps=m+1-lc;
 if (lc==m+1) {initQDpoly(v,0); QD_copy(wmax,QD_one,(*v).coeff[0]); return;}
 if ((bpt==1) || (bpt>=29))
 {ap=ec_ap(Ec4,Ec6,p); initQDpoly(v,1); QD_copy(wmax,QD_one,(*v).coeff[0]);
  if ((ap>=0) || !(m&1)) QD_neg(wmax,(*v).coeff[0],(*v).coeff[1]);
  else QD_copy(wmax,QD_one,(*v).coeff[1]); return;}
 if ((bpt>=2) && (bpt<=6))
 {if ((p%bpt)!=1) cyclic_nonabelian(p,m,eps,v);
  else cyclic_abelian(p,m,bpt,v,bpt); return;}
 if (((bpt>=8) && (bpt<=10)) || ((bpt>=24) && (bpt<=27)))
 {if (((m&7)==2) || ((m&7)==4)) generic_plusminus(2,m,eps/2,eps/2,v);
  if ((m&7)==6) generic_plusminus(2,m,(eps+1)/2,(eps-1)/2,v);
  if ((m&7)==0) generic_plusminus(2,m,(eps-1)/2,(eps+1)/2,v); return;
 }
 if (bpt==22) {cyclic_abelian(2,m,4,v,bpt); return;}
 if (bpt==23) {cyclic_nonabelian(p,m,eps,v); return;}
 if ((bpt==12) || (bpt==13))
 {if ((eps&1)==0) generic_plusminus(3,m,eps/2,eps/2,v);
  else if ((m&3)==0) generic_plusminus(3,m,(eps-1)/2,(eps+1)/2,v);
  else if ((m&3)==2) generic_plusminus(3,m,(eps+1)/2,(eps-1)/2,v);}
 if ((bpt==16) || (bpt==17)) cyclic_abelian(2,m,2,v,bpt);
 if (bpt==14) cyclic_abelian(3,m,3,v,bpt);
 if (bpt==15) cyclic_abelian(3,m,6,v,bpt);
 if ((bpt==18) || (bpt==19)) cyclic_nonabelian(p,m,eps,v);
 if ((bpt==20) || (bpt==21)) cyclic_nonabelian(p,m,eps,v); return;}

static int deflate(int CM)
{if ((CM==-27) || (CM== -12)) return 3; if (CM==-28) return 7;
 if (CM==-16) return 4; return -CM;}

void euler_factor_hecke_bad(llint p,int bpt,int m,QDpoly *v)
{QDpoly ef1,ef2,poly,R; int i,k; QD P,T; int A[8]={0,1,0,-1,0,-1,0,1};
 if (DEBUG) printf("euler_factor_hecke_bad p:%lli bpt:%i sp:%i\n",p,bpt,m);
 QD_copy(wmax,QD_zero,P); P[0]=(double) p; euler_factor_bad(p,bpt,m,&ef1);
 if (m>2) euler_factor_bad(p,bpt,m-2,&ef2);
 else {initQDpoly(&ef2,0); QD_copy(wmax,QD_one,ef2.coeff[0]);}
 if ((m&1)==0)
 {initQDpoly(&poly,1); QD_copy(wmax,QD_one,poly.coeff[0]);
  k=-deflate(CM_CASE); if (p>2) k=kronll((llint) k,p); else k=A[k&7];
  if ((m&3)==2)
  {QD_powi(wmax,P,m/2-1,T); if (k==1) QD_neg(wmax,T,T);
   if (k==0) poly.deg=0; else QD_copy(wmax,T,poly.coeff[1]);
   QDpoly_mul(ef2,poly,&R,-1); delQDpoly(&ef2); delQDpoly(&poly); ef2=R;}
  if ((m&3)==0)
  {QD_powi(wmax,P,m/2,T); if (k==1) QD_neg(wmax,T,T);
   if (k==0) poly.deg=0; else QD_copy(wmax,T,poly.coeff[1]);
   QDpoly_mul(ef1,poly,&R,-1); delQDpoly(&ef1); delQDpoly(&poly); ef1=R;}
  initQDpoly(&poly,1); QD_copy(wmax,QD_one,poly.coeff[0]);
  if ((m&3)==0)
  {QD_powi(wmax,P,m/2-1,T); QD_neg(wmax,T,poly.coeff[1]);
   QDpoly_mul(ef2,poly,&R,-1); delQDpoly(&ef2); delQDpoly(&poly); ef2=R;}
  if ((m&3)==2)
  {if (m<=2) delQDpoly(&poly);
   else
   {QD_powi(wmax,P,m/2,T); QD_neg(wmax,T,poly.coeff[1]);
    QDpoly_mul(ef1,poly,&R,-1); delQDpoly(&ef1); delQDpoly(&poly); ef1=R;}}}
 QD_copy(wmax,P,T);
 for (i=1;i<=ef2.deg;i++)
 {QD_mul(wmax,ef2.coeff[i],T,ef2.coeff[i]); QD_mul(wmax,T,P,T);}
 QDpoly_inv(ef2,ef1.deg-ef2.deg,&R); QDpoly_mul(ef1,R,v,ef1.deg-ef2.deg);
 delQDpoly(&R); delQDpoly(&ef1); delQDpoly(&ef2);}

void euler_factor(llint p,int m,QDpoly *v)
{int i;
 if ((COND0%p)!=0) {euler_factor_good(p,m,v);}
 else
 {for (i=0;p!=badprimes[i];i++);
  if (HECKE) euler_factor_hecke_bad(p,badprimetype[i],m,v);
  else euler_factor_bad(p,badprimetype[i],m,v);}
 QDpoly_intround(v);}

static void localinfo(llint p,int sp)
{QDpoly v; printf("Euler factor at %lli is ",p);
 euler_factor(p,sp,&v); QDpoly_intout(v); delQDpoly(&v);
 if ((!HECKE) && (p<=3))
 {printf("Tame cond exponent at %lli is %i\n",p,get_tame_conductor(p,sp));
  printf("Wild cond exponent at %lli is %i\n",p,get_wild_conductor(p,sp));}
 else printf("Conductor exponent at %lli is %i\n",p,get_conductor(p,sp));}

static int primetest(llint x)
{int i; llint p; if (x==1) return FALSE;
 for(i=0;PRIMES[i]!=0;i++)
 {p=PRIMES[i]; if (x==p) return 1; if ((x%p)==0) return 0;} return 1;}

void localinfos(char *p,char *sp)
{int k=0; llint prl=0,prh=0,pr; int ml=0,mh=0,m;
 while ((p[k]!='-') && (p[k]!=0))
 {prl*=10; ASSERT(ISA_NUMBER(p[k])); prl+=p[k]-'0'; k++;}
 if (p[k]=='-') k++; else prh=prl;
 while (p[k]!=0)
 {prh*=10; ASSERT(ISA_NUMBER(p[k])); prh+=p[k]-'0'; k++;}
 k=0; while ((sp[k]!='-') && (sp[k]!=0))
 {ml*=10; ASSERT(ISA_NUMBER(sp[k])); ml+=sp[k]-'0'; k++;}
 if (sp[k]=='-') k++; else mh=ml;
 while (sp[k]!=0)
 {mh*=10; ASSERT(ISA_NUMBER(sp[k])); mh+=sp[k]-'0'; k++;}
 if ((ml<=0) || (mh>99)) errorit("Symmetric power range invalid");
 for (m=ml;m<=mh;m++)
 {if (HECKE) printf("Hecke "); printf("Symmetric power %i\n",m);
  for (pr=prl;pr<=prh;pr++) if (primetest(pr)) localinfo(pr,m);}}
