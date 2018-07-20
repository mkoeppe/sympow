#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static void real_root(QD r1)
{QD U,R,A3,sU,T,cS,cT,M,isU,O,V,S;
 if (DEBUG) printf("real root\n");
 QD_div1(wmax,Edisc,-1728.0,U); QD_div1(wmax,Ec6,1728.0,R);
 QD_div1(wmax,Eb2,12.0,A3);
 if (U[0]>0.0) /* negative discriminant */
 {QD_sqrt(wmax,U,sU); QD_add(wmax,R,sU,S); QD_sub(wmax,R,sU,T);
  QD_cbrt(wmax,S,cS); QD_cbrt(wmax,T,cT);
  QD_add(wmax,cS,cT,T); QD_sub(wmax,T,A3,r1); return;}
 QD_neg(wmax,U,U); QD_sqr(wmax,R,T);
 QD_add(wmax,T,U,T); QD_sqrt(wmax,T,M);
 if (R[0]==0.0) {QD_copy(wmax,QD_pi,O); QD_mul2n(wmax,O,-1,O);}
 else
 {QD_sqrt(wmax,U,isU); QD_div(wmax,isU,R,T); QD_atan(wmax,T,O);
  if (O[0]<0.0) QD_add(wmax,O,QD_pi,O);}
 QD_div1(wmax,O,3.0,O); QD_cos(wmax,O,T); QD_cbrt(wmax,M,U);
 QD_mul(wmax,T,U,V); QD_mul2n(wmax,V,1,T); QD_sub(wmax,T,A3,r1);}

static void more_roots(QD r1,QD r2,QD r3)
{QD A,b,Z,T,c,U,D,F,S;
 if (DEBUG) printf("more roots\n");
 QD_mul2n(wmax,Eb2,-2,A); QD_add(wmax,A,r1,b); QD_mul2n(wmax,Eb4,-1,Z);
 QD_mul(wmax,b,r1,T); QD_add(wmax,T,Z,c); QD_sqr(wmax,b,T);
 QD_mul2n(wmax,c,2,U); QD_sub(wmax,T,U,D); QD_sqrt(wmax,D,S); QD_neg(wmax,b,F);
 QD_add(wmax,F,S,r2); QD_mul2n(wmax,r2,-1,r2);
 QD_sub(wmax,F,S,r3); QD_mul2n(wmax,r3,-1,r3);}

void do_periods()
{QD a,r1,r2,r3,r12,r13,r23,p1,p2,t; QD T,A,B,NB,P1,P2,E1,E2,E3;
 if (Edisc[0]>0.0)
 {real_root(r1); more_roots(r1,r2,r3);
  if (r1[0]<r2[0])
  {QD_copy(wmax,r1,t); QD_copy(wmax,r1,r2); QD_copy(wmax,r2,t);}
  if (r1[0]<r3[0])
  {QD_copy(wmax,r1,t); QD_copy(wmax,r1,r3); QD_copy(wmax,r3,t);}
  if (r2[0]<r3[0])
  {QD_copy(wmax,r2,t); QD_copy(wmax,r2,r3); QD_copy(wmax,r3,t);}
  QD_sub(wmax,r1,r2,a); QD_sqrt(wmax,a,r12);
  QD_sub(wmax,r1,r3,a); QD_sqrt(wmax,a,r13);
  QD_sub(wmax,r2,r3,a); QD_sqrt(wmax,a,r23);
  QD_agm(wmax,r12,r13,p1); QD_agm(wmax,r23,r13,p2);
  QD_div(wmax,QD_pi,p1,REAL_PERIOD); QD_div(wmax,QD_pi,p2,IMAG_PERIOD);}
 else
 {real_root(r1); 
  QD_mul2n(wmax,Eb2,-2,T); QD_mul1(wmax,r1,3.0,A); QD_add(wmax,A,T,A);
   /* A=3*r+b2/4 */  /* B=sqrt(3*r^2+r*b2/2+b4/2) */
  QD_sqr(wmax,r1,T); QD_mul1(wmax,T,3,B); QD_mul2n(wmax,Eb2,-1,T);
  QD_mul(wmax,r1,T,T); QD_add(wmax,B,T,B); QD_mul2n(wmax,Eb4,-1,T);
  QD_add(wmax,T,B,B); QD_sqrt(wmax,B,B);
   /* agm(2*sqrt(B),sqrt(2*B+A)) and agm(2*sqrt(B),sqrt(2*B-A)) */
  QD_sqrt(wmax,B,T); QD_mul2n(wmax,T,1,E1); QD_mul2n(wmax,B,1,NB);
  QD_add(wmax,A,NB,T); QD_sqrt(wmax,T,E2);
  QD_sub(wmax,NB,A,T); QD_sqrt(wmax,T,E3);
  QD_agm(wmax,E1,E2,P1); QD_agm(wmax,E1,E3,P2);
  QD_div(wmax,QD_twopi,P1,REAL_PERIOD); QD_div(wmax,QD_pi,P2,IMAG_PERIOD);}
 if (DEBUG)
 {printf("real: "); QD_output(wmax,16*wmax,REAL_PERIOD);
  printf("imag: "); QD_output(wmax,16*wmax,IMAG_PERIOD);}}
