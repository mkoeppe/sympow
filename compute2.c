#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

void special_value(int p,QD s,int sp)
{QD t,u,rp,ip,C; int k,rpp=0,ipp=0;
 QD_copy(wmax,REAL_PERIOD,rp); QD_copy(wmax,IMAG_PERIOD,ip);
 QD_copy(wmax,COND[sp],C);
 if (VERBOSE>=2)
 {printf("rp: "); QD_output(wmax,16*wmax,rp);
  printf("ip: "); QD_output(wmax,16*wmax,ip);}
 if (HECKE)
 {switch (sp&3) /* L*(2*Pi)^[(sp-1)/2]/(ip)^(sp) for 3mod4 and ip->rp for 1*/
  {case 0: QD_powi(p,rp,sp+2,t); QD_powi(p,QD_twopi,sp/2-1,u);
           QD_div(p,u,t,t); QD_mul(p,rp,t,t); QD_mul(p,ip,t,t);
	   QD_powi(p,COND[sp],1+sp/4,u); QD_mul(p,t,u,t); QD_mul(p,s,t,s);
	   return; /* L*(2*Pi)^(sp/2-1)/rp^(sp+2) * (rp*ip)*N^(1+sp/4) */
   case 1: QD_powi(p,rp,sp,t); QD_powi(p,QD_twopi,(sp-1)/2,u);
           QD_div(p,u,t,t); QD_mul(p,s,t,s); return;
   case 2: QD_powi(p,ip,sp,t); QD_powi(p,QD_twopi,sp/2-1,u);
           QD_div(p,u,t,t); QD_mul(p,s,t,s); return;
   case 3: QD_powi(p,ip,sp,t); QD_powi(p,QD_twopi,(sp-1)/2,u);
           QD_div(p,u,t,t); QD_mul(p,s,t,s); return;}}
 switch(sp)
 {case 2: QD_mul(p,rp,ip,t); QD_mul(p,t,QD_twopi,t);
          QD_div(p,s,t,s); QD_mul1(p,s,(double) COND0,s); return;
  case 3: QD_powi(p,ip,3,t); QD_mul(p,t,rp,t);
          QD_div(p,s,t,s); QD_mul1(p,s,(double) COND0,s);
	  QD_mul(p,s,QD_twopi,s); return;
  case 5: QD_sqr(p,rp,t); QD_mul(p,t,ip,t);
          QD_mul1(p,QD_twopi,(double) COND0,t);
          QD_sqr(p,rp,u); QD_mul(p,u,ip,u);
          QD_div(p,t,u,t); QD_powi(p,t,3,t); QD_mul(p,s,t,s); return;
  case 6: QD_mul(p,s,QD_twopi,s); QD_mul(p,s,QD_twopi,s);
          QD_mul(p,rp,ip,t); QD_div1(p,t,(double) COND0,t);
          QD_powi(p,t,6,t); QD_div(p,s,t,s); return;
  case 7: QD_mul1(p,QD_twopi,(double) COND0,t); QD_div(p,t,rp,t);
          QD_powi(p,t,6,t); QD_powi(p,ip,10,u);
	  QD_div(p,t,u,t); QD_mul(p,s,t,s); return;
  case 9: QD_sqr(p,ip,t); QD_sqr(p,rp,u); QD_mul(p,u,rp,u); QD_mul(p,t,u,t);
          QD_mul1(p,QD_twopi,(double) COND0,u); QD_sqr(p,u,u);
	  QD_div(p,u,t,t); QD_powi(p,t,5,t); QD_mul(p,s,t,s); return;
  default:
   switch(sp&3)
   {case 0: return;
    case 1: k=(sp+1)/2; rpp=(k*(k+1))/2; ipp=(k*(k-1))/2;
            QD_mul1(p,QD_twopi,(double) COND0,u);
	    QD_div(p,u,ip,u); QD_powi(p,u,ipp,u); QD_powi(p,rp,rpp,t);
	    QD_div(p,u,t,t); QD_mul(p,s,t,s); return;
    case 2: k=((sp/2)*(sp/2+1))/2; QD_mul(p,rp,ip,t);
            QD_div1(p,t,(double) COND0,t); QD_powi(p,t,k,t);
            QD_powi(p,QD_twopi,3*k-(sp/2+1)*(sp/2+1),u);
	    QD_div(p,u,t,t); QD_mul(p,s,t,s); return;
    case 3: k=(sp+1)/2; rpp=(k*(k-1))/2; ipp=(k*(k+1))/2;
            QD_mul1(p,QD_twopi,(double) COND0,u);
	    QD_div(p,u,rp,u); QD_powi(p,u,rpp,u); QD_powi(p,ip,ipp,t);
	    QD_div(p,u,t,t); QD_mul(p,s,t,s); return;}}}

static int sbinom(int a,int b)
{int i; llint s; if (b&1) s=-1; else s=1; if (a==b) return (int) s;
 if (b>a/2) b=a-b; for (i=0;i<b;i++) s*=(a-i);
 for (i=2;i<=b;i++) s/=i; return (int) s;}

void pth_coeff_from_ap(int W,llint p,int ap,int sp,QD x) /* p<2^50 */
{QD t={0.0,0.0,0.0,0.0},u={0.0,0.0,0.0,0.0},v={0.0,0.0,0.0,0.0}; int i;
 QD y={0.0,0.0,0.0,0.0},z={0.0,0.0,0.0,0.0}; int W2=2,W3=3,W4=4; QD ppow[16];
 if (W<4) {W4=W; if (W<3) {W3=W; if (W<2) W2=W;}}
 QD_copy(wmax,QD_zero,x); if ((HECKE) && (ap==0)) return;
 if (HECKE) return power_trace_from_ap(W,p,ap,sp,x);
 if (ap==0) {if (!(sp&1)) {t[0]=(double) (-p); QD_powi(W,t,sp/2,x);} return;}
 switch(sp)
 {case 1: x[0]=(double) ap; return;
  case 2: x[0]=(double) ap; x[0]=x[0]*x[0]-(double) p; return; /* a^2-p */
  case 3: y[0]=(double) ap; x[0]=y[0]*y[0]-(double) (2*p);
          QD_mul1(W2,x,y[0],x); return; /* a^3-2pa */
  case 4: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
          QD_mul1(W2,y,y[0]-3.0*z[0],x); QD_sqr(W2,z,z);
	  QD_add(W2,x,z,x); return; /* a^4-3*a^2*p+p^2 */
  case 5: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
          QD_mul1(W2,z,3.0*z[0]-4.0*y[0],x);
	  QD_sqr(W2,y,t); QD_add(W2,t,x,x); QD_mul1(W3,x,(double) ap,x);
	  return; /* a^5-4*a^3*p+3*a*p^2 */
  case 6: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
          QD_sqr(W2,y,t); QD_mul1(W3,t,y[0]-5.0*z[0],t);
	  QD_sqr(W2,z,x); QD_mul1(W3,x,6.0*y[0]-z[0],x); QD_add(W3,x,t,x);
	  return; /* a^6-5*a^4*p+6*a^2*p^2-p^3 */
  case 7: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
          QD_sqr(W2,y,t); QD_mul1(W3,t,y[0]-6.0*z[0],t);
	  QD_sqr(W2,z,x); QD_mul1(W3,x,10.0*y[0]-4.0*z[0],x); QD_add(W3,x,t,x);
	  QD_mul1(W4,x,(double) ap,x); return;
	  /* a^7-6*a^5*p+10*a^3*p^2-4*a*p^3 [(7,0)+(6,1)+(5,2)+(4,3)]*/
  case 8: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
          QD_sqr(W2,y,t); QD_mul1(W3,t,-7*y[0]+15.0*z[0],u);
	  QD_sqr(W2,z,x); QD_mul1(W3,x,-10.0*y[0]+z[0],x);
	  QD_add(W3,x,u,x); QD_mul1(W4,x,z[0],x); QD_sqr(W4,t,t);
	  QD_add(W4,x,t,x); return; /* a^8-7*a^6*p+15*a^4*p^2-10*a^2*p^3+p^4 */
  case 9: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
          QD_sqr(W2,y,t); QD_mul1(W3,t,-8*y[0]+21.0*z[0],u);
	  QD_sqr(W2,z,x); QD_mul1(W3,x,-20.0*y[0]+5.0*z[0],x);
	  QD_add(W3,x,u,x); QD_mul1(W4,x,z[0],x); QD_sqr(W4,t,t);
	  QD_add(W4,x,t,x); QD_mul1(W4,x,(double) ap,x); return;
	  /* a^9-8*a^7*p+21*a^5*p^2-20*a^3*p^3+5*a*p^4 [1,-8,21,-20,5] */
  case 10: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
           QD_sqr(W2,y,t); QD_mul1(W3,t,y[0]-9*z[0],u);
	   QD_sqr(W2,z,x); QD_mul1(W3,x,28*y[0]-35*z[0],v);
	   QD_add(W4,v,u,v); QD_mul(W4,t,v,v);
	   QD_sqr(W4,x,x); QD_mul1(W4,x,15*y[0]-z[0],u); QD_add(W4,u,v,x);
	   return; /* [1,-9,28,-35,15,-1] */
  case 11: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
           QD_sqr(W2,y,t); QD_mul1(W3,t,y[0]-10*z[0],u);
	   QD_sqr(W2,z,x); QD_mul1(W3,x,36*y[0]-56*z[0],v);
	   QD_add(W4,v,u,v); QD_mul(W4,t,v,v); /* [1,-10,36,-56,35,-6] */
	   QD_sqr(W4,x,x); QD_mul1(W4,x,35*y[0]-6*z[0],u);
	   QD_add(W4,u,v,x); QD_mul1(W4,x,(double) ap,x); return;
  case 12: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
           QD_sqr(W2,y,t); QD_mul1(W3,t,y[0]-11*z[0],u);
	   QD_sqr(W2,z,x); QD_mul1(W3,x,45*y[0]-84*z[0],v);
	   QD_add(W4,v,u,v); QD_mul(W4,t,v,v);
	   QD_sqr(W4,x,t); QD_mul1(W4,t,70*y[0]-21*z[0],u);
	   QD_add(W4,u,v,u); QD_mul(W4,u,y,u); /* [1,-11,45,-84,70,-21,1] */
	   QD_mul(W4,t,x,y); QD_add(W4,u,y,x); return;
  case 13: y[0]=(double) ap; y[0]=y[0]*y[0]; z[0]=(double) p;
           QD_sqr(W2,y,t); QD_mul1(W3,t,y[0]-12*z[0],u);
	   QD_sqr(W2,z,x); QD_mul1(W3,x,55*y[0]-120*z[0],v);
	   QD_add(W4,v,u,v); QD_mul(W4,t,v,v);
	   QD_sqr(W4,x,t); QD_mul1(W4,t,126*y[0]-56*z[0],u);
	   QD_add(W4,u,v,u); QD_mul(W4,u,y,u);
	   QD_mul(W4,t,x,y); QD_mul1(W4,y,7.0,y); QD_add(W4,u,y,t);
	   QD_mul1(W4,t,(double) ap,x); return; /* [1,-12,55,-120,126,-56,7] */
  default: y[0]=(double) ap; y[0]=y[0]*y[0];
	   QD_copy(W4,QD_zero,ppow[1]); ppow[1][0]=(double) p;
	   for (i=2;i<=sp/2;i++) QD_mul(W4,ppow[1],ppow[i-1],ppow[i]);
	   QD_copy(W4,y,u); QD_copy(W4,ppow[sp/2],z);
	   QD_mul1(W4,z,(double) sbinom(sp-sp/2,sp/2),z);
	   for (i=1;i<sp/2;i++)
	   {QD_mul(W4,u,ppow[sp/2-i],t); QD_mul(W4,u,y,u);
 	    QD_mul1(W4,t,(double) sbinom(sp-sp/2+i,sp/2-i),v);
	    QD_add(W4,v,z,z);}
	   QD_add(W4,z,u,x); if (sp&1) QD_mul1(W4,x,(double) ap,x);
	   return;}}

static int binom(int a,int b)
{int i; llint s=1; if ((a<b) || (b<0)) return 0; if (a==b) return (int) s;
 if (b>a/2) b=a-b; for (i=0;i<b;i++) s*=(a-i);
 for (i=2;i<=b;i++) s/=i; return (int) s;}

static int tbinom(int a,int b)
{int s=binom(a-b+1,b)-binom(a-b-1,b-2); if (b&1) s=-s; return s;}

void power_trace_from_ap(int W,llint p,int ap,int sp,QD x) /* p<2^50 */
{QD t={0.0,0.0,0.0,0.0},u={0.0,0.0,0.0,0.0},v={0.0,0.0,0.0,0.0}; int i;
 QD y={0.0,0.0,0.0,0.0},z={0.0,0.0,0.0,0.0}; int W2=2,W3=3,W4=4; QD ppow[16];
 if (DEBUG) printf("power_trace_from_ap %lli %i %i\n",p,ap,sp);
 if (W<4) {W4=W; if (W<3) {W3=W; if (W<2) W2=W;}}
 QD_copy(wmax,QD_zero,x); if (sp==1) {x[0]=(double) ap; return;}
 y[0]=(double) ap; y[0]=y[0]*y[0];
 QD_copy(W4,QD_one,ppow[0]); QD_copy(W4,QD_zero,ppow[1]);
 ppow[1][0]=(double) p;
 for (i=2;i<=sp/2;i++) QD_mul(W4,ppow[1],ppow[i-1],ppow[i]);
 QD_copy(W4,y,u); QD_copy(W4,ppow[sp/2],z);
 QD_mul1(W4,z,(double) tbinom(sp,sp/2),z);
 for (i=1;i<sp/2;i++)
 {QD_mul(W4,u,ppow[sp/2-i],t); QD_mul(W4,u,y,u);
  QD_mul1(W4,t,(double) tbinom(sp,sp/2-i),v); QD_add(W4,v,z,z);}
 QD_add(W4,z,u,x); if (sp&1) QD_mul1(W4,x,(double) ap,x);}

