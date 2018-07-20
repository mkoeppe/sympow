#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static void compute_local_powser(QD h,QD *P,QD R,int NUM,int W)
{int i; if (W==1) for (i=NUM;i>=0;i--) {R[0]*=h[0]; R[0]+=P[i][0];}
 else for (i=NUM;i>=0;i--) {QD_mul(W,R,h,R); QD_add(W,R,P[i],R);}}

static void loopit(x,left_point,step,index,which,result,W,S)
double *x,left_point,step; int index,which; double *result; int W,S;
{double Y; QD z,*P; int n;
 QD_copy(W,x,z); z[0]-=left_point;
 Y=Floor(0.5+z[0]/step); z[0]-=step*Y; n=(int) Y;
 QD_self_renorm(W,z); P=TABLE[which][index+n];
 QD_copy(W,QD_zero,result); compute_local_powser(z,P,result,3+wprec[S]/2,W);}

void get_wt_large(int which,QD x,QD R,int W,int S)
{int i=0; double y,tb=TOO_BIG[which],B=EXPAND0_LIM[which];
 double STEP=STEP_SIZE[which],LP=STEP*32.0;
 y=x[0]; if (y>tb) {QD_copy(W,QD_zero,R); return;}
 while (1)
 {if (y<B+B) {loopit(x,LP,STEP,STEPS*i,which,R,W,S); break;}
  i++; B+=B; LP+=LP; STEP+=STEP;}}

static void ps_tackon(int W,QD R,QD L,int which)
{QD M,T; int i;
 QD_powi(W,L,NUM_LOGS[which],M);
 for (i=0;i<TACKON[which];i++)
 {QD_mul(W,M,L,M); QD_mul(W,TACKS[which][i],M,T); QD_add(W,T,R,R);}}

static void get_wt_powser(int which,QD x,QD L,QD R,int W)
{int i,l; QD T,xs,**P=POWSER[which],*Q;
 QD_copy(W,QD_zero,R); QD_sqr(W,x,xs);
 for(l=NUM_LOGS[which];l>=0;l--)
 {Q=P[l]; QD_mul(W,R,L,R); QD_copy(W,QD_zero,T);
  if ((l==NUM_LOGS[which]) && HALF_ZERO[which])
  {for(i=A0PT/2-1;i>=0;i--) {QD_mul(W,T,xs,T); QD_add(W,T,Q[i],T);}
   if (HALF_ZERO[which]==2) QD_mul(W,T,x,T);
  }
  else {for(i=A0PT-1;i>=0;i--) {QD_mul(W,T,x,T); QD_add(W,T,Q[i],T);}}
  QD_add(W,R,T,R);}
 if (TACKON[which]) ps_tackon(W,R,L,which);}

#define FCHECK FALSE

void get_weight(llint n,QD WEIGHT,int S)
{QD nQD,nQDi,x1,x2,L1,L2,r1,r2,t;
 int wh1=whi[S],wh2=wlo[S],sp=SYMPOW[wh1],W=w[S];
 QD_copy(W,QD_zero,nQD); nQD[0]=(double) n; QD_div(W,QD_one,nQD,nQDi);
 QD_mul(W,nQD,DECAY[sp],x1); if (WIGGLE[S][0]!=1.0) QD_mul(W,x1,WIGGLE[S],x1);
 if (x1[0]<=EXPAND0_LIM[wh1]) {QD_log(W,x1,L1); get_wt_powser(wh1,x1,L1,r1,W);}
 else get_wt_large(wh1,x1,r1,W,S);
#if 0
 if (n<10) {printf("wt %d ",n); QD_output(W,16*W,r1);}
 if (n<10) {printf("in %d ",n); QD_output(W,16*W,x1);}
#endif
 if ((FCHECK) && (r1[0]!=0.0))
 {printf("v="); QD_output(W,16*W,x1); printf("r="); QD_output(W,16*W,r1);
  printf("u=ev(PH,v,600)-r; if (abs(u)>1.0e-50,[v,r,u]);\n");}
 QD_powi(W,nQDi,evalpt[wh1],t); QD_mul(W,r1,t,r1);
 if ((HECKE) && (sp&1)==0) QD_mul1(W,r1,(double) (sp/2),r1);
 if (DEBUG>=2) printf("MID %f %f\n",r1[0],WIGSQI[S][0]);
 if ((sp&1) && (WIGGLE[S][0]==1.0)) {QD_mulall(W,r1,2.0,WEIGHT); return;}
 if (WIGGLE[S][0]!=1.0) QD_mul(W,x1,WIGSQI[S],x2); else QD_copy(W,x1,x2);
 if (x2[0]<=EXPAND0_LIM[wh2]) {QD_log(W,x2,L2); get_wt_powser(wh2,x2,L2,r2,W);}
 else get_wt_large(wh2,x2,r2,W,S);
 if ((FCHECK) && (r2[0]!=0.0))
 {printf("v="); QD_output(W,16*W,x2); printf("r="); QD_output(W,16*W,r2);
  printf("u=ev(PL,v,600)-r; if (abs(u)>1.0e-50,[v,r,u]);\n");}
 if ((sp&1)==0) QD_mul(W,r2,DECAY[sp],r2);
 QD_powi(W,nQDi,evalpt[wh2],t); QD_mul(W,r2,t,r2); QD_add(W,r1,r2,WEIGHT);
}
