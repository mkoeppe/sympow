#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

#define tlc tame_local_conductor
#define wlc wild_local_conductor

int tame_local_conductor(int type,int sympow)
{if (sympow<=0) return 0;
 switch(type)
 {case 0: errorit("case 0 in tame_local_conductor"); /* good */
  case 1: return sympow; /* multiplicative */
  case 29: case 30: case 31: /* potentially multiplicative */
   if (sympow&1) return sympow+1; else return sympow;
  case 2: case 16: case 17: /* C2 */
   if (sympow&1) return sympow+1; else return 0;
  case 3: case 14: case 18:  /* C3 */
  switch(sympow%3)
  {case 0: return (sympow*2)/3; case 1: return sympow+1-(sympow-1)/3;
   case 2: return sympow+1-(sympow+1)/3;}
  case 4: case 22: case 23: /* C4 */
  switch(sympow&3)
  {case 1: case 3: return sympow+1;
   case 0: return sympow/2; case 2: return 1+sympow/2;}
  case 6: case 15: case 19: case 20: case 21: /* C6 */
  switch(sympow%6)
  {case 1: case 3: case 5: return sympow+1;
   case 0: return (sympow*2)/3; case 4: return sympow+1-(sympow-1)/3;
   case 2: return sympow+1-(sympow+1)/3;}
  case 8: case 9: case 10: /* Q8 */
  switch(sympow&3)
  {case 1: case 3: return (sympow+1);
   case 0: return (3*sympow)/4; case 2: return (3*sympow+6)/4;}
  case 12: case 13: /* C3x|C4 */
  switch(sympow%12)
  {case 1: case 3: case 5: case 7: case 9: case 11: return sympow+1;
   case 0: return (5*sympow)/6; case 2: return (5*sympow+8)/6;
   case 4: return (5*sympow+4)/6; case 6: return (5*sympow)/6+1;
   case 8: return (5*sympow+2)/6; case 10: return (5*sympow+10)/6;}
  case 24: case 25: case 26: case 27: /* SL2F3 */
  switch(sympow%12)
  {case 1: case 3: case 5: case 7: case 9: case 11: return sympow+1;
   case 0: return (11*sympow)/12; case 2: return (11*sympow+14)/12;
   case 4: return (11*sympow+16)/12; case 6: return (11*sympow+6)/12;
   case 8: return (11*sympow+8)/12; case 10: return (11*sympow+22)/12;}
  default: printf("%i ",type); errorit("Bad type in tame_local_conductor");}
 return -1;}

int wild_local_conductor(int type,int m)
{if (m<=0) return 0;
 switch(type)
 {case 0: case 1: case 2: case 3: case 4: case 6: case 31: return 0; /* p>5 */
  case 8: return tlc(8,m)+tlc(2,m)/2;
  case 9: return tlc(8,m)+tlc(2,m);
  case 10: return tlc(8,m)+tlc(4,m)+tlc(2,m);
  case 12: return tlc(3,m)/2;
  case 13: return 3*tlc(3,m)/2;
  case 14: case 15: case 18: case 19: return tlc(3,m);
  case 16: case 20: return tlc(2,m);
  case 17: case 21: return 2*tlc(2,m);
  case 22: case 23: return 2*tlc(4,m)+tlc(2,m);
  case 24: return (2*tlc(8,m)+tlc(2,m))/6;
  case 25: return (tlc(8,m)+2*tlc(2,m))/3;
  case 26: return (tlc(8,m)+5*tlc(2,m))/3;
  case 27: return (10*tlc(8,m)+5*tlc(2,m))/6;
  case 29: return 2*tlc(2,m);
  case 30: return tlc(2,m);
  default: printf("%i ",type); errorit("Bad type in wild_local_conductor");}
 return -1;}

void badprimetype_output(int x,llint p)
{switch(x)
 {case 0: printf("C1 ADDITIVE REDUCTION\n"); break;
  case 1: printf("C1 MULTIPLICATIVE REDUCTION\n"); break;
  case 2: case 3: case 4: case 6:
  {if ((p%x)==1) printf("C%i abelian\n",x);
   else printf("C%i nonabelian\n",x); break;}
  case 8: printf("Q8 v2(N)=5\n"); break;
  case 9: printf("Q8 v2(N)=6\n"); break;
  case 10: printf("Q8 v2(N)=8\n"); break;
  case 12: printf("C3x|C4 v3(N)=3\n"); break;
  case 13: printf("C3x|C4 v3(N)=5\n"); break;
  case 14: printf("C3 p=3 abelian\n"); break;
  case 15: printf("C6 p=3 abelian\n"); break;
  case 16: printf("C2 p=2 v2(N)=4\n"); break;
  case 17: printf("C2 p=2 v2(N)=6\n"); break;
  case 18: printf("C3 p=3 nonabelian\n"); break;
  case 19: printf("C6 p=3 nonabelian\n"); break;
  case 20: printf("C6 p=2 v2(N)=4\n"); break;
  case 21: printf("C6 p=2 v2(N)=6\n"); break;
  case 22: printf("C4 p=2 abelian\n"); break;
  case 23: printf("C4 p=2 nonabelian\n"); break;
  case 24: printf("SL2F3 v2(N)=3\n"); break;
  case 25: printf("SL2F3 v2(N)=4\n"); break;
  case 26: printf("SL2F3 v2(N)=6\n"); break;
  case 27: printf("SL2F3 v2(N)=7\n"); break;
  case 29: printf("C2 p=2 MULTIPLICATIVE v2(N)=6\n"); break;
  case 30: printf("C2 p=2 MULTIPLICATIVE v2(N)=4\n"); break;
  case 31: printf("C2 MULTIPLICATIVE\n"); break;
  default: errorit("Unknown badprime type");}}

static int bpt_at(double p)/* p>3.0 */
{int v4,v6,vd; if (p<=3.0) errorit("prime too small in bpt_at");
 if (!QD_is_divisible(Edisc,p)) return 0;
 if (!QD_is_divisible(Ec4,p)) return 1;
 v4=QD_valuation(Ec4,p); v6=QD_valuation(Ec6,p); vd=QD_valuation(Edisc,p);
 if ((v4==2) && (v6==3) && (vd>6)) return 31;
 if (vd==6) return 2; if ((vd==2) || (vd==10)) return 6;
 if ((vd==3) || (vd==9)) return 4; if ((vd==4) || (vd==8)) return 3;
 errorit("bad valuation in bpt_ap"); return -1;}

static int ecfp3(int c4t,int c6t,int Dt)
{int c4r,c6r; if (DEBUG) printf("ecfp3 %i %i %i\n",c4t,c6t,Dt);
 if ((Dt%3)!=0) return 0; if ((c6t%27)!=0) return 1;
 c4r=(c4t/9)%3; c6r=(c6t/27)%9;
 if (c4r==0)
 {if (c6r==0) {if ((c4t%81)==0) return 5; return 4;}
  switch (c6r)
  {case 1: case 2: case 7: case 8: return 3;
   case 3: case 6: return 5; case 4: case 5: return 2;}}
 if (c4r==1)
 {switch(c6r)
  {case 0: return 2; case 1: case 3: case 6: case 8: return 3;
   case 2: case 4: case 5: case 7: return 4;}}
 if (c4r==2)
 {switch(c6r)
  {case 0: return 2; case 2: case 7: return 2;
  case 1: case 3: case 4: case 5: case 6: case 8: return 3;}}
 return 0;
}

static int bpt_at3()
{int v4,v6,vd,t,u,c4,c6,D,c4t,c6t,Dt;
 if (DEBUG) printf("bpt_at3\n");
 if (!QD_is_divisible(Edisc,3.0)) return 0;
 if (!QD_is_divisible(Ec4,3.0)) {fp3=1; return 1;}
 v4=QD_valuation(Ec4,3.0); v6=QD_valuation(Ec6,3.0);
 vd=QD_valuation(Edisc,3.0); c4=QD_modi(Ec4,387420489);
 c6=QD_modi(Ec6,387420489); D=QD_modi(Edisc,387420489);
 if (DEBUG) printf("%i %i %i\n",c4,c6,D);
 if ((v4==2) && (v6==3) && (v6!=5) && (vd>6)) {fp3=2; return 31;}
 if ((v4>=2) && ((v6>=3) && (v6!=5)) && (vd==6)) {fp3=2; return 2;}
 if ((v4<2) || (v6<3) || (v6==5) || (vd<6)) {c4t=c4; c6t=c6; Dt=D;}
 else {c4t=c4/9; c6t=(c6-387420489)/(-27); Dt=D/729;}
 fp3=ecfp3(c4t,c6t,Dt); ASSERT(fp3>=2);
 if (fp3==2) return 4; if (fp3==3) return 12; if (fp3==5) return 13;
 if (vd&3) t=15; else t=14;
 if (((c6t%243)==0) && ((c4t%81)==54)) t+=4;
 if ((c4t%27)==9) {u=c6t%243; if ((u==54) || (u==189)) t+=4;}
 return t;
}

static int ecfp2(int c4,int c6)
{int c4r,c6r; c4r=c4&15; c6r=c6&15;
 if ((c4r&3)==1)
 {switch(c6r&3)
  {case 0: return 6; case 1: return 4; case 3: return 3;
  case 2:
  {if (c4r==1) switch(c6r)
   {case 2: return 4; case 6: case 10: return 5; case 14: return 3;}
   if (c4r==5) switch(c6r)
   {case 2: return 3; case 6: return 2; case 10: case 14: return 4;}
   if (c4r==9) switch(c6r)
   {case 2: case 14: return 5; case 6: return 3; case 10: return 4;}
   if (c4r==13) switch(c6r)
   {case 2: case 6: return 4; case 10: return 3; case 14: return 2;}}}}
 if ((c4r&3)==2)
 {switch(c6r&3)
  {case 1: return 3; case 2: return 6; case 3: return 4;
   case 0: if (c6r&7) return 7; else return 8;}}
 if ((c4r&3)==3)
 {switch(c6r&3)
  {case 0: return 5; case 1: return 2; case 2: return 7; case 3: return 4;}}
 if (((c4r&3)==0) && ((c6r&3)!=0))
 {switch(c6r&3) {case 1: return 2; case 3: return 4;}}
 c4r=(c4/4)&1; c6r=(c6/4)&3;
 if (c6r==3) return 4; if (c4r) return 3; return 2;
}

static int bpt_at2()
{int v4,v6,vd,e,c4,c6,D,c4t,c6t,Dt,t1=FALSE;
 if (!QD_is_divisible(Edisc,2.0)) {fp2=0; return 0;}
 if (!QD_is_divisible(Ec4,2.0)) {fp2=1; return 1;}
 v6=QD_valuation(Ec6,2.0); if (v6==3) {fp2=0; return 0;}
 vd=QD_valuation(Edisc,2.0); e=minimum(v6/3,vd/6); v4=QD_valuation(Ec4,2.0);
 c4=QD_modi(Ec4,1<<27); c6=QD_modi(Ec6,1<<27); D=QD_modi(Edisc,1<<27);
 v4-=2*e; v6-=3*e; if ((v4==2) || (v4==3) || (v6==4)) e--;
 c4t=c4>>(2*e); c6t=c6>>(3*e); Dt=D>>(6*e);
 if ((c6t&31)!=0)
 {if ((Dt&1)==0) {if (e&1) {fp2=6; return 29;} fp2=4; return 30;}
  else {if (e&1) {fp2=6; return 17;} fp2=4; return 16;}}
 fp2=ecfp2(c4t/16,c6t/32); if (fp2==4) {fp2=ecfp2(c4t/16,-c6t/32); t1=TRUE;}
 ASSERT(fp2>=2); ASSERT(fp2!=4);
 if (fp2==2)
 {if ((e&1)==0) {if (t1) {fp2=4; return 20;} return 3;} fp2=6; return 21;}
 if (fp2==3)
 {if ((e&1)==0) {if (t1) {fp2=4; return 25;} return 24;} fp2=6; return 26;}
 if (fp2==5) {if ((e&1)==0) return 8; fp2=6; return 9;}
 if (fp2==6) {if ((e&1)==1) {fp2=5; return 8;} return 9;}
 if (fp2==7) return 27;
 if (fp2==8)
 {if ((c6t&511)==0) return 10; if ((c4t&127)==32) return 23; return 22;}
 return -1;
}

static void CMUL(llint p) /* multiply COND0 by p and check if overflow */
{if (((COND0*p)/p)!=COND0) errorit("Curve conductor is too large"); COND0*=p;}

int do_badprime(llint p) /* E is a minimal model */
{int f=0,i;
 if (p>3) {f=bpt_at((double) p); if (f==1) CMUL(p); else CMUL(p*p);}
 if (p==3) {f=bpt_at3(); for (i=0;i<fp3;i++) CMUL(3);}
 if (p==2) {f=bpt_at2(); for (i=0;i<fp2;i++) CMUL(2);} return f;}

static int lkup(int m)
{if ((m>=10) && (m<=20)) return 0; m=m%10; if (m>=4) return 0; return m;}

void compute_conductor(int m)
{int i,T,W; llint P; QD A; QDpoly v; char STR[4][4]={"th","st","nd","rd"};
 QD_copy(wmax,QD_one,COND[m]);
 for (i=0;badprimes[i]!=0;i++)
 {T=tlc(badprimetype[i],m); P=badprimes[i]; W=wlc(badprimetype[i],m);
  if (VERBOSE) printf("sp %i: Conductor at %lli is %i+%i, root number is %i\n",
		      m,P,T,W,local_rootno(i,badprimes[i],m));
  QD_copy(wmax,QD_zero,A); A[0]=(double) P; QD_powi(wmax,A,T+W,A);
  QD_mul(wmax,COND[m],A,COND[m]);
  if (VERBOSE)
  {printf("sp %i: Euler factor at %lli is ",m,P);
   euler_factor(P,m,&v); QDpoly_intout(v); delQDpoly(&v);}}
 if (VERBOSE)
 {printf("%i%s sym power conductor is ",m,STR[lkup(m)]);
  QD_intout_noplus(COND[m]);
  printf(", global root number is %i\n",global_rootno(m));}}

static int prim_part(int CM)
{if ((CM==-27) || (CM== -12)) return 3; if (CM==-28) return 7;
 if ((CM==-16) || (CM==-4) || (CM==-8)) return 2; return -CM;}

static int f(int CM)
{if ((CM==-16) || (CM==-4)) return 2; if (CM==-8) return 3; return 1;}

void compute_conductor_hecke(int m)
{int i,T,W,u; llint P=0; QD A; QDpoly v; char STR[4][4]={"th","st","nd","rd"};
 if (!CM_CASE) errorit("not CM case in compute_conductor_hecke\n");
 QD_copy(wmax,QD_one,COND[m]);
 for (i=0;badprimes[i]!=0;i++)
 {T=tlc(badprimetype[i],m)-tlc(badprimetype[i],m-2);
  W=wlc(badprimetype[i],m)-wlc(badprimetype[i],m-2); P=badprimes[i];
  if ((m&3)==0) {u=prim_part(CM_CASE); if (P==u) {T+=1; W+=f(CM_CASE)-1;}}
  if ((m&3)==2) {u=prim_part(CM_CASE); if (P==u) {T-=1; W-=f(CM_CASE)-1;}}
  if (VERBOSE) printf("sp %i: Conductor at %lli is %i+%i\n",m,P,T,W);
  QD_copy(wmax,QD_zero,A); A[0]=(double) P; QD_powi(wmax,A,T+W,A);
  QD_mul(wmax,COND[m],A,COND[m]);
  if (VERBOSE)
  {printf("sp %i: Euler factor at %lli is ",m,P);
   euler_factor(P,m,&v); QDpoly_intout(v); delQDpoly(&v);}}
 if (VERBOSE)
 {printf("%i%s sym power conductor is ",m,STR[lkup(m)]);
  QD_intout_noplus(COND[m]); printf("\n");}}

int get_tame_conductor(llint p,int m)
{int i; if ((COND0%p)!=0) return 0;
 if (HECKE) errorit("Tame conductor for Hecke?");
 for (i=0;badprimes[i]!=p;i++); return(tlc(badprimetype[i],m));}

int get_wild_conductor(llint p,int m)
{int i; if ((COND0%p)!=0) return 0;
 if (HECKE) errorit("Wild conductor for Hecke?");
 for (i=0;badprimes[i]!=p;i++); return(wlc(badprimetype[i],m));}

int get_conductor(llint p,int m)
{int T,W,i,u; if ((COND0%p)!=0) return 0;
 for (i=0;badprimes[i]!=p;i++);
 if (!HECKE) return(tlc(badprimetype[i],m)+wlc(badprimetype[i],m));
 T=tlc(badprimetype[i],m)-tlc(badprimetype[i],m-2);
 W=wlc(badprimetype[i],m)-wlc(badprimetype[i],m-2);
 if ((m&3)==0) {u=prim_part(CM_CASE); if (p==u) T+=f(CM_CASE);}
 if ((m&3)==2) {u=prim_part(CM_CASE); if (p==u) T-=f(CM_CASE);} return T+W;}
