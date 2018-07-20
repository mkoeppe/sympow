#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static int local_rootno_at2(int b)
{int c4,c6,D,v4,v6,vd,t,i; QD T; if (DEBUG) printf("local_rootno_at2\n");
 if (b==1) return -ec_ap(Ec4,Ec6,2);
 v4=QD_valuation(Ec4,2.0); v6=QD_valuation(Ec6,2.0);
 vd=QD_valuation(Edisc,2.0);
 if (v4>100) v4=100; if (v6>100) v6=100; if (vd>100) vd=100;
 QD_copy(wmax,Ec6,T); for (i=0;i<v6;i++) QD_mul1(wmax,T,0.5,T);
 c6=QD_modi(T,64.0); if (b>=29) return ((c6%4)==3) ? 1:-1;
 QD_copy(wmax,Ec4,T); for (i=0;i<v4;i++) QD_mul1(wmax,T,0.5,T);
 c4=QD_modi(T,64.0);
 QD_copy(wmax,Edisc,T); for (i=0;i<vd;i++) QD_mul1(wmax,T,0.5,T);
 D=QD_modi(T,64.0);
 if (DEBUG) printf("r2 c4:%i c6:%i D:%i v4:%i v6:%i vD:%i\n",c4,c6,D,v4,v6,vd);
 if ((v4==4) && (v6==5) && (vd==4)) /* cases 1,7,13 */
 {if ((c4%4)==(c6%4)) return ((c4%4)==1) ? 1:-1; /* case 1 */
  if ((c4%4)==1) return (((c4*c6)%8)==3) ? 1:-1; return -1;} /* case 7,13 */
 if ((v4==5) && (v6==5) && (vd==4)) /* cases 2,8 */
 {if ((c6%4)==3) return 1; if ((c6%8)==5) return 1; return -1;}
 if ((v4>=6) && (v6==5) && (vd==4)) return -1; /* cases 2,14 */
 if ((v4==4) && (v6>=7) && (vd==6)) /* case 3 and 9 it seems */
 {if ((COND0%64)==0) return (v6==7) ? 1:-1; /* case 3 */
  QD_mulall(wmax,Ec6,0.0078125,T); t=(c4-4*QD_modi(T,16))%16; if (t<0) t+=16;
  if ((t==7) || (t==11)) return 1; return -1;} /* case 9 */
 if ((v4==4) && (v6==6) && (vd==7)) /* case 4 */
 {if (((c6%8)==5) || ((c6%8)==((5*c4)%8))) return 1; return -1;}
 if ((v4==5) && (v6==6) && (vd==6)) return ((c4%4)==3) ? 1:-1; /* case 5 */
 if ((v4>=6) && (v6==6) && (vd==6)) return ((c6%4)==1) ? 1:-1; /* case 6 */
 if ((v4==5) && (v6==7) && (vd==8)) /* case 10 */
 {if (((c6%8)==3) || (((2*c4+c6)%8)==7)) return 1; return -1;}
 if ((v4==5) && (v6==8) && (vd==9)) /* case 11 */
 {if ((((2*c6+c4)%8)==1) || (((2*c6+c4)%8)==7)) return 1; return -1;}
 if ((v4==5) && (v6>=9) && (vd==9)) return (((c4%8)==1) || ((c4%8)==3)) ? 1:-1;
 if ((v4==4) && (v6==6) && (vd==8)) /* cases 15,25,30 */
 {t=(2*c6+c4)%16; if (t==11) return -1; /* case 15 */
  if (t==7) return (((2*c6+c4)%32)==23) ? 1:-1; /* case 30 */
  return (((2*c6+c4)%16)==3) ? 1:-1;} /* case 25 */
 if ((v4==6) && (v6==7) && (vd==8)) /* cases 27,31 */
 {if ((c6%4)==3) return 1; return (((2*c4+c6)%8)==3) ? 1:-1;}
 if ((v4>=7) && (v6==7) && (vd==8)) return -1; /* cases 16,27 */
 if ((v4==4) && (v6==6) && (vd==10)) /* cases 17,32 */
 {if ((c6%4)==1) return 1; t=(c4-2*c6)%64; if (t<0) t+=64;
  if ((t==3) || (t==19)) return 1; return -1;}
 if ((v4==7) && (v6==9) && (vd==12)) /* case 18 */
 {if (((c4%4)==1) && (((c6%8)==1)||((c6%8)==7))) return 1;
  if (((c4%4)==3) && (((c6%8)==1)||((c6%8)==3))) return 1; return -1;}
 if ((v4==7) && (v6==10) && (vd==14)) /* case 19 */
 {if (((c4%4)==1) && ((((c4*c6)%8)==5)||(((c4*c6)%8)==7))) return 1;
  return (((c4%4)==3) && ((((c4*c6)%8)==1)||(((c4*c6)%8)==7))) ? 1:-1;}
 if ((v4==7) && (v6==11) && (vd==15)) /* case 20 */
 {if ((((2*c6+c4)%8)==1) || (((2*c6+c4)%8)==3)) return 1; return -1;}
 if ((v4==7) && (v6>=12) && (vd==15)) return (((c4%8)==5)||((c4%8)==7)) ? 1:-1;
 if ((v4==4) && (v6==6) && (vd==11)) return ((c6%8)!=7) ? 1:-1; /* 22,35! */
 if ((v4>=8) && (v6==9) && (vd==12)) return ((c6%4)==1) ? 1:-1; /* case 23 */
 if ((v4>=8) && (v6==10) && (vd==14)) return ((c6%4)==1) ? 1:-1; /* case 24 */
 if ((v4==4) && (v6==6) && (vd==9)) /* case 26 */
 {if ((c6%8)==7) return 1; if (((2*c6+c4)%32)==11) return 1; return -1;}
 if ((v4==6) && (v6==8) && (vd==10)) return (((c4*c6)%4)==3) ? 1:-1; /* 28 */
 if ((v4>=7) && (v6==8) && (vd==10)) return ((c6%4)==1) ? 1:-1; /* case 29 */
 if ((v4==6) && (v6==10) && (vd==12)) /* cases 33,36 */
 {if ((c4%4)==3) return 1; t=(c4+4*c6)%16;
  if ((t==9) || (t==13)) return 1; return -1;}
 if ((v4==6) && (v6>=11) && (vd==12)) /* cases 33,37 */
 {if ((c4%4)==3) return -1; /* case 33 */
  QD_mulall(wmax,Ec6,0.0009765625,T); t=(c4-4*QD_modi(T,16))%16;
  if (t<0) t+=16; if ((t==5) || (t==9)) return 1; return -1;}
 if ((v4==6) && (v6==9) && (vd==13)) /* case 34 */
 {if ((c4%16)==11) return 1; if (((c4+4*c6)%16)==3) return 1; return -1;}
 if ((v4==4) && (v6==6) && (vd==12) && ((c6%4)==1)) return -1; /* case 38 */
 if ((v4==6) && (v6==9) && (vd==14)) return ((D%4)==(c6%4)) ? 1:-1; /* 39 */
 if ((v4==6) && (v6==9) && (vd==15)) return ((D%4)==3) ? 1:-1; /* case 40 */
 if ((v4==6) && (v6==9) && (vd>=16)) return ((c6%4)==3) ? 1:-1; /* case 41 */
 errorit("Problem in local_rootno_at2"); return 0;}

static int local_rootno_at3(int b)
{int c4,c6,D,v4,v6,vd,i; QD T; if (DEBUG) printf("local_rootno_at3\n");
 if (b==1) return -kronll(-QD_modll(Ec6,3.0),3); if (b==31) return -1;
 v4=QD_valuation(Ec4,3.0); v6=QD_valuation(Ec6,3.0);
 vd=QD_valuation(Edisc,3.0);
 if (v4>100) v4=100; if (v6>100) v6=100; if (vd>100) vd=100;
 QD_copy(wmax,Ec4,T); for (i=0;i<v4;i++) QD_div1(wmax,T,3.0,T);
 QD_round(wmax,T,T); c4=QD_modi(T,9.0);
 QD_copy(wmax,Ec6,T); for (i=0;i<v6;i++) QD_div1(wmax,T,3.0,T);
 QD_round(wmax,T,T); c6=QD_modi(T,9.0);
 QD_copy(wmax,Edisc,T); for (i=0;i<vd;i++) QD_div1(wmax,T,3.0,T);
 QD_round(wmax,T,T); D=QD_modi(T,9.0);
 if (DEBUG) printf("r3 c4:%i c6:%i D:%i v4:%i v6:%i vD:%i\n",c4,c6,D,v4,v6,vd);
 if ((v4>=2) && (v6==3) && (vd==3)) /* cases 1 and 5 */
 {if ((COND0%27)!=0) return 1;
  if ((c6==4) || (c6==7) || (c6==8)) return 1; return -1;}
 if ((v4==2) && (v6==4) && (vd==3)) {if ((c4%3)!=(c6%3)) return 1; return -1;}
 if ((v4==2) && (v6==3) && (vd==4)) return 1; /* case 3 */
 if ((v4>=2) && (v6==4) && (vd==5)) {if ((c6%3)==2) return 1; return -1;}
 if ((v4==2) && (v6>=5) && (vd==3)) return 1; /* case 6 */
 if ((v4==2) && (v6==3) && (vd==5))  {if ((D%3)==(c6%3)) return 1; return -1;}
 if ((v4==3) && (v6==5) && (vd==6)) {if ((c4%3)==2) return 1; return -1;}
 if ((v4>=4) && (v6==5) && (vd==7)) {if ((c6%3)==2) return 1; return -1;}
 if ((v4==4) && (v6==6) && (vd==9)) /* cases 10,15 */
 {if ((COND0%27)!=0) return 1; /* case 15 */
  if (((c6%9)==4) || ((c6%9)==8)) return 1; return -1;} /* case 10 */
 if ((v4>=5) && (v6==6) && (vd==9)) /* cases 11,15 */
 {if ((COND0%27)!=0) return 1; /* case 15 */
  if (((c6%9)==1) || ((c6%9)==2)) return 1; return -1;} /* case 11 */
 if ((v4==4) && (v6==7) && (vd==9)) {if ((c6%3)==2) return 1; return -1;}
 if ((v4==4) && (v6==6) && (vd==10)) /* case 13 */
 {if (((c6%9)==2) || ((c6%9)==7)) return 1; return -1;}
 if ((v4>=5) && (v6==7) && (vd==11)) {if ((c6%3)==1) return 1; return -1;}
 if ((v4==4) && (v6>=8) && (vd==9)) return 1; /* case 16 */
 if ((v4==4) && (v6==6) && (vd==11)) {if ((c6%3)==1) return 1; return -1;}
 if ((v4==5) && (v6==8) && (vd==12)) return 1; /* case 18 */
 if ((v4>=6) && (v6==8) && (vd==13)) {if ((c6%3)==1) return 1; return -1;}
 if ((v4==2) && (v6==3) && (vd==6)) return -1;
 if ((v4==3) && (v6>=6) && (vd==6)) return -1;
 errorit("Problem in local_rootno_at3"); return 0;}

static int local_rootno1(int i,llint p)
{int b=badprimetype[i],e;
 if (p==2) return local_rootno_at2(b); if (p==3) return local_rootno_at3(b);
 if (b==1) return -kronll(-QD_modll(Ec6,p),p); if (b==31) return kronll(-1,p);
 e=QD_valuation(Edisc,p)%12; if (e>6) e=12-e; e=12/e; /* e is 2,3,4,6 */
 if (e==4) return kronll(-2,p); if (e==3) return kronll(-3,p);
 return kronll(-1,p);}

int local_rootno(int i,llint p,int sp)
{int w1,b=badprimetype[i];
 if ((sp&1)==0) return 1; if (sp==1) return local_rootno1(i,p);
 w1=local_rootno1(i,p);
 if (p==2)
 {if (b==1) return w1; if ((b==29) || (b==30)) return ((sp&3)==3) ? 1:w1;
  if (((tame_local_conductor(b,sp)/2)&1)==0) w1=1;
  if (((sp&7)==3) && (fp2&1)) return -w1; return w1;}
 if (b==31) return ((sp&3)==3) ? 1:w1;
 if (p==3)
 {if ((b==4) || (b==14) || (b==18)) return 1;
  if (b==2) return ((sp&3)==3) ? 1:-1; if (b==1) return w1;
  if (((b==15)||(b==19)) && (w1==-1)) return ((sp&3)==3) ? 1:-1;
  sp=sp%12; if ((b==15)||(b==19)) return ((sp==3)||(sp==5)||(sp==7)) ? -1:1;
  if (w1==1) return (sp==5) ? -1:1; return (sp==11) ? 1:-1;}
 if (w1==1) return 1; if (b==1) return -1;
 if (b!=3) {if ((sp&3)==1) return -1; return 1;}
 if ((tame_local_conductor(b,sp)/2)&1) return -1; return 1;}

int global_rootno(int m)
{int i,e; if (((m&7)==1) || ((m&7)==3)) e=-1; else e=1;
 for (i=0;badprimes[i]!=0;i++) e*=local_rootno(i,badprimes[i],m); return e;}
