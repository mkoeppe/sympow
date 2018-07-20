#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static void from_c4c6()
{int b2; QD T,U; if (DEBUG) printf("from_c4c6\n"); QD_copy(wmax,QD_zero,Eb2);
 b2=-QD_modi(Ec6,12); if (b2<-5) b2+=12; Eb2[0]=(double) b2;
 QD_sqr(wmax,Eb2,T); QD_sub(wmax,T,Ec4,T); QD_div1(wmax,T,24.0,Eb4);
 QD_powi(wmax,Eb2,3,T); QD_add(wmax,T,Ec6,T);
 QD_mul(wmax,Eb2,Eb4,U); QD_mul1(wmax,U,36.0,U);
 QD_sub(wmax,U,T,T); QD_div1(wmax,T,216.0,Eb6);
 QD_copy(wmax,QD_zero,Ea1); QD_copy(wmax,QD_zero,Ea3);
 Ea1[0]=(double) QD_modi(Eb2,2); Ea3[0]=(double) QD_modi(Eb6,2);
 QD_sub(wmax,Eb2,Ea1,Ea2); QD_mul2n(wmax,Ea2,-2,Ea2);
 QD_sub(wmax,Eb6,Ea3,Ea6); QD_mul2n(wmax,Ea6,-2,Ea6);
 QD_mul(wmax,Ea1,Ea3,T); QD_sub(wmax,Eb4,T,Ea4); QD_mul2n(wmax,Ea4,-1,Ea4);
 QD_powi(wmax,Ec4,3,T); QD_sqr(wmax,Ec6,U);
 QD_sub(wmax,T,U,T); QD_div1(wmax,T,1728.0,U); QD_round(wmax,U,Edisc);
 if (Edisc[0]==0.0) errorit("Curve is singular");}

static void from_a1a2a3a4a6()
{QD T,U,V; if (DEBUG) printf("from_a1a2a3a4a6\n");
 QD_sqr(wmax,Ea1,T); QD_mul2n(wmax,Ea2,2,U); QD_add(wmax,T,U,Eb2);
 QD_mul(wmax,Ea1,Ea3,T); QD_mul2n(wmax,Ea4,1,U); QD_add(wmax,T,U,Eb4);
 QD_sqr(wmax,Ea3,T); QD_mul2n(wmax,Ea6,2,U); QD_add(wmax,T,U,Eb6);
 QD_sqr(wmax,Eb2,T); QD_mul1(wmax,Eb4,24.0,U); QD_sub(wmax,T,U,Ec4);
 QD_mul(wmax,Eb2,Eb4,T); QD_mul1(wmax,T,36.0,V); QD_powi(wmax,Eb2,3,T);
 QD_mul1(wmax,Eb6,216.0,U); QD_add(wmax,T,U,T); QD_sub(wmax,V,T,Ec6);
 QD_powi(wmax,Ec4,3,T); QD_sqr(wmax,Ec6,U);
 QD_sub(wmax,T,U,T); QD_div1(wmax,T,1728.0,U); QD_round(wmax,U,Edisc);
 if (Edisc[0]==0.0) errorit("Curve is singular");}

static void minimise_at(QD g,llint p)
{QD s4,s6,P; int d,a,b; if (DEBUG) printf("minimise at %lli\n",p);
 QD_copy(wmax,QD_zero,P); P[0]=(double) p; d=QD_valuation(g,(double) p)/12;
 if ((p==3) && (QD_valuation(Ec6,3)==6*d+2)) d--;
 QD_powi(wmax,P,4*d,s4); QD_powi(wmax,P,6*d,s6);
 QD_div(wmax,Ec4,s4,Ec4); QD_div(wmax,Ec6,s6,Ec6);
 QD_round(wmax,Ec4,Ec4); QD_round(wmax,Ec6,Ec6);
 if (p!=2) return; a=QD_modi(Ec4,16); b=QD_modi(Ec6,32);
 if (((b&3)!=3) && !((a==0) && ((b==0) || (b==8))))
 {QD_mul1(wmax,Ec4,16,Ec4); QD_mul1(wmax,Ec6,64,Ec6);}}

static void minimal_model_c4c6()
{QD T,g; LIST L; int i=0; if (DEBUG) printf("minimal_model_c4c6\n");
 QD_sqr(wmax,Ec6,T); QD_intgcd(T,Edisc,g); L.p[0]=0; QD_factor(g,&L);
 for (i=0;L.p[i]!=0;i++) minimise_at(g,L.p[i]);
 C4C6LL=FALSE;
 if ((Abs(Ec4[0])<3.4e17) && (Abs(Ec6[0])<1.7e17))
 {C4C6LL=TRUE; Ec4ll=((llint) Ec4[0])+((llint) Ec4[1]);
  Ec6ll=((llint) Ec6[0])+((llint) Ec6[1]);}}

static int is_nonmintwist()
{llint p; int v4,v6,vd; int i;
 for (i=0;badprimes[i]!=0;i++)
 {p=badprimes[i];
  if (p>3)
  {if ((QD_valuation(Ec4,p)>=2) && (QD_valuation(Ec6,p)>=3)) return TRUE;}
  if (p==3)
  {v4=QD_valuation(Ec4,3.0); v6=QD_valuation(Ec6,3.0);
   vd=QD_valuation(Edisc,3.0);
   if ((v4>=2) && (v6>=3) && (vd>=6) && (v6!=5)) return TRUE;}
  if (p==2)
  {v4=QD_valuation(Ec4,2.0); v6=QD_valuation(Ec6,2.0);
   vd=QD_valuation(Edisc,2.0);
   if ((v4>=4) && (v6>=5) && (Ec6[0]<0.0)) return TRUE; /* twist by -1 */
   if ((v4==4) && (v6==6) && (vd>=12)) return TRUE; /* twist by -1 drops */
   if ((v4>=6) && (v6==6) && (vd>=6)) return TRUE; /* twist by pm2 */
   if ((v4>=6) && (v6>=8) && (vd>=6)) return TRUE;}} /* twist by pm2 */
 return FALSE;}

void mintwist()
{llint p; int v4,v6,vd; int i; double P;
 QD_copy(wmax,Ec4,Etw4); QD_copy(wmax,Ec6,Etw6); QD_copy(wmax,Edisc,EtwD);
 for (i=0;badprimes[i]!=0;i++)
 {p=badprimes[i]; P=(double) p; if ((p&3)==3) P=-P;
  if (p>3)
  {if ((QD_valuation(Ec4,p)>=2) && (QD_valuation(Ec6,p)>=3))
   {QD_div1(wmax,Etw4,P,Etw4); QD_div1(wmax,Etw4,P,Etw4);
    QD_div1(wmax,Etw6,P,Etw6); QD_div1(wmax,Etw6,P,Etw6);
    QD_div1(wmax,Etw6,P,Etw6);
    QD_div1(wmax,EtwD,P,EtwD); QD_div1(wmax,EtwD,P,EtwD);
    QD_div1(wmax,EtwD,P,EtwD); QD_div1(wmax,EtwD,P,EtwD);
    QD_div1(wmax,EtwD,P,EtwD); QD_div1(wmax,EtwD,P,EtwD);}}
  if (p==3)
  {v4=QD_valuation(Ec4,3.0); v6=QD_valuation(Ec6,3.0);
   vd=QD_valuation(Edisc,3.0);
   if ((v4>=2) && (v6>=3) && (vd>=6) && (v6!=5))
   {QD_div1(wmax,Etw4,9.0,Etw4); QD_div1(wmax,Etw6,-27.0,Etw6);
    QD_div1(wmax,EtwD,729.0,EtwD);}}
  if (p==2)
  {v4=QD_valuation(Ec4,2.0); v6=QD_valuation(Ec6,2.0);
   vd=QD_valuation(Edisc,2.0);
   if ((v4==4) && (v6==6) && (vd>=12)) /* twist by -1 to c4 odd */
   {QD_div1(wmax,Etw4,16.0,Etw4); QD_div1(wmax,Etw6,-64.0,Etw6);
    if (QD_modi(Etw6,4)!=3) QD_neg(wmax,Etw6,Etw6);
    QD_div1(wmax,EtwD,4096.0,EtwD);}
   else if ((v4>=8) && (v6==9) && (vd==12)) /* tw by -1 to c6 8 mod32 drop */
   {QD_div1(wmax,Etw4,16.0,Etw4); QD_div1(wmax,Etw6,-64.0,Etw6);
    if (QD_modi(Etw6,32)!=8) QD_neg(wmax,Etw6,Etw6);
    QD_div1(wmax,EtwD,4096.0,EtwD);}
   else if ((v4==6) && (v6==9) && (vd>=18)) /* twist by +-2 to c4 odd */
   {QD_div1(wmax,Etw4,64.0,Etw4); QD_div1(wmax,Etw6,512.0,Etw6);
    if (QD_modi(Etw6,4)!=3) QD_neg(wmax,Etw6,Etw6);
    QD_div1(wmax,EtwD,262144.0,EtwD);}
   else if ((v4>=6) && (v6==6) && (vd==6)) /* tw by +-2 to c6 8 mod 32 dir */
   {QD_div1(wmax,Etw4,4.0,Etw4); QD_div1(wmax,Etw6,8.0,Etw6);
    if (QD_modi(Etw6,32)!=8) QD_neg(wmax,Etw6,Etw6);
    QD_div1(wmax,EtwD,64.0,EtwD);}
   else if ((v4>=6) && (v6>=8)) /* twist by +-2 */
   {QD_div1(wmax,Etw4,4.0,Etw4); QD_div1(wmax,Etw6,8.0,Etw6);
    QD_div1(wmax,EtwD,64.0,EtwD);}}}
 QD_round(wmax,Etw4,Etw4); QD_round(wmax,Etw6,Etw6); QD_round(wmax,EtwD,EtwD);}

static void check_cm()
{QD T,q,U; int n=0; if (DEBUG) printf("check_cm\n");
 QD_powi(wmax,Ec4,3,T); QD_round(wmax,T,T);
 QD_div(wmax,T,Edisc,q); QD_round(wmax,q,q); QD_mul(wmax,q,Edisc,U);
 QD_round(wmax,U,U); QD_sub(wmax,T,U,U); if (U[0]!=0.0) return;
 QD_div(wmax,T,Edisc,T); QD_round(wmax,T,T);
 if (T[0]==0.0) CM_CASE=-3; else if (T[0]==1728.0) CM_CASE=-4;
 else if (T[0]==-3375.0)
 {CM_CASE=-7; QD_div1(wmax,Ec6,1323.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==8000.0)
 {CM_CASE=-8; QD_div1(wmax,Ec6,-1792.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==-32768.0)
 {CM_CASE=-11; QD_div1(wmax,Ec6,-6776.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==54000.0)
 {CM_CASE=-12; QD_div1(wmax,Ec6,-19008.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==287496.0)
 {CM_CASE=-16; QD_div1(wmax,Ec6,12096.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==-884736.0)
 {CM_CASE=-19; QD_div1(wmax,Ec6,-77976.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==-12288000.0)
 {CM_CASE=-27; QD_div1(wmax,Ec6,-54648.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==16581375.0)
 {CM_CASE=-28; QD_div1(wmax,Ec6,75411.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==-884736000.0)
 {CM_CASE=-43; QD_div1(wmax,Ec6,-8387064.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if (T[0]==-147197952000.0)
 {CM_CASE=-67; QD_div1(wmax,Ec6,-210408408.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}
 else if ((T[0]==-262537412640768000.0) && (T[1]==0.0))
 {CM_CASE=-163; QD_div1(wmax,Ec6,-1066294102104.0,T);
  if (T[0]<0.0) {n=1; QD_neg(wmax,T,T);}
  CM_TWIST=(int) Round(Root(T[0],3)); if (n==1) CM_TWIST=-CM_TWIST;}}

static void read_curve(char *S,int tw)
{int i=0,n=0; QD T; if (DEBUG) printf("read_curve %s %i\n",S,tw);
 QD_copy(wmax,QD_zero,Ea1); QD_copy(wmax,QD_zero,T); ASSERT(S[0]=='[');
 QD_copy(wmax,QD_zero,Ea2); QD_copy(wmax,QD_zero,Ea3);
 QD_copy(wmax,QD_zero,Ea4); QD_copy(wmax,QD_zero,Ea6);
 i++; if (S[i]=='-') {n=1; i++;}
 while (S[i]!=',')
 {ASSERT(ISA_NUMBER(S[i])); T[0]=(double) (S[i]-'0');
  QD_mul1(wmax,Ea1,10,Ea1); QD_add(wmax,Ea1,T,Ea1); i++;}
 if (n==1) QD_neg(wmax,Ea1,Ea1); n=0;
 i++; if (S[i]=='-') {n=1; i++;}
 while (S[i]!=',')
 {ASSERT(ISA_NUMBER(S[i])); T[0]=(double) (S[i]-'0');
  QD_mul1(wmax,Ea2,10,Ea2); QD_add(wmax,Ea2,T,Ea2); i++;}
 if (n==1) QD_neg(wmax,Ea2,Ea2); n=0;
 i++; if (S[i]=='-') {n=1; i++;}
 while (S[i]!=',')
 {ASSERT(ISA_NUMBER(S[i])); T[0]=(double) (S[i]-'0');
  QD_mul1(wmax,Ea3,10,Ea3); QD_add(wmax,Ea3,T,Ea3); i++;}
 if (n==1) QD_neg(wmax,Ea3,Ea3); n=0;
 i++; if (S[i]=='-') {n=1; i++;}
 while (S[i]!=',')
 {ASSERT(ISA_NUMBER(S[i])); T[0]=(double) (S[i]-'0');
  QD_mul1(wmax,Ea4,10,Ea4); QD_add(wmax,Ea4,T,Ea4); i++;}
 if (n==1) QD_neg(wmax,Ea4,Ea4); n=0;
 i++; if (S[i]=='-') {n=1; i++;}
 while (S[i]!=']')
 {ASSERT(ISA_NUMBER(S[i])); T[0]=(double) (S[i]-'0');
  QD_mul1(wmax,Ea6,10,Ea6); QD_add(wmax,Ea6,T,Ea6); i++;}
 if (n==1) QD_neg(wmax,Ea6,Ea6); n=0; from_a1a2a3a4a6();
 QD_mul1(wmax,Ec4,(double) tw,Ec4); QD_mul1(wmax,Ec4,(double) tw,Ec4);
 QD_mul1(wmax,Ec6,(double) tw,Ec6); QD_mul1(wmax,Ec6,(double) tw,Ec6);
 QD_mul1(wmax,Ec6,(double) tw,Ec6);
 QD_mul1(wmax,Edisc,(double) tw,Edisc); QD_mul1(wmax,Edisc,(double) tw,Edisc);
 QD_mul1(wmax,Edisc,(double) tw,Edisc); QD_mul1(wmax,Edisc,(double) tw,Edisc);
 QD_mul1(wmax,Edisc,(double) tw,Edisc); QD_mul1(wmax,Edisc,(double) tw,Edisc);
 minimal_model_c4c6(); from_c4c6();
 if (Abs(Ec4[0])>1e18) errorit("c4 invariant is too large");
 if (Abs(Ec6[0])>1e27) errorit("c6 invariant is too large");}

static void sn_curve(char *S,int tw)
{char I[256],J[96],K[256]; int wh; llint sn; QD R,T,U;
 if (DEBUG) printf("sn_curve tw:%i %s",tw,S);
 sscanf(S,"%s %s %s",J,I,K); wh=atoi(I); sn=atoll(K);
 if ((sn>10000000) || (sn<-10000000)) errorit("Too large in sn_curve");
 QD_copy(wmax,QD_zero,Ec4); QD_copy(wmax,QD_zero,Ec6);
 switch(wh)
 {case 1: Ec4[0]=(double) (sn*sn+48); if ((sn&3)==1) sn=-sn;
   Ec6[0]=(double) (sn*sn+72); QD_mul1(wmax,Ec6,(double) sn,Ec6); break;
  case 2: Ec4[0]=(double) (16*sn*sn+48); if ((sn&3)==1) sn=-sn;
   Ec6[0]=(double) (16*sn*sn+72); QD_mul1(wmax,Ec6,(double) (4*sn),Ec6); break;
  case 3: if ((sn&3)==1) sn=-sn; Ec4[0]=(double) (sn*sn+16*sn+16);
   Ec6[0]=(double) (sn*sn+16*sn-8); QD_mul1(wmax,Ec6,(double) (sn+8),Ec6);
   break;
  case 4: QD_copy(wmax,QD_zero,R); R[0]=(sn+3);
   QD_powi(wmax,R,4,T); QD_mul1(wmax,R,24.0,U); QD_sub(wmax,T,U,Ec4);
   QD_powi(wmax,R,6,T); QD_powi(wmax,R,3,U); QD_mul1(wmax,U,36.0,U);
   QD_sub(wmax,U,T,T); QD_copy(wmax,QD_zero,U); U[0]=216.0;
   QD_sub(wmax,T,U,Ec6); break;
  default: errorit("Unknown Neumann-Setzer type");}
 QD_round(wmax,Ec4,Ec4); QD_round(wmax,Ec6,Ec6); QD_round(wmax,Edisc,Edisc);
 QD_mul1(wmax,Ec4,(double) tw,Ec4); QD_mul1(wmax,Ec4,(double) tw,Ec4);
 QD_mul1(wmax,Ec6,(double) tw,Ec6); QD_mul1(wmax,Ec6,(double) tw,Ec6);
 QD_mul1(wmax,Ec6,(double) tw,Ec6);
 from_c4c6(); minimal_model_c4c6(); from_c4c6();
 if (Abs(Ec4[0])>1e18) errorit("c4 invariant is too large");
 if (Abs(Ec6[0])>1e27) errorit("c6 invariant is too large");}

static void cm_curve(char *S,int tw)
{char I[256],J[96]; int cm; if (DEBUG) printf("cm_curve tw:%i %s",tw,S);
 sscanf(S,"%s %s",J,I); cm=atoi(I); if (cm<0) cm=-cm;
 QD_copy(wmax,QD_zero,Ec4); QD_copy(wmax,QD_zero,Ec6);
 switch(cm)
 {case 3: Ec6[0]=-864.0; break; case 4: Ec4[0]=48.0; break;
  case 7: Ec4[0]=105.0; Ec6[0]=1323.0; break;
  case 8: Ec4[0]=160.0; Ec6[0]=-1792.0; break;
  case 11: Ec4[0]=352.0; Ec6[0]=-6776.0; break;
  case 12: Ec4[0]=720.0; Ec6[0]=-19008.0; break;
  case 16: Ec4[0]=528.0; Ec6[0]=-12096.0; break;
  case 19: Ec4[0]=1824.0; Ec6[0]=-77976.0; break;
  case 27: Ec4[0]=1440.0; Ec6[0]=-54648.0; break;
  case 28: Ec4[0]=1785.0; Ec6[0]=75411.0; break;
  case 43: Ec4[0]=41280.0; Ec6[0]=-8387064.0; break;
  case 67: Ec4[0]=353760.0; Ec6[0]=-210408408.0; break;
  case 163: Ec4[0]=104372160.0; Ec6[0]=-1066294102104.0; break;
  default: errorit("Unknown CM-type");}
 QD_mul1(wmax,Ec4,(double) tw,Ec4); QD_mul1(wmax,Ec4,(double) tw,Ec4);
 QD_mul1(wmax,Ec6,(double) tw,Ec6); QD_mul1(wmax,Ec6,(double) tw,Ec6);
 QD_mul1(wmax,Ec6,(double) tw,Ec6);
 from_c4c6(); minimal_model_c4c6(); from_c4c6();
 if (Abs(Ec4[0])>1e18) errorit("c4 invariant is too large");
 if (Abs(Ec6[0])>1e27) errorit("c6 invariant is too large");}

static int get_twist(char *S)
{char N[8],T[1024]; int u; sscanf(S,"%s %s",N,T); u=atoi(T);
 if (((u&3)==2) || ((u&3)==3)) return (4*u); return u;}

void curve_init(char *S,char *LSTR)
{char *T,*To,*U; int u,i=0,tw=1; QD X; LIST L;
 if (DEBUG) printf("curve_init %s %s\n",S,LSTR);
 T=malloc(1024); U=malloc(1024); strcpy(T,S); To=T;
 if (T[0]=='t') {tw=get_twist(T); sscanf(T+3,"%s",U); T+=(4+strlen(U));}
 if (T[0]=='[') read_curve(T,tw); else if (T[0]=='c') cm_curve(T,tw);
 else if (T[0]=='s') sn_curve(T,tw); else if (T[0]=='n') sn_curve(T,tw);
 else errorit("Input appears to be bogus"); free(To); free(U);
 printf("Minimal model of curve %s is [",LSTR);
 QD_intout_noplus(Ea1); printf(","); QD_intout_noplus(Ea2); printf(",");
 QD_intout_noplus(Ea3); printf(","); QD_intout_noplus(Ea4); printf(",");
 QD_intout_noplus(Ea6); printf("]\n");
 if (Edisc[0]<0.0) QD_neg(wmax,Edisc,X); else QD_copy(wmax,Edisc,X);
 L.p[0]=0; COND0=1; QD_factor(X,&L);
 for (u=0;u<32;u++) badprimes[u]=0; for (u=0;u<32;u++) badprimetype[u]=0;
 for (i=0;L.p[i]!=0;i++)
 {badprimes[i]=L.p[i]; badprimetype[i]=do_badprime(L.p[i]);
  if (VERBOSE)
  {printf("At %lli: Inertia Group is  ",badprimes[i]);
   badprimetype_output(badprimetype[i],badprimes[i]);}}
 if (VERBOSE) printf("Conductor is %lli\n",COND0);
 CM_CASE=0; check_cm(); mintwist(); NO_QT&=is_nonmintwist(); do_periods();
 if (VERBOSE>=2)
 {printf("rp: "); QD_output(4,64,REAL_PERIOD);
  printf("ip: "); QD_output(4,64,IMAG_PERIOD);}}
