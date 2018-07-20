#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

static void prep_moddeg_at2(int bpt)
{double ap; int c4,c6,D,M=1<<29; QDpoly v;
 if ((bpt==3) || (bpt==20) || (bpt==21) || (bpt==23))
 {EVEN_TOP=3; EVEN_BOTTOM=2;}
 if (bpt==22) {EVEN_TOP=1; EVEN_BOTTOM=2;}
 if ((bpt==1) || (bpt==3) || (bpt==8) || (bpt==9) || (bpt==24)) return;
 c4=QD_modi(Ec4,M); c6=QD_modi(Ec6,M); D=QD_modi(Edisc,M);
 if ((bpt==10) || (bpt==22) || (bpt==23))
 {if (((c4&63)!=31) || ((c6&255)!=0))
  {QD_mul1(wmax,VOLUME,2.0,VOLUME); QD_mul1(wmax,TW_EFF,2.0,TW_EFF);} return;}
 if (bpt==16)
 {CONDTW/=16; QD_mul1(wmax,VOLUME,4.0,VOLUME);
  euler_factor_bad(2,bpt,2,&v); ap=Round(Sqrt(-v.coeff[1][0]+2.0));
  QD_mul1(wmax,TW_EFF,72.0-8.0*ap*ap,TW_EFF); delQDpoly(&v); return;}
 if (bpt==17)
 {CONDTW/=64;
  if (((c6&127)==64) && ((c4&63)==0)) QD_mul1(wmax,VOLUME,2.0,VOLUME);
  else QD_mul1(wmax,VOLUME,8.0,VOLUME);
  euler_factor_bad(2,bpt,2,&v); ap=Round(Sqrt(-v.coeff[1][0]+2.0));
  QD_mul1(wmax,TW_EFF,9.0-ap*ap,TW_EFF); delQDpoly(&v);
  if (((c6&127)==64) && ((c4&63)==0)) QD_mul1(wmax,TW_EFF,16.0,TW_EFF);
  else QD_mul1(wmax,TW_EFF,64.0,TW_EFF); return;}
 if (bpt==20) {CONDTW/=4; QD_mul1(wmax,TW_EFF,4.0,TW_EFF); return;}
 if (bpt==21)
 {QD_mul1(wmax,VOLUME,2.0,VOLUME); CONDTW/=16;
  QD_mul1(wmax,TW_EFF,32.0,TW_EFF); return;}
 if (bpt==25) {CONDTW/=2; QD_mul1(wmax,TW_EFF,2.0,TW_EFF); return;}
 if (bpt==26)
 {QD_mul1(wmax,VOLUME,2.0,VOLUME); CONDTW/=8;
  QD_mul1(wmax,TW_EFF,16.0,TW_EFF); return;}
 if (bpt==27)
 {if (((c4&63)==0) && ((c6&255)==0))
  {QD_mul1(wmax,VOLUME,2.0,VOLUME); QD_mul1(wmax,TW_EFF,2.0,TW_EFF);} return;}
 if (bpt==29) /* v2(N)=6 multiplicative twist */
 {CONDTW/=32; QD_mul1(wmax,VOLUME,4.0,VOLUME);
  QD_mul1(wmax,TW_EFF,96.0,TW_EFF); return;}
 if (bpt==30) /* v2(N)=4 multiplicative twist */
 {CONDTW/=8; QD_mul1(wmax,VOLUME,4.0,VOLUME);
  QD_mul1(wmax,TW_EFF,24.0,TW_EFF); return;}
 errorit("bad prime type at 2 in prep_moddeg_at2");}

static llint is_square(llint x)
{llint y; if (x<0) return 0;
 y=(llint) Round(Sqrt((double) x)); if (y*y==x) return y; return 0;}

static void isogcheck3(int n)
{int i,u=(n+3),v=0; QD U,V,W,R; llint p,q;
 if (DEBUG) printf("isogcheck3 %i\n",n); q=(llint) n; q=q*q+9*q+27;
 for (i=0;badprimes[i]!=0;i++)
 {p=badprimes[i]; if (((CONDTW%p)==0) && ((n%p)==0) && ((p%6)==1)) return;
  if ((q%p)==0) v++;} if (v!=1) return;
 QD_copy(wmax,QD_zero,U); U[0]=(double) u; QD_powi(wmax,U,4,V);
 QD_mul1(wmax,U,24.0,W); QD_sub(wmax,V,W,R); QD_round(wmax,R,R);
 QD_sub(wmax,Etw4,R,R); if (R[0]!=0.0) return;
 QD_powi(wmax,U,6,V); QD_powi(wmax,U,3,W); QD_mul1(wmax,W,36.0,W);
 QD_sub(wmax,W,V,W); QD_copy(wmax,QD_zero,V); V[0]=-216.0; QD_add(wmax,W,V,R);
 QD_sub(wmax,Etw6,R,R); if (R[0]!=0.0) return;
 MANIN_TWIST=9; if (Ec6[0]==Etw6[0]) MANIN=9; return;}

static void check_manin()
{QD T,U,R; llint p,bp[32]; int i,n,tau=0; if (DEBUG) printf("check_manin\n");
 MANIN=1; MANIN_TWIST=1;
 if (CONDTW==11)
 {if (Etw4[0]==16.0) {MANIN_TWIST=25; if (Ec4[0]==16.0) MANIN=25;} return;}
 if (CONDTW==15)
 {if (Etw4[0]==1.0) {MANIN_TWIST=16; if (Ec4[0]==1.0) MANIN=16;}
  if (Etw4[0]==241.0) {MANIN_TWIST=4; if (Ec4[0]==241.0) MANIN=4;}
  if (Etw4[0]==3841.0) {MANIN_TWIST=4; if (Ec4[0]==3841.0) MANIN=4;} return;}
 if (CONDTW==17)
 {if (Etw6[0]==-81.0) {MANIN_TWIST=16; if (Ec6[0]==-81.0) MANIN=16;}
  if (Etw4[0]==273.0) {MANIN_TWIST=4; if (Ec4[0]==273.0) MANIN=4;}
  if (Etw4[0]==4353.0) {MANIN_TWIST=4; if (Ec4[0]==4353.0) MANIN=4;} return;}
 if (CONDTW==20)
 {if (Etw4[0]==64.0) {MANIN_TWIST=4; if (Ec4[0]==64.0) MANIN=4;}
  if (Etw4[0]==1984.0) {MANIN_TWIST=4; if (Ec4[0]==1984.0) MANIN=4;} return;}
 if (CONDTW==24)
 {if (Etw4[0]==-32.0) {MANIN_TWIST=4; if (Ec4[0]==-32.0) MANIN=4;} return;}
 if (CONDTW==32)
 {if (Ec6[0]==0.0)
  {MANIN_TWIST=4; if ((COND0==32) && (QD_modi(Ec4,32)==16)) MANIN=4;}
  else {MANIN_TWIST=4; if (COND0==32) MANIN=4;} return;}
 if (CONDTW==40)
 {if (Etw4[0]==96.0) {MANIN_TWIST=4; if (Ec4[0]==96.0) MANIN=4;} return;}
 if (CONDTW==64)
 {if (Ec6[0]==0.0)
  {MANIN_TWIST=4; if ((COND0==64) && (QD_modi(Ec4,32)==16)) MANIN=4;} return;}
 if (CONDTW==128)
 {if (Etw4[0]==112.0)
  {MANIN_TWIST=4; if ((COND0==128) && (QD_modi(Ec4,32)==16)) MANIN=4;} return;}
 /* 3-isog case */
 for (n=-10;n<=10;n++) {if (EtwD[0]==(double) (n*(n*n+9*n+27))) isogcheck3(n);}
 QD_copy(wmax,QD_zero,R);
 if (EtwD[0]>0.0) R[0]=Round(Root(EtwD[0],3)-3);
 else R[0]=-Round(Root(-EtwD[0],3)+3);
 QD_sqr(wmax,R,T); QD_mul1(wmax,R,9.0,U); QD_add(wmax,T,U,T);
 QD_copy(wmax,QD_zero,U); U[0]=27.0; QD_add(wmax,T,U,T);
 QD_mul(wmax,T,R,T); QD_round(wmax,T,T); QD_sub(wmax,T,EtwD,U);
 if (U[0]==0.0) isogcheck3((int) R[0]);
 for (i=0;badprimes[i]!=0;i++)
   if ((CONDTW%badprimes[i])==0) bp[tau++]=badprimes[i];
 if (tau>2) return;
 if (tau==1) /* original neumann-setzer case, don't check c6 */
 {if (!(is_square(CONDTW-64))) return;
  if ((CONDTW-16)==(llint) Etw4[0])
  {MANIN_TWIST=4; if (Ec4[0]==Etw4[0]) MANIN=4; return;} return;}
 if (CONDTW&1) /* Third NS case bad for curves with big Etw6 */
 {if (Etw4[0]<-48.0) return; p=(llint) Round(Sqrt(Etw4[0]+48.0));
  if ((p&3)!=3) p=-p; p-=8; if ((p*p+16*p+16)!=(llint) Etw4[0]) return;
  if ((p*p*p+24*p*p+120*p-64)!=(llint) Etw6[0]) return;
  if (((bp[0]&3)!=3) && ((bp[1]&3)!=3)) return;
  if (((p%bp[0])!=0) && ((p%bp[1])!=0)) return;
  if ((((p+16)%bp[0])!=0) && (((p+16)%bp[1])!=0)) return;
  MANIN_TWIST=4; if (Ec4[0]==Etw4[0]) MANIN=4; return;}
 if ((CONDTW&7)!=4) return; if (!(is_square(CONDTW/4-4))) return;
 if ((4*CONDTW-16)==(llint) Etw4[0]) /* second NS case, again ignore c6 */
 {MANIN_TWIST=4; if (Ec4[0]==Etw4[0]) MANIN=4;}}

void prepare_moddeg(char *IN)
{int i,s,bpt,v4,v6,vd; llint p; QDpoly v; double ap,P; QD C;
 QD_mul(wmax,REAL_PERIOD,IMAG_PERIOD,VOLUME);
 QD_mul(wmax,VOLUME,QD_twopi,VOLUME); EVEN_TOP=1; EVEN_BOTTOM=1;
 QD_copy(wmax,QD_one,TW_EFF); CONDTW=COND0;
 for (i=0;badprimes[i]!=0;i++)
 {p=badprimes[i]; bpt=badprimetype[i]; P=(double) p;
  if (p>=3)
  {if (bpt==2)
   {QD_mul1(wmax,VOLUME,P,VOLUME); CONDTW/=(p*p);
    euler_factor_bad(p,bpt,2,&v); ap=Round(Sqrt(-v.coeff[1][0]+P));
    QD_mul1(wmax,TW_EFF,P-1.0,TW_EFF); QD_mul1(wmax,TW_EFF,P+1.0-ap,TW_EFF);
    QD_mul1(wmax,TW_EFF,P+1.0+ap,TW_EFF); delQDpoly(&v);}
   else if (bpt==31)
   {QD_mul1(wmax,VOLUME,P,VOLUME); CONDTW/=p;
    QD_mul1(wmax,TW_EFF,P+1.0,TW_EFF); QD_mul1(wmax,TW_EFF,P-1.0,TW_EFF);}
   else if (bpt!=1)
   {v4=QD_valuation(Ec4,P); v6=QD_valuation(Ec6,P); vd=QD_valuation(Edisc,P);
    if ((v4>=2) && (v6>=3) && ((p!=3) || ((vd>=6) && (v6!=5))))
    {QD_mul1(wmax,VOLUME,P,VOLUME); QD_mul1(wmax,TW_EFF,P,TW_EFF);}}
   if ((bpt>=3) && (bpt<=6))
   {if ((p%bpt)==1) {EVEN_TOP*=(p-1); EVEN_BOTTOM*=p;}
    else {EVEN_TOP*=(p+1); EVEN_BOTTOM*=p;}}
   if ((bpt==14) || (bpt==15)) {EVEN_TOP*=2; EVEN_BOTTOM*=3;}
   if ((bpt==18) || (bpt==19)) {EVEN_TOP*=4; EVEN_BOTTOM*=3;}}
  if (p==2) prep_moddeg_at2(bpt);}
 IN[0]='2'; IN[1]='w'; IN[2]='0'; IN[3]='p';
 QD_copy(wmax,QD_zero,C); C[0]=(double) CONDTW; QD_div(wmax,C,VOLUME,SCALE);
 if (CM_CASE==-3)
 {QD_mul(wmax,SCALE,QD_pi,SCALE); C[0]=27.0;
  QD_sqrt(wmax,C,C); QD_div(wmax,SCALE,C,SCALE); HECKE=TRUE;}
 if (CM_CASE==-4)
 {QD_mul(wmax,SCALE,QD_pi,SCALE); QD_mul2n(wmax,SCALE,-2,SCALE); HECKE=TRUE;}
 s=4+Ceil(Log10(SCALE[0])); if (s<6) s=6; if ((s>12) && (s<=16)) s=17;
 if (s>28) {s=28; printf("*WARNING* reducing moddeg precision to 28");}
 if (s<10) {IN[4]='0'+s; IN[5]=0;}
 else {IN[4]='0'+(s/10); IN[5]='0'+(s%10); IN[6]=0;}
 if (VERBOSE)
 {printf("Twist effect is "); QD_intout_noplus(TW_EFF);
  printf(", Even valuation gives %lli/%lli\n",EVEN_TOP,EVEN_BOTTOM);}
 if (VERBOSE) printf("Twist conductor is %lli\n",CONDTW);
 if (VERBOSE) printf("Using precision of %i digits\n",s);
 check_manin(); if (MANIN!=1) printf("Manin effect for curve is %i\n",MANIN);
 if (MANIN_TWIST!=1) printf("Manin effect for twist is %i\n",MANIN_TWIST);}

llint postpare_moddeg()
{QD C,D,W; double L; llint N=CONDTW/10;
 QD_copy(wmax,QD_zero,C); QD_copy(wmax,QD_zero,W);
 C[0]=(double) CONDTW; QD_div(wmax,C,VOLUME,C);
 while (1)
 {get_weight(N,W,0); QD_mul(wmax,C,W,D); L=Log((double) N); L*=L;
  QD_mul1(wmax,D,L,D); QD_mul1(wmax,D,(double) N,D);
  QD_mul1(wmax,D,(double) N,D); if (D[0]<0.5) break; N=(11*N)/10;} return(N);}
