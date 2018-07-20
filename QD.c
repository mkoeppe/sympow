#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

extern double QD_pi[4],QD_twopi[4],QD_sqrtpi[4];
extern double QD_log2[4],QD_one[4],QD_e[4],QD_zero[4];
double QD_log2[4]={6.931471805599452862e-01,2.319046813846299558e-17,
		   5.707708438416212066e-34,-3.582432210601811423e-50};
double QD_e[4]={2.718281828459045091e+00,1.445646891729250158e-16,
		-2.127717108038176765e-33,1.515630159841218954e-49};
double QD_one[4]={1.0,0.0,0.0,0.0};
double QD_zero[4]={0.0,0.0,0.0,0.0};
double QD_pi[4]=
{3.141592653589793116e+00,1.224646799147353207e-16,
 -2.994769809718339666e-33,1.112454220863365282e-49};
double QD_twopi[4]=
{6.283185307179586232e+00,2.449293598294706414e-16,
 -5.989539619436679332e-33,2.224908441726730563e-49};
double QD_sqrtpi[4]=
{1.772453850905516104e+00,-7.66658649982579883e-17,
 -1.305833490794542932e-33,-2.61101420878271562e-50};

// int MY_DEBUG=0;

#define QD_qsum(a,b,s,e) {double z1; s=a+b; z1=s-a; e=b-z1;}
#define QD_sum(a,b,s,e)\
 {double bb,z1,z2; s=a+b; bb=s-a; z1=s-bb; z2=a-z1; z1=b-bb; e=z1+z2;}
#define QD_diff(a,b,s,e)\
 {double bb,z1,z2; s=a-b; bb=s-a; z1=s-bb; z2=a-z1; z1=b+bb; e=z2-z1;}
#define QD_split(a,hi,lo)\
{double xt=134217729.0*a,z1; z1=xt-a; hi=xt-z1; lo=a-hi;}
#define QD_prod(a,b,p,e)\
{double ah,al,bh,bl,z1,z2,z3,z4; p=a*b; QD_split(a,ah,al); QD_split(b,bh,bl);\
 z1=ah*bh; z2=z1-p; z1=ah*bl; z3=al*bh; z4=z1+z3; z1=z2+z4; z3=al*bl; e=z1+z3;}

static void QD_renorm_3(double h,double m,double l,double *o)
{double l1,m1,m2,o0; /* renorm_2 is just QD_qsum */
 QD_qsum(m,l,m1,l1); QD_qsum(h,m1,o0,m2); o[0]=o0; QD_qsum(m2,l1,o[1],o[2]);
 if (o[0]==0.0) {o[0]=o[1]; o[1]=o[2]; o[2]=0.0;}
 if (o[0]==0.0) {o[0]=o[1]; o[1]=0.0;}}

static void QD_renorm_4(double i0,double i1,double i2,double i3,double *o)
{double a2,a3,b1,b2,c1,d2,o0;
 QD_qsum(i2,i3,a2,a3); QD_qsum(i1,a2,b1,b2); QD_qsum(i0,b1,o0,c1);
 o[0]=o0; QD_qsum(c1,b2,o[1],d2); QD_qsum(d2,a3,o[2],o[3]);
 if (o[0]==0.0) {o[0]=o[1]; o[1]=o[2]; o[2]=o[3]; o[3]=0.0;}
 if (o[0]==0.0) {o[0]=o[1]; o[1]=o[2]; o[2]=0.0;}
 if (o[0]==0.0) {o[0]=o[1]; o[1]=0.0;}}
 
static void QD_3sum
(double a,double b,double c,double *o0,double *o1,double *o2)
{double t1,t2,t3;
 QD_sum(a,b,t1,t2); QD_sum(c,t1,*o0,t3); QD_sum(t2,t3,*o1,*o2);}

static void QD_add_22(QD a,QD b,QD c)
{double e,s,f,g;
 QD_sum(a[0],b[0],s,e); g=e+a[1]; f=g+b[1]; QD_qsum(s,f,c[0],c[1]);}
 
static void QD_sub_22(QD a,QD b,QD c)
{double e,s,f,g;
 QD_diff(a[0],b[0],s,e); g=e+a[1]; f=g-b[1]; QD_qsum(s,f,c[0],c[1]);}
 
static void QD_mul_21(QD a,double b,QD c)
{double q,r,s,t;
 QD_prod(a[0],b,s,q); t=a[1]*b; r=q+t; QD_qsum(s,r,c[0],c[1]);}
 
static void QD_div_21(QD a,double b,QD c)
{double s,t,u,v,w,r[2],x,y;
 s=a[0]/b; QD_prod(b,s,r[0],r[1]); QD_diff(a[0],r[0],t,u);
 x=u-r[1]; v=x+a[1]; y=t+v; w=y/b; QD_qsum(s,w,c[0],c[1]);}
 
static void QD_mul_22(QD a,QD b,QD c)
{double q,r,s,t1,t2,t3; QD_prod(a[0],b[0],s,q);
 t1=a[0]*b[1]; t2=a[1]*b[0]; t3=q+t1; r=t3+t2; QD_qsum(s,r,c[0],c[1]);}
 
static void QD_div_22(QD a,QD b,QD c)
{double s,t,u,v,w,r[2],x,y;
 s=a[0]/b[0]; QD_mul_21(b,s,r); QD_diff(a[0],r[0],t,u);
 x=u-r[1]; v=x+a[1]; y=t+v; w=y/b[0]; QD_qsum(s,w,c[0],c[1]);}

static void QD_sqr_2(QD i,QD o)
{double q,r,s,l;
 QD_prod(i[0],i[0],s,q); l=i[0]*i[1]; r=q+l+l; QD_qsum(s,r,o[0],o[1]);}

static void QD_add_33(QD a,QD b,QD c)
{double h,m,l,m1,m2,l1,l2;
 QD_sum(a[0],b[0],h,m1); QD_sum(a[1],b[1],m2,l1); QD_sum(m1,m2,m,l2);
 l=l1+l2+a[2]+b[2]; QD_renorm_3(h,m,l,c);}
 
static void QD_sub_33(QD a,QD b,QD c)
{double h,m,l,m1,m2,l1,l2;
 QD_diff(a[0],b[0],h,m1); QD_diff(a[1],b[1],m2,l1); QD_sum(m1,m2,m,l2);
 l=l1+l2+a[2]-b[2]; QD_renorm_3(h,m,l,c);}
 
static void QD_mul_31(QD a,double b,QD c)
{double h,m,l,m1,m2,l1,l2,t1,t2;
 QD_prod(a[0],b,h,m1); QD_prod(a[1],b,m2,l1);
 QD_sum(m1,m2,m,l2); t1=a[2]*b; t2=l1+l2; l=t1+t2; QD_renorm_3(h,m,l,c);}
 
static void QD_mul_33(QD a,QD b,QD c)
{double h,m,l,m1,m2,m3,m4,l1,l2,l3,l4,lf,ls,t1,t2;
 QD_prod(a[0],b[0],h,m1); QD_prod(a[0],b[1],m2,l1);
 QD_prod(a[1],b[0],m3,l2); t1=a[2]*b[0]; t2=a[1]*b[1]; lf=t1+t2;
 QD_sum(m1,m2,m4,l3); t1=lf+l1; t2=a[0]*b[2]; ls=t1+t2;
 QD_sum(m3,m4,m,l4); t1=ls+l2; t2=l3+l4; l=t1+t2; QD_renorm_3(h,m,l,c);}
 
static void QD_div_31(QD a,double b,QD c)
{double q0,q1,q2,t[3],r[3];
 q0=a[0]/b; QD_prod(b,q0,t[0],t[1]); t[2]=0.0; QD_sub_33(a,t,r);
 q1=r[0]/b; QD_prod(b,q1,t[0],t[1]); QD_sub_33(r,t,r);
 q2=r[0]/b; QD_renorm_3(q0,q1,q2,c);}
 
static void QD_div_33(QD a,QD b,QD c)
{double q0,q1,q2,t[3],r[3];
 q0=a[0]/b[0]; QD_mul_31(b,q0,t); QD_sub_33(a,t,r);
 q1=r[0]/b[0]; QD_mul_31(b,q1,t); QD_sub_33(r,t,r);
 q2=r[0]/b[0]; QD_renorm_3(q0,q1,q2,c);}
 
static void QD_sqr_3(QD i,QD o)
{double h,m1,m2,l1,l2,l3,m,l,t1,t2;
 QD_prod(i[0],i[0],h,m1); QD_prod(i[0],i[1],m2,l1);
 m2=m2+m2; l1=l1+l1; l2=i[2]*i[0]; t1=l2+l2; t2=i[1]*i[1]; l2=t1+t2;
 QD_sum(m1,m2,m,l3); t1=l1+l2; l=t1+l3; QD_renorm_3(h,m,l,o);}

static void QD_add_44(QD i,QD j,QD o)
{double a0,a1,b1,b2,c2,c3,d1,d2,e2,e3,f2,f3,g3,t1,t2;
 QD_sum(i[0],j[0],a0,a1); QD_sum(i[1],j[1],b1,b2); QD_sum(i[2],j[2],c2,c3);
 QD_sum(a1,b1,d1,d2); QD_sum(b2,c2,e2,e3); QD_sum(d2,e2,f2,f3);
 t1=c3+e3; t2=t1+f3; t1=i[3]+j[3]; g3=t1+t2; QD_renorm_4(a0,d1,f2,g3,o);}
 
static void QD_sub_44(QD i,QD j,QD o)
{double a0,a1,b1,b2,c2,c3,d1,d2,e2,e3,f2,f3,g3,t1,t2;
 QD_diff(i[0],j[0],a0,a1); QD_diff(i[1],j[1],b1,b2);
 QD_diff(i[2],j[2],c2,c3); QD_sum(a1,b1,d1,d2); QD_sum(b2,c2,e2,e3);
 QD_sum(d2,e2,f2,f3); t1=c3+e3; t2=t1+f3; t1=i[3]-j[3]; g3=t1+t2;
 QD_renorm_4(a0,d1,f2,g3,o);}

static void QD_mul_41(QD i,double j,QD o)
{double a0,a1,c1,c2,e2,e3,g3,h1,h2,k2,k3,l2,l3,t1;
 QD_prod(i[0],j,a0,a1); QD_prod(i[1],j,c1,c2); QD_prod(i[2],j,e2,e3);
 g3=i[3]*j; QD_sum(a1,c1,h1,h2); QD_sum(c2,e2,k2,k3);
 g3+=e3; QD_sum(h2,k2,l2,l3); t1=k3+l3; g3+=t1; QD_renorm_4(a0,h1,l2,g3,o);}
 
static void QD_mul_44(QD i,QD j,QD o)
{double a0,a1,b1,b2,c1,c2,d2,d3,e2,e3,f2,f3;
 double g3,h1,h2,h3,k2,k3,l2,l3,m2,m3,n2,n3,p2,p3;
 QD_prod(i[0],j[0],a0,a1); QD_prod(i[0],j[1],b1,b2);
 QD_prod(i[1],j[0],c1,c2); QD_prod(i[0],j[2],d2,d3);
 QD_prod(i[2],j[0],e2,e3); QD_prod(i[1],j[1],f2,f3);
 g3=i[0]*j[3]+i[3]*j[0]; QD_3sum(a1,b1,c1,&h1,&h2,&h3);
 QD_sum(b2,c2,k2,k3); QD_sum(d2,e2,l2,l3); QD_sum(f2,h2,m2,m3);
 g3+=(i[1]*j[2]+i[2]*j[1]); QD_sum(k2,l2,n2,n3); QD_sum(m2,n2,p2,p3);
 g3+=d3+e3+f3+h3+k3+l3+m3+n3+p3; QD_renorm_4(a0,h1,p2,g3,o);}
 
static void QD_div_41(QD i,double j,QD o)
{double q0,q1,q2,q3,t[4],r[4];
 q0=i[0]/j; QD_prod(j,q0,t[0],t[1]); t[2]=0.0; t[3]=0.0; QD_sub_44(i,t,r);
 q1=r[0]/j; QD_prod(j,q1,t[0],t[1]); QD_sub_44(r,t,r);
 q2=r[0]/j; QD_prod(j,q2,t[0],t[1]); QD_sub_44(r,t,r);
 q3=r[0]/j; QD_renorm_4(q0,q1,q2,q3,o);}
 
static void QD_div_44(QD i,QD j,QD o)
{double q0,q1,q2,q3,t[4],r[4];
 q0=i[0]/j[0]; QD_mul_41(j,q0,t); QD_sub_44(i,t,r);
 q1=r[0]/j[0]; QD_mul_41(j,q1,t); QD_sub_44(r,t,r);
 q2=r[0]/j[0]; QD_mul_41(j,q2,t); QD_sub_44(r,t,r);
 q3=r[0]/j[0]; QD_renorm_4(q0,q1,q2,q3,o);}

static void QD_sqr_4(QD i,QD o)
{double a0,a1,b1,b2,c2,c3,d2,d3,e3,f1,f2,g2,g3,h2,h3,k2,k3,l3,t1,t2,t3;
 QD_prod(i[0],i[0],a0,a1); QD_prod(i[0],i[1],b1,b2); b1=b1+b1; b2=b2+b2;
 QD_prod(i[0],i[2],c2,c3); c2=c2+c2; c3=c3+c3;
 QD_prod(i[1],i[1],d2,d3); t1=i[0]*i[3]; t2=i[1]*i[2]; e3=t1+t2; e3=e3+e3;
 QD_sum(a1,b1,f1,f2); QD_sum(b2,c2,g2,g3);
 QD_sum(d2,f2,h2,h3); QD_sum(g2,h2,k2,k3);
 t1=c3+d3; t2=e3+g3; t3=t2+t1; t1=h3+k3; l3=t3+t1; QD_renorm_4(a0,f1,k2,l3,o);}

double QD_2pow(int l)
{double M,C; int ll; M=1.0;
 if (l>0)
 {ll=l; C=2.0; for (l>>=1;l>0;l>>=1) {C*=C; if (l&1) M=M*C;} if (ll&1) M*=2.0;}
 else if (l<0)
 {l=-l; ll=l; C=0.5;
  for (l>>=1;l>0;l>>=1) {C*=C; if (l&1) M=M*C;} if (ll&1) M/=2.0;} return M;}

void QD_add(int n,QD a,QD b,QD c)
{switch(n)
 {case 1: {c[0]=a[0]+b[0]; return;} case 2: {QD_add_22(a,b,c); return;}
  case 3: {QD_add_33(a,b,c); return;} case 4: {QD_add_44(a,b,c); return;}}}

void QD_sub(int n,QD a,QD b,QD c)
{switch(n)
 {case 1: {c[0]=a[0]-b[0]; return;} case 2: {QD_sub_22(a,b,c); return;}
  case 3: {QD_sub_33(a,b,c); return;} case 4: {QD_sub_44(a,b,c); return;}}}

void QD_mul(int n,QD a,QD b,QD c)
{switch(n)
 {case 1: {c[0]=a[0]*b[0]; return;} case 2: {QD_mul_22(a,b,c); return;}
  case 3: {QD_mul_33(a,b,c); return;} case 4: {QD_mul_44(a,b,c); return;}}}

void QD_div(int n,QD a,QD b,QD c)
{switch(n)
 {case 1: {c[0]=a[0]/b[0]; return;} case 2: {QD_div_22(a,b,c); return;}
  case 3: {QD_div_33(a,b,c); return;} case 4: {QD_div_44(a,b,c); return;}}}

void QD_mul1(int n,QD a,double b,QD c)
{switch(n)
 {case 1: {c[0]=a[0]*b; return;} case 2: {QD_mul_21(a,b,c); return;}
  case 3: {QD_mul_31(a,b,c); return;} case 4: {QD_mul_41(a,b,c); return;}}}

void QD_div1(int n,QD a,double b,QD c)
{switch(n)
 {case 1: {c[0]=a[0]/b; return;} case 2: {QD_div_21(a,b,c); return;}
  case 3: {QD_div_31(a,b,c); return;} case 4: {QD_div_41(a,b,c); return;}}}

void QD_sqr(int n,QD i,QD o)
{switch(n)
 {case 1: {o[0]=i[0]*i[0]; return;} case 2: {QD_sqr_2(i,o); return;}
  case 3: {QD_sqr_3(i,o); return;} case 4: {QD_sqr_4(i,o); return;}}}

void QD_copy(int n,QD i,QD o) {int j; for (j=0;j<n;j++) o[j]=i[j];}
void QD_neg(int n,QD i,QD o) {int j; for (j=0;j<n;j++) o[j]=-i[j];}

void QD_mulall(int n,QD a,double b,QD c)
{int i; for(i=0;i<n;i++) c[i]=a[i]*b;}

void QD_mul2n(int n,QD a,int k,QD b) {QD_mulall(n,a,QD_2pow(k),b);}

void QD_exp(int n,QD a,QD b)
{double s[4],r[4],f,m,p[4],t[4],M,C; int i,l;
 if (a[0]==0.0) {QD_copy(n,QD_one,b); return;}
 if (a[0]<=-709.0) {QD_copy(n,QD_zero,b); return;}
 if (a[0]>=709.0) {QD_copy(n,QD_zero,b); return;}
 if (n==1) b[0]=Exp(a[0]);
 QD_div(n,a,QD_log2,t); l=(int) Round(t[0]);
 QD_mul1(n,QD_log2,(double) l,t); QD_sub(n,a,t,s);
 C=QD_2pow(2*n+2); for (i=0;i<n;i++) r[i]=s[i]/C;
 C=1.0; for (i=0;i<n;i++) C*=1.0e-16;
 QD_sqr(n,r,p); QD_add(n,r,QD_one,s);
 m=2.0; f=2.0; QD_div1(n,p,f,t);
 do {m+=1.0; QD_add(n,s,t,s); QD_mul(n,p,r,p); f*=m; QD_div1(n,p,f,t);}
 while (Abs(t[0])>C); for (i=0;i<2*n+2;i++) QD_sqr(n,s,s);
 M=QD_2pow(l); for (i=0;i<n;i++) b[i]=s[i]*M;}

void QD_log(int n,QD a,QD b)
{double q[4],r[4],s[4],t[4],u[4]; int i;
 if (a[0]<=0.0) errorit("non-positive number in QD_log");
 if (n==1) {b[0]=Log(a[0]); return;}
 QD_copy(n,QD_zero,q); q[0]=Log(a[0]);
 for (i=0;i<2;i++) r[i]=-q[i];
 QD_exp(2,r,s); QD_sub(2,q,QD_one,t); QD_mul(2,a,s,u); QD_add(2,t,u,q);
 for (i=0;i<n;i++) r[i]=-q[i];
 QD_exp(n,r,s); QD_sub(n,q,QD_one,t); QD_mul(n,a,s,u); QD_add(n,t,u,q);
 if (n<4) {QD_copy(n,q,b); return;}
 for (i=0;i<n;i++) r[i]=-q[i];
 QD_exp(n,r,s); QD_sub(n,q,QD_one,t); QD_mul(n,a,s,u); QD_add(n,t,u,b);}

static int QD_outputi(int n,int d,QD a,char *S)
{double f,*p,*b,*T,*U; int e,i,j,*Z;
 const char QD_digit[16]="0123456789abcdef";
 f=Abs(a[0]); if (f==0.0) {strcpy(S,"0.0"); return 0;}
 p=malloc(n*sizeof(double)); b=malloc(n*sizeof(double));
 T=malloc(n*sizeof(double)); U=malloc(n*sizeof(double));
 Z=malloc(sizeof(int)*(d+8));
 if (a[0]>0.0) for (i=0;i<n;i++) b[i]=a[i];
 else for (i=0;i<n;i++) b[i]=-a[i];
 T[0]=10.0; p[0]=1.0; for (i=1;i<n;i++) {p[i]=0.0; T[i]=0.0; U[i]=0.0;}
 e=(int) Floor(Log10(f));
 for (i=0;i<e;i++) QD_mul(n,p,T,p); for (i=0;i<-e;i++) QD_mul(n,p,T,p);
 if (e>0) QD_div(n,b,p,b); if (e<0) QD_mul(n,b,p,b);
 if (b[0]>10.0) {QD_div(n,b,T,b); e++;}
 else if (b[0]<1.0) {QD_mul(n,b,T,b); e--;}
 for (i=0;i<=d;i++)
 {Z[i]=(int) Floor(b[0]); U[0]=(double) Z[i];
  QD_sub(n,b,U,b); QD_mul(n,b,T,b);}
 for (i=d;i>0;i--) if (Z[i]<0) {Z[i-1]--; Z[i]+=10;}
 if (Z[d]>=5) {Z[d-1]++; i=d-1; while (i>0 && Z[i]>=10) {Z[i]-=10; Z[--i]++;}}
 i=0; if (a[0]<0.0) S[i++]='-';
 if (Z[0]==10) {S[i++]='1'; S[i++]='.'; S[i++]='0'; e++;}
 else {S[i++]=QD_digit[Z[0]]; S[i++]='.';}
 for (j=1;j<d;j++,i++) S[i]=QD_digit[Z[j]]; S[i++]='\0';
 free(p); free(b); free(U); free(T); free(Z); return e;}

void QD_output(int n,int d,QD a)
{char S[128]; int e,u; e=QD_outputi(n,d,a,S); u=strlen(S); S[u++]='E';
 if (e<0) {S[u++]='-'; e=-e;} else S[u++]='+';
 S[u++]='0'+(e/10); S[u++]='0'+(e%10); S[u]=0; printf("%s\n",S);}

void QD_ddump53(double a)
{int i,e=0,l; double y; const char QD_digit[16]="0123456789abcdef";
 if (a==0.0) {printf("0\n"); return;} if (a<0.0) {a=-a; printf("-");}
 while (a<1.0) {a*=2.0; e--;} while (a>=2.0) {a/=2.0; e++;}
 printf("1."); y=a-1.0; y=y*16.0;
 for (i=0;i<13;i++)
 {l=(int) Floor(y); printf("%c",QD_digit[l]); y=y-(double) l; y=y*16.0;}
 printf("E%i\n",e);}

void QD_floorQD(int n,QD i,QD o)
{double x0,x1=0.0,x2=0.0,x3=0.0;
 x0=Floor(i[0]); if (x0==i[0] && n>1)
 {x1=Floor(i[1]); if (x1==i[1] && n>2)
  {x2=Floor(i[2]); if (x2==i[2] && n>3) x3=Floor(i[3]);}}
 if (n==1) o[0]=x0; if (n==2) QD_qsum(x0,x1,o[0],o[1]);
 if (n==3) QD_renorm_3(x0,x1,x2,o); if (n==4) QD_renorm_4(x0,x1,x2,x3,o);} 

void QD_round(int n,QD i,QD o)
{double x0,x1=0.0,x2=0.0,x3=0.0;
 x0=Round(i[0]); if (x0==i[0] && n>1)
 {x1=Round(i[1]); if (x1==i[1] && n>2)
  {x2=Round(i[2]); if (x2==i[2] && n>3) x3=Round(i[3]);}}
 if (n==1) o[0]=x0; if (n==2) QD_qsum(x0,x1,o[0],o[1]);
 if (n==3) QD_renorm_3(x0,x1,x2,o); if (n==4) QD_renorm_4(x0,x1,x2,x3,o);}

void QD_self_renorm(int w,QD A)
{double t; if (w==2) {QD_qsum(A[0],A[1],t,A[1]); A[0]=t;}
 if (w==3) QD_renorm_3(A[0],A[1],A[2],A);
 if (w==4) QD_renorm_4(A[0],A[1],A[2],A[3],A);}

void QD_powi(int n,QD i,int p,QD o)
{QD c,t; int l=p,ll;
 switch(p)
 {case 0: QD_copy(n,QD_one,o); return; case 1: QD_copy(n,i,o); return;
  case 2: QD_sqr(n,i,o); return;
  case 3: QD_sqr(n,i,c); QD_mul(n,c,i,o); return;
  case 4: QD_sqr(n,i,c); QD_sqr(n,c,o); return;
  case 5: QD_sqr(n,i,c); QD_sqr(n,c,c); QD_mul(n,c,i,o); return;
  case 6: QD_sqr(n,i,c); QD_mul(n,c,i,c); QD_sqr(n,c,o); return;
  default:
    QD_copy(n,i,c); QD_copy(n,QD_one,t); ll=l;
    for (l>>=1;l>0;l>>=1) {QD_sqr(n,c,c); if (l&1) QD_mul(n,t,c,t);}
    if (ll&1) QD_mul(n,t,i,o); else QD_copy(n,t,o);}}

void QD_sqrt(int n,QD a,QD b)
{QD q,s,t,u,v; if (DEBUG) printf("QD_sqrt %i %f\n",n,a[0]);
 if (a[0]<0.0) errorit("negative number in QD_sqrt");
 if (a[0]<1e-60) {QD_copy(n,QD_zero,b); return;}
 if (n==1) {b[0]=Sqrt(a[0]); return;}
 QD_copy(n,QD_zero,q); q[0]=1.0/Sqrt(a[0]);
 QD_sqr(2,q,s); QD_mul(2,a,s,t); QD_sub(2,QD_one,t,u);
 QD_mul(2,q,u,v); QD_mulall(2,v,0.5,v); QD_add(2,q,v,q);
 QD_sqr(n,q,s); QD_mul(n,a,s,t); QD_sub(n,QD_one,t,u);
 QD_mul(n,q,u,v); QD_mulall(n,v,0.5,v); QD_add(n,q,v,q);
 if (n==4)
 {QD_sqr(n,q,s); QD_mul(n,a,s,t); QD_sub(n,QD_one,t,u);
  QD_mul(n,q,u,v); QD_mulall(n,v,0.5,v); QD_add(n,q,v,q);}
 QD_mul(n,q,a,b);}

void QD_cbrt(int n,QD i,QD o)
{QD q,s,t,u,v; int b=FALSE; if (DEBUG) printf("QD_cbrt %i %.16f\n",n,i[0]);
 if (i[0]<1e-60 && i[0]>-1e-60) {QD_copy(n,QD_zero,o); return;}
 if (n==1) {o[0]=Root(i[0],3); return;} if (i[0]<0.0) {QD_neg(n,i,i); b=TRUE;}
 QD_copy(n,QD_zero,q); q[0]=1.0/Root(i[0],3); /* Bad when i[0] is small! */
 QD_powi(2,q,3,s); QD_mul(2,i,s,t); QD_sub(2,QD_one,t,u);
 QD_mul(2,q,u,v); QD_div1(2,v,3.0,v); QD_add(2,q,v,q);
 QD_powi(n,q,3,s); QD_mul(n,i,s,t); QD_sub(n,QD_one,t,u);
 QD_mul(n,q,u,v); QD_div1(n,v,3.0,v); QD_add(n,q,v,q);
 if (n==4)
 {QD_powi(n,q,3,s); QD_mul(n,i,s,t); QD_sub(n,QD_one,t,u);
  QD_mul(n,q,u,v); QD_div1(n,v,3.0,v); QD_add(n,q,v,q);}
 QD_div(n,QD_one,q,o); if (b==TRUE) {QD_neg(n,i,i); QD_neg(n,o,o);}}

void QD_cos(int n,QD i,QD o)
{QD t,f,s,I,q,N,D,T; int k; if (DEBUG) printf("QD_cos %i %f\n",n,i[0]);
 if (n==1) {o[0]=Cos(i[0]); return;}
 QD_div(n,i,QD_pi,t); QD_floorQD(n,t,f); QD_mul(n,i,f,s); QD_sub(n,i,s,I);
 QD_sqr(n,I,I); QD_neg(n,I,I); QD_copy(n,I,N);
 QD_copy(n,QD_one,q); QD_copy(n,QD_one,D); D[0]=2.0;
 for (k=3;k<48;k+=2)
 {QD_div(n,N,D,T); QD_mul(n,N,I,N); QD_add(n,q,T,q);
  QD_mul1(n,D,(double) k,D); QD_mul1(n,D,(double) (k+1),D);}
 QD_copy(n,q,o);}

void QD_agm(int n,QD i1,QD i2,QD o)
{QD a,g,ap,gp; int i,D; if (DEBUG) printf("QD_agm %i %f %f\n",n,i1[0],i2[0]);
 if (i1[0]<=0.0) errorit("non-positive first input in QD_agm");
 if (i2[0]<=0.0) errorit("non-positive second input in QD_agm");
 QD_add(n,i1,i2,a); QD_mulall(n,a,0.5,a);
 QD_mul(n,i1,i2,g); QD_sqrt(n,g,g);
 while (1)
 {D=TRUE; for (i=0;i<n-1;i++) if (a[i]!=g[i]) D=FALSE; if (D) break;
  QD_copy(n,a,ap); QD_copy(n,g,gp);
  QD_add(n,ap,gp,a); QD_mulall(n,a,0.5,a); QD_mul(n,ap,gp,g); QD_sqrt(n,g,g);}
 QD_add(n,a,g,o); QD_mulall(n,o,0.5,o);}

void QD_atan(int n,QD i,QD o)
{QD q,c,c2,T,f,N,t,u,hpi; if (DEBUG) printf("QD_atan %i %f\n",n,i[0]);
 if (n==1) {o[0]=Atan(i[0]); return;}
 if (Abs(i[0])>1.0) QD_div(n,QD_one,i,u); else QD_copy(n,i,u);
 if (i[0]<0.0) QD_neg(n,u,u); QD_copy(n,QD_zero,q); q[0]=Atan(u[0]);
 QD_cos(2,q,c); QD_sqr(2,c,c2); QD_sub(2,QD_one,c2,T);
 QD_sqrt(2,T,N); QD_div(2,N,c,t); QD_sub(2,t,u,f);
 QD_mul(2,f,c2,T); QD_sub(2,q,T,q);
 QD_cos(n,q,c); QD_sqr(n,c,c2); QD_sub(n,QD_one,c2,T);
 QD_sqrt(n,T,N); QD_div(n,N,c,t); QD_sub(n,t,u,f);
 QD_mul(n,f,c2,T); QD_sub(n,q,T,q);
 if (n==4)
 {QD_cos(n,q,c); QD_sqr(n,c,c2); QD_sub(n,QD_one,c2,T);
  QD_sqrt(n,T,N); QD_div(n,N,c,t); QD_sub(n,t,u,f);
  QD_mul(n,f,c2,T); QD_sub(n,q,T,q);}
 if (Abs(i[0])>1.0) {QD_mulall(n,QD_pi,0.5,hpi); QD_sub(n,hpi,q,o);}
 else QD_copy(n,q,o); if (i[0]<0.0) QD_neg(n,o,o);}

int QD_modi(QD A,double p)
{QD T; int w=4; if (A[0]==0.0) return 0; if (Abs(A[0])<4.5e15) w=1;
 else if (Abs(A[0])<2e31) w=2; else if (Abs(A[0])<9e46) w=3;
 QD_div1(w,A,p,T); QD_floorQD(w,T,T); QD_mul1(w,T,p,T);
 QD_sub(w,A,T,T); return (int) Round(T[0]);}

llint QD_modll(QD A,double p)
{QD T; QD_div1(wmax,A,p,T); QD_floorQD(wmax,T,T); QD_mul1(wmax,T,p,T);
 QD_sub(wmax,A,T,T); return (llint) Round(T[0]);}

int QD_is_divisible(QD B,double p) {if (QD_modi(B,p)==0) return 1; return 0;}

int QD_valuation(QD A,double p)
{QD B; int v=0; if (A[0]==0.0) return (1<<28); QD_copy(wmax,A,B);
 while (QD_is_divisible(B,p)) {QD_div1(wmax,B,p,B); v++;} return v;}

void QD_mod(QD x,QD m,QD r)
{QD T; QD_div(wmax,x,m,T); QD_floorQD(wmax,T,T);
 QD_mul(wmax,T,m,T); QD_sub(wmax,x,T,T); QD_round(wmax,T,r);}

void QD_intout(QD x)
{char S[128]; int e,i=0,j; QD y,z;
 if (x[0]==0.0) {printf("+0"); return;}
 if (x[0]>1e60) {printf("+BIG"); return;}
 if (x[0]<-1e60) {printf("-BIG"); return;}
 QD_copy(wmax,QD_zero,y); y[0]=0.5; if (x[0]<0) y[0]=-0.5; QD_add(wmax,x,y,z);
 e=QD_outputi(wmax,16*wmax,z,S); if (S[0]!='-') printf("+");
 while (S[i]!='.') printf("%c",S[i++]); j=i+1;
 while (j-i<=e) printf("%c",S[j++]);}

void QD_intout_noplus(QD x)
{char S[128]; int e,i=0,j; QD y,z; if (x[0]==0.0) {printf("0"); return;}
 if (x[0]>1e60) {printf("+BIG"); return;}
 if (x[0]<-1e60) {printf("-BIG"); return;}
 QD_copy(wmax,QD_zero,y); y[0]=0.5; if (x[0]<0) y[0]=-0.5; QD_add(wmax,x,y,z);
 e=QD_outputi(wmax,16*wmax,z,S);
 while (S[i]!='.') printf("%c",S[i++]); j=i+1;
 while (j-i<=e) printf("%c",S[j++]);}

void initQDpoly(QDpoly *v,int d)
{int i; (*v).coeff=malloc((d+1)*sizeof(QD)); (*v).deg=d;
 for (i=0;i<=d;i++) QD_copy(wmax,QD_zero,(*v).coeff[i]);}
void delQDpoly(QDpoly *v) {free((*v).coeff);}

void QDpoly_add(QDpoly a,QDpoly b,QDpoly *c)
{int d,maxd,mind;
 if (a.deg>b.deg) {maxd=a.deg; mind=b.deg;} else {maxd=b.deg; mind=a.deg;}
 initQDpoly(c,maxd);
 for (d=0;d<=mind;d++) QD_add(wmax,a.coeff[d],b.coeff[d],(*c).coeff[d]);
 if (a.deg==maxd)
   for (d=mind+1;d<=maxd;d++) QD_copy(wmax,a.coeff[d],(*c).coeff[d]);
 else for (d=mind+1;d<=maxd;d++) QD_copy(wmax,b.coeff[d],(*c).coeff[d]);}

void QDpoly_mul(QDpoly a,QDpoly b,QDpoly *c,int maxd)
{int at,bt; QD T; /* if (DEBUG) printf("QDpoly_mul %i\n",maxd); */
 if (maxd==-1) maxd=a.deg+b.deg; initQDpoly(c,maxd);
 for (at=0;at<=a.deg;at++) for (bt=0;at+bt<=maxd && bt<=b.deg;bt++)
 {QD_mul(wmax,a.coeff[at],b.coeff[bt],T);
  QD_add(wmax,(*c).coeff[at+bt],T,(*c).coeff[at+bt]);}}

void QDpoly_pow(QDpoly a,int m,QDpoly *c,int maxd) /* lazy */
{int i; QDpoly s; /* if (DEBUG) printf("QDpoly_pow %i %i\n",m,maxd); */
 if (maxd==-1) maxd=a.deg*m; initQDpoly(&s,0); QD_copy(wmax,QD_one,s.coeff[0]);
 for (i=1;i<=m;i++) {QDpoly_mul(a,s,c,maxd); delQDpoly(&s); s=*c;}
 if (m==0) *c=s;}

void QDpoly_inv(QDpoly a,int m,QDpoly *c)
{int i,j; QDpoly s,t,A;
 if (m==0) {initQDpoly(c,0); QD_copy(wmax,QD_one,(*c).coeff[0]); return;}
 initQDpoly(&s,a.deg); initQDpoly(&A,0); QD_copy(wmax,QD_one,A.coeff[0]);
 for (i=1;i<=a.deg;i++) QD_neg(wmax,a.coeff[i],s.coeff[i]);
 for (j=1;j<=m;j++)
 {QDpoly_pow(s,j,&t,m); QDpoly_add(A,t,c); delQDpoly(&A); delQDpoly(&t); A=*c;}
 delQDpoly(&s);}

void QDpoly_intout(QDpoly v)
{int i;
 for (i=0;i<=v.deg;i++)
 {if (v.coeff[i][0]!=0.0)
  {if (i==0) QD_intout_noplus(v.coeff[i]); else QD_intout(v.coeff[i]);
   if (i==1) printf("*x"); if (i>=2) printf("*x^%i",i);}}
 printf("\n");}

void QDpoly_intround(QDpoly *v)
{int i; for (i=0;i<=(*v).deg;i++) QD_round(wmax,(*v).coeff[i],(*v).coeff[i]);}

static void QD_intgcdi(QD i1,QD i2,QD o)
{QD T; if (DEBUG>=2) printf("QD_intgcdi %f %f\n",i1[0],i2[0]);
 if (i1[0]==0.0) return QD_copy(wmax,i2,o);
 if (i2[0]==0.0) return QD_copy(wmax,i1,o);
 if (i1[0]>i2[0]) return QD_intgcdi(i2,i1,o);
 if (i1[0]==i2[0])
 {if (i1[1]>i2[1]) return QD_intgcdi(i2,i1,o);
  if (i1[1]==i2[1])
  {if (i1[2]>i2[2]) return QD_intgcdi(i2,i1,o);
   if (i1[2]==i2[2])
   {if (i1[3]>i2[3]) return QD_intgcdi(i2,i1,o);
    if (i1[3]==i2[3]) return QD_copy(wmax,i1,o);}}}
 QD_div(wmax,i2,i1,T); QD_floorQD(wmax,T,T);
 QD_mul(wmax,i1,T,T); QD_sub(wmax,i2,T,T); QD_intgcd(T,i1,o);}

void QD_intgcd(QD i1,QD i2,QD o)
{QD I1,I2; if (DEBUG) printf("QD_intgcd\n");
 QD_round(wmax,i1,I1); QD_round(wmax,i2,I2);
 if (I1[0]<0.0) QD_neg(wmax,I1,I1); if (I2[0]<0.0) QD_neg(wmax,I2,I2);
 QD_intgcdi(I1,I2,o);}

int QD_is_power(QD y,int m,QD z) /* only works for small-sized QD */
{int i; QD t; QD_copy(wmax,QD_zero,t); QD_copy(wmax,QD_zero,z);
 z[0]=t[0]=Round(Root(y[0],m)); QD_powi(wmax,t,m,t);
 QD_round(wmax,y,y); QD_round(wmax,t,t);
 for (i=0;i<4;i++) if (t[0]!=y[0]) return 0; return 1;}

void QD_check()
{QD x; QD_mul(wmax,QD_pi,QD_e,x);
 if (x[0]!=8.53973422267356774285) errorit("QD_check failed at x[0]");
 if (x[1]!=-6.7738152905024242803e-16) errorit("QD_check failed at x[1]");
 if (x[2]!=1.6082340642907151632e-32) errorit("QD_check failed at x[2]");}

/********************* MATH LIB REPLACEMENTS ***********************/

#define PI 3.14159265358979323846264338
#define TWOPI 6.283185307179586476925286767
#define HALFPI 1.570796326794896619231321692
#define QUARTERPI 0.7853981633974483096156608458
#define LOG2 0.6931471805599453094172321215
#define LOG10 2.302585092994045684017991455
#define SQRT2 1.414213562373095048801688724
#define SQRT2I 0.7071067811865475244008443622

double Abs(double x) {if (x<0.0) return -x; return x;}
double Round(double x) /* to even */
{if (x>=0)
 {if (x>=4503599627370496.0) return x;
  return ((x+4503599627370496.0)-(4503599627370496.0));}
 if (x<=-4503599627370496.0) return x;
  return ((x-4503599627370496.0)+(4503599627370496.0));}
double Floor(double x)
{double y; if (x<0) return -Ceil(-x);
 y=Round(x-0.5); if (x-y==1.0) return x; return y;}
double Ceil(double x)
{double y; if (x<0) return -Floor(-x);
 y=Round(x+0.5); if (y-x==1.0) return x; return y;}

double Cos(double x)
{int n=0; double s,p,b,u,E,S,d; printf("Cos %.16e\n",x);
 x=x-TWOPI*Floor(x/TWOPI); if (x>PI) x=TWOPI-x; if (x>HALFPI) {n=1; x=PI-x;}
 if (x>QUARTERPI) /* expand cos(Pi/2-x) as s-s^3/3!+s^5/120 + ... */
 {s=HALFPI-x; u=-s*s; p=s; d=1.0; b=2.0; S=s; E=Abs(s)*(1e-17);
  while ((p>E) || (p<-E)) {d*=b; b+=1.0; d*=b; b+=1.0; p*=u; S+=p/d;}
  if (n==1) return -S; return S;}
 s=x; u=-s*s; p=1.0; d=1.0; b=1.0; S=1.0; E=Abs(s)*(1e-17);
 while ((p>E) || (p<-E)) {d*=b; b+=1.0; d*=b; b+=1.0; p*=u; S+=p/d;}
 if (n==1) return -S; return S;}

double Exp(double x)
{double l,s,p,d,b,S,E;
 if (DEBUG>=2) printf("Exp %.16e\n",x); if (x==0) return 1;
 if (x>709.0) errorit("Overflow in Exp"); if (x<-709.0) return 0.0;
 l=Round(x/LOG2); x-=l*LOG2; s=x; p=1.0; d=1.0; b=1.0; S=1.0; E=Abs(s)*(1e-17);
 while ((p>E) || (p<-E)) {d*=b; b+=1.0; p*=s; S+=p/d;} return S*QD_2pow(l);}

double Log(double x)
{int e=0; double s,p,d,S,E; if (DEBUG>=2) printf("Log %.16e\n",x);
 if (x==0.0) return -1e300; if (x<=0.0) errorit("must be positive in Log");
 if (x==1.0) return 0.0;
 while (x>4.0) {x*=0.25; e+=4;} while (x>2.0) {x*=0.5; e+=2;}
 while (x>SQRT2) {x*=SQRT2I; e+=1;} while (x<0.25) {x*=4.0; e-=4;}
 while (x<0.5) {x*=2.0; e-=2;} while (x<SQRT2I) {x*=SQRT2; e-=1;}
 s=1-x; p=-1.0; d=0.0; S=0.0; E=Abs(s)*(1e-17);
 while ((p>E) || (p<-E)) {d+=1.0; p*=s; S+=p/d;} return (S+0.5*e*LOG2);}

double Log10(double x) {return Log(x)/LOG10;}

double Atan(double x)
{int r=0; double s,u,p,d,S,E; if (DEBUG>=2) printf("Atan %.16e\n",x);
 if (x==0.0) return x; if ((x>1.0) || (x<-1.0)) {x=1.0/x; r=1;}
 x=(Sqrt(4.0+4.0*x*x)-2)/(2.0*x); /* reduction */
 s=x; u=-s*s; p=s; d=1.0; S=s; E=Abs(s)*(1e-17);
 while ((p>E) || (p<-E)) {p*=u; d+=2.0; S+=p/d;} S*=2.0;
 if (r==1) {if (x>0.0) return (HALFPI-S); return (-HALFPI-S);} return S;}

double Root(double x,int n)
{double c,z; int i; if (DEBUG>=2) printf("Root %.16e %i\n",x,n);
 if (n<=0) errorit("bad Root"); z=Pow(x,1.0/(double) n);
 while (1)
 {c=1.0; for (i=0;i<n;i++) c*=z; if (c<x) z+=z/4503599627370496.0; else break;}
 while (1)
 {c=1.0; for (i=0;i<n;i++) c*=z; if (c>x) z-=z/4503599627370496.0; else break;}
 return z;}
double Sqrt(double y)
{double z,r=1.0,g=1.0,x=y; int i; if (x<0.0) errorit("Sqrt of neg");
 if (x==0.0) return 0.0;
 while (x>2.0) {x=x*0.25; r=r+r;} while (x<0.5) {x=x*4.0; r=r*0.5;}
 for (i=0;i<5;i++) g=g-(g*g-x)/(g+g); z=g*r;
 while (1) {if (z*z<y) z+=z/4503599627370496.0; else break;}
 while (1) {if (z*z>y) z-=z/4503599627370496.0; else break;} return z;}
double Pow(double x,double y)
{if (x==0.0) {if (y==0.0) return 1.0; if (y>0.0) return 0.0;}
 if (x<=0.0) errorit("must be nonnegative in Pow"); return Exp(y*Log(x));}
