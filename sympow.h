#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <strings.h>
#include "config.h"

/* flavour/type */
#ifndef FLAVOUR
#define FLAVOUR "RELEASE"
#endif

/* user level usage */
#define SYMPOW_ENV_CACHEDIR    "SYMPOW_CACHEDIR"

/* mainly for developping usage */
#define SYMPOW_ENV_PKGDATADIR  "SYMPOW_PKGDATADIR"
#define SYMPOW_ENV_PKGLIBDIR   "SYMPOW_PKGLIBDIR"
#define SYMPOW_ENV_PKGCACHEDIR "SYMPOW_PKGCACHEDIR"

/* internal usage */
#define SYMPOW_ENV_SYMPOW_PROG           "SYMPOW_INVOCATIONNAME"
#define SYMPOW_ENV_SYMPOW_OPTS_VERBOSITY "SYMPOW_OPTS_VERBOSITY"

#ifndef PREFIX
#define PREFIX "/usr/local"
#endif

#ifndef VARPREFIX
#define VARPREFIX "/var/local"
#endif

#ifndef PKGFOLDER
#define PKGFOLDER "sympow"
#endif
#ifndef PKGDATADIR
#define PKGDATADIR PREFIX "/share/" PKGFOLDER
#endif
#ifndef PKGLIBDIR
#define PKGLIBDIR PREFIX "/lib/" PKGFOLDER
#endif
#ifndef PKGCACHEDIR
#define PKGCACHEDIR VARPREFIX "/cache/" PKGFOLDER
#endif
#ifndef PKGHOMECACHEDIR
#define PKGHOMECACHEDIR "." PKGFOLDER
#endif

#ifndef ENDIANTUPLE
#define ENDIANTUPLE "bin"
#endif

#define TRUE 1
#define FALSE 0
#define VERBOSE_DEFAULT 1
#define GLOBAL_DEBUG FALSE
#define llint long long int
#define ASSERT(x) {if (!(x)) errorit("Assertion failed");}
#define wmax 4
#define CURR_MAX 60
#define MIN_TERMS 2000
#define MIN_SQRT 1000
#define STEPS 32
#define LPT 35 /* Local power series terms */
#define A0PT 20 /* about zero power series terms */
#define minimum(a,b) ((a)<(b))?(a):(b)
#define ISA_NUMBER(x) ((x>='0') && (x<='9'))
#define MASK (S_IXUSR|S_IXGRP|S_IXOTH|S_ISVTX|S_ISUID|S_ISGID)

typedef double QD[4];
extern QD QD_pi,QD_twopi,QD_sqrtpi,QD_log2,QD_one,QD_e,QD_zero;
typedef struct {int deg; QD *coeff;} QDpoly;
typedef struct {llint p[16]; int e[16];} LIST;

extern QD ***TABLE,***POWSER,*DECAY,**TACKS,*WIGGLE,*WIGSQI;
extern double *TOO_BIG,*EXPAND0_LIM,*STEP_SIZE,MD_SPEED;
extern int *evalpt,*derivative,*NUM_LOGS,*HALF_ZERO;
extern int *TACKON,SLOPPY[CURR_MAX],*MESH_COUNT;
extern QD VOLUME,SCALE,TW_EFF;
extern llint CONDTW,EVEN_TOP,EVEN_BOTTOM;
extern int MANIN,MANIN_TWIST,ZEROCHECK,RERUN,MODDEG;

extern int VERBOSE,GLOBAL,HECKE,NO_QT,TWIST,AP_SAVE,CM_CASE,CM_TWIST,ANAL_RANK;
extern int *w,*wprec; /* precision indicator */
extern llint *badprimes;
extern int *badprimetype;
extern QD COND[CURR_MAX],REAL_PERIOD,IMAG_PERIOD;
extern llint COND0;
extern QD Ea1,Ea2,Ea3,Ea4,Ea6,Eb2,Eb4,Eb6,Ec4,Ec6,Edisc,Etw4,Etw6,EtwD;
extern int C4C6LL;
extern llint Ec4ll,Ec6ll;

extern int WHICH,NUM_SUMS,fp3,fp2,*PRIMES,*RN;
extern int *whi,*wlo,*SYMPOW,*NUM_WIGS,*BLOCH_KATO,*apsave;
extern char *GET;
extern int MAX_TABLE;

extern const char *VERBOSE2option[3];

extern mode_t pkgdatamode;
extern mode_t datamode;
extern uid_t datauid;
extern gid_t datagid;
extern char *invocationname;
extern char *pkgdatadir;
extern char *pkglibdir;
extern char *pkgcachedir;
extern char *pkgdatafilesdir;
extern char *pkgdatafilesbindir;
extern char *cachedir;
extern char *datafilesdir;
extern char *datafilesbindir;
extern char *newdatascript;
extern char *paramdatafile;
extern char *bindatafiletemplate;
extern char *txtdatafiletemplate;

/* analrank.c */
void prep_analrank(llint,int);

/* analytic.c */
void get_wt_large(int,QD,QD,int,int);
void get_weight(llint,QD,int);

/* compute.c */
double go(llint,llint);
llint ec_do(llint);

/* compute2.c */
void special_value(int,QD,int);
void pth_coeff_from_ap(int,llint,int,int,QD);
void power_trace_from_ap(int,llint,int,int,QD);

/* conductors.c */
int tame_local_conductor(int,int);
int wild_local_conductor(int,int);
void badprimetype_output(int,llint);
int do_badprime(llint);
void compute_conductor(int);
void compute_conductor_hecke(int);
int get_tame_conductor(llint,int);
int get_wild_conductor(llint,int);
int get_conductor(llint,int);

/* disk.c */
void load_files(int,int,int,int);
void load_files_hecke(int,int,int,int);
int getline0(FILE*,char*,int);

/* ec_ap_bsgs.c */
llint ec_bsgs_ap(QD,QD,llint);
llint ec_bsgs_ap_AB(int,int,llint);

/* ec_ap.c */
llint ec_ap(QD,QD,llint);
llint ec_ap_with_disc(QD,QD,QD,llint);

/* ec_ap_large.c */
llint ec_bsgs_ap_large(QD,QD,llint);

/* eulerfactors.c */
void euler_factor_bad(llint,int,int,QDpoly*);
void euler_factor_hecke_bad(llint,int,int,QDpoly*);
void euler_factor(llint,int,QDpoly*);
void localinfos(char*,char*);

/* factor.c */
void init_primes();
void QD_factor(QD,LIST*);
void modular_exponentiation(llint,QD,llint,QD);
void IFACT_INIT(llint);
void ifactor(llint,LIST*,int,int);

/* fpu.s */
void fpu_53bits();

/* generate.c */
int assure_line(char*);
void new_sympow_s1(char*);
void new_sympow_pari(char*);
void new_sympow_s2(char*);
void rewarp_params();
void txt2bin(int,char*,FILE*,mode_t);
void new_data(char*);
int fork_new_data(char*);

/* help.c */
void help_message();

/* init_curve.c */
void curve_init(char*,char*);

/* moddeg.c */
void prepare_moddeg(char*);
llint postpare_moddeg();

/* periods.c */
void do_periods();

/* prepare.c */
llint prepare_decay(int,int,int);
llint process_string(char*,llint);
llint preparation(int,char*,llint);

/* QD.c */
void QD_add(int,QD,QD,QD);
void QD_sub(int,QD,QD,QD);
void QD_mul(int,QD,QD,QD);
void QD_div(int,QD,QD,QD);
void QD_mul1(int,QD,double,QD);
void QD_div1(int,QD,double,QD);
void QD_sqr(int,QD,QD);
void QD_copy(int,QD,QD);
void QD_neg(int,QD,QD);
void QD_mulall(int,QD,double,QD);
void QD_mul2n(int,QD,int,QD);
void QD_exp(int,QD,QD);
void QD_log(int,QD,QD);
void QD_round(int,QD,QD);
void QD_floorQD(int,QD,QD);
void QD_self_renorm(int,QD);
void QD_powi(int,QD,int,QD);
void QD_sqrt(int,QD,QD);
void QD_output(int,int,QD);
void QD_ddump53(double);
void errorit(char*);
void QD_agm(int,QD,QD,QD);
void QD_cos(int,QD,QD);
void QD_cbrt(int,QD,QD);
void QD_atan(int,QD,QD);
void QD_intout(QD);
void QD_intout_noplus(QD);
double QD_2pow(int);

int QD_modi(QD,double);
llint QD_modll(QD,double);
int QD_is_divisible(QD,double);
int QD_valuation(QD,double);
void QD_mod(QD,QD,QD);

void initQDpoly(QDpoly*,int);
void delQDpoly(QDpoly*);
void QDpoly_add(QDpoly,QDpoly,QDpoly*);
void QDpoly_mul(QDpoly,QDpoly,QDpoly*,int);
void QDpoly_pow(QDpoly,int,QDpoly*,int);
void QDpoly_inv(QDpoly,int,QDpoly*);
void QDpoly_intout(QDpoly);
void QDpoly_intround(QDpoly*);
void QD_intgcd(QD,QD,QD);
int QD_is_power(QD,int,QD);
void QD_check();

double Abs(double);
double Atan(double);
double Ceil(double);
double Cos(double);
double Exp(double);
double Floor(double);
double Log(double);
double Log10(double);
double Pow(double,double);
double Root(double,int);
double Round(double);
double Sqrt(double);

/* rootno.c */
int local_rootno(int,llint,int);
int global_rootno(int);

/* util.c */
void errorit(char*);
int u8(int);
void get_primes_ll(llint,llint,llint*);
void free_data();
int kron(int,int);
int kronll(llint,llint);
int gcd(int,int);
