#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

#include <grp.h>
#include <pwd.h>
#include <sys/stat.h>
#include <errno.h>

#define PARAMSIZE (1024*64)

mode_t pkgdatamode=0666;
mode_t datamode=0666;
uid_t datauid=0;
gid_t datagid=0;
char *invocationname=NULL;
char *pkgdatadir=NULL;
char *pkglibdir=NULL;
char *pkgcachedir=NULL;
char *pkgdatafilesdir=NULL;
char *pkgdatafilesbindir=NULL;
char *cachedir=NULL;
char *datafilesdir=NULL;
char *datafilesbindir=NULL;
char *newdatascript=NULL;
char *paramdatafile=NULL;
char *bindatafiletemplate=NULL;
char *txtdatafiletemplate=NULL;

static const mode_t ownership2mask[4]={
 /* unique user oriented */ (S_IWGRP|S_IWOTH),
 /* one group oriented */   (S_IWOTH),
 /* everybody oriented */   00,
 /* everybody oriented */   00};

static int ownership(struct stat info)
{int status=0;
 if ((info.st_mode & S_IRWXG)==S_IRWXG)
 {gid_t gid=info.st_gid;
  if (gid==getgid()) status=1;
	{uid_t uid=getuid(); struct passwd *pwduid=getpwuid(uid);
	 struct group *grgid=getgrgid(gid);
	 char *pwdname=pwduid->pw_name; char **grmem=grgid->gr_mem;
	 while(*grmem) {if (!strcmp(pwdname,*grmem)) {status=1; break;} grmem++;}
	 }
  }
 if ((info.st_mode & S_IRWXO)==S_IRWXO) status|=2;
 return status;}

static char *dupxgetxenvdir(const char *name, size_t extra, int *flag)
{const char *ename=getenv(name); char * tname=NULL;
 if ((ename!=NULL)&&(*ename!='\0'))
 {const size_t lenghtof_ename=strlen(ename);
  char *magichead=NULL; size_t lentgthof_magichead=0; char *dum=NULL;
	tname=(char *)(malloc(lenghtof_ename+extra+1)); strcpy(tname,ename);
  magichead=tname+lenghtof_ename;
  while ((magichead!=tname)&&(*(dum=magichead-1)=='/')) {magichead=dum; ++lentgthof_magichead;}
  if (magichead!=tname) {*magichead='\0'; if (flag!=NULL) *flag=(int)(lentgthof_magichead/3);}
  else {*(tname+1)='\0'; if (flag!=NULL) *flag=0;}
  }
 else {if (flag!=NULL) *flag=0;}
 return tname;}

static const char *strlastitem(const char *string, int c) {
 const char *item=NULL;
 if (string!=NULL) {item=strrchr(string,c); if (item == NULL) item=string; else ++item;}
 return (item); }

char * dupdirname(const char *name) {
 char * foldername=NULL;
 if ((name != NULL) && (*name != '\0')) {
  const char *eofn=name+strlen(name);
  /*while ((eofn!=name)&&(*(--eofn)=='/')) ;*/
  while ((eofn!=name)&&(*eofn!='/')) --eofn;
  if (eofn!=name) {
   const size_t ublofn=eofn-name;
   size_t cdx=0; const char * dim=NULL; char * dam=NULL; int flag=0;
   foldername=(char *)(malloc(ublofn)); memset(foldername,'\0',ublofn);
   for(cdx=0,dim=name,dam=foldername;cdx<ublofn;++dim,++cdx) {
    if (*dim == '/') flag=1;
    else {if (flag) {*dam='/'; ++dam; flag=0;} *dam=*dim; ++dam;}}
   if (flag) *dam='/';
   }
  else {foldername=strdup((*name != '/')?".":"/");}}
	return foldername;}

static void prepare_main(char *argv0)
{char *env_pkgdir; struct stat infod, infodb, infof; int flag;
 invocationname=argv0;
 env_pkgdir=getenv(SYMPOW_ENV_PKGDATADIR); pkgdatadir=(env_pkgdir!=NULL)?env_pkgdir:PKGDATADIR;
 env_pkgdir=getenv(SYMPOW_ENV_PKGLIBDIR); pkglibdir=(env_pkgdir!=NULL)?env_pkgdir:PKGLIBDIR;
 env_pkgdir=getenv(SYMPOW_ENV_PKGCACHEDIR); pkgcachedir=(env_pkgdir!=NULL)?env_pkgdir:PKGCACHEDIR;
 datauid=getuid();
 if ((cachedir=dupxgetxenvdir(SYMPOW_ENV_CACHEDIR,sizeof(PKGFOLDER),&flag))!=NULL) {
  if (strncmp(strlastitem(cachedir,'/'),PKGFOLDER,sizeof(PKGFOLDER)-1)) {strcat(cachedir,"/"PKGFOLDER);};
  if ((flag)&&(stat(cachedir,&infod))) {char *parcachedir=dupdirname(cachedir);
   if (stat(parcachedir,&infod)) {fprintf(stderr,"**ERROR** parent cache directory %s does not exist\n",parcachedir); exit(-1);}
   if (!S_ISDIR(infod.st_mode)) {fprintf(stderr,"**ERROR** parent cache directory %s exists but is not a directory\n",parcachedir); exit(-1);}
	 if (access(parcachedir,(R_OK|W_OK|X_OK))) {fprintf(stderr,"**ERROR** parent cache directory %s yields insufficient permissions\n",parcachedir); exit(-1);}
	 if (1<flag) {int os=ownership(infod); mode_t mask=umask(ownership2mask[os]);
    if (mkdir(cachedir,infod.st_mode)) {fprintf(stderr,"**ERROR** failed to create cache folder %s\n",cachedir); exit(-1);}
    if (os&1) {chown(cachedir,datauid,infod.st_gid);}
    umask(mask);}
   else {if (mkdir(cachedir,(S_IRWXU|S_IRGRP|S_IXGRP))) {fprintf(stderr,"**ERROR** failed to create cache folder %s\n",cachedir); exit(-1);}}
   free(parcachedir);}
  }
 else if ((env_pkgdir=getenv("HOME"))!=NULL) {asprintf(&cachedir,"%s/"PKGHOMECACHEDIR,env_pkgdir);
  if (stat(cachedir,&infod)) {if (mkdir(cachedir,(S_IRWXU|S_IRGRP|S_IXGRP)))
	 {fprintf(stderr,"**ERROR** failed to create local cache folder %s\n",cachedir); exit(-1);}}
 	}
 else {fprintf(stderr,"**ERROR** unset environment variable HOME\n"); exit(-1);}
 asprintf(&pkgdatafilesdir,"%s/datafiles",pkgdatadir);
 if (stat(pkgdatafilesdir,&infod)) {free(pkgdatafilesdir); pkgdatafilesdir=NULL;}
 asprintf(&pkgdatafilesbindir,"%s/datafiles/"ENDIANTUPLE,pkgcachedir);
 if (stat(pkgdatafilesbindir,&infodb)) {mode_t mask=umask(0);
  if (mkdir(pkgdatafilesbindir,(S_IRWXU|S_IRWXG|S_IRWXO|S_ISVTX)))
	{if (VERBOSE>=1) fprintf(stderr,"**WARNING** failed to create data bin package cache folder %s\n",pkgdatafilesbindir);
   free(pkgdatafilesbindir); pkgdatafilesbindir=NULL;}
  else
  {stat(pkgdatafilesbindir,&infodb); pkgdatamode= infodb.st_mode & ~MASK;}
	umask(mask);}
 else
 {if (!S_ISDIR(infodb.st_mode))
  {if (VERBOSE>=1) fprintf(stderr,"**WARNING** %s exists but is not a directory\n",pkgdatafilesbindir);
   free(pkgdatafilesbindir); pkgdatafilesbindir=NULL;}
  else if (access(pkgdatafilesbindir,(R_OK|W_OK|X_OK)))
  {if (VERBOSE>=1) fprintf(stderr,"**WARNING** %s yields insufficient permissions\n",pkgdatafilesbindir);
   free(pkgdatafilesbindir); pkgdatafilesbindir=NULL;}
	else {pkgdatamode= infodb.st_mode & ~MASK;}}
 asprintf(&datafilesdir,"%s/datafiles",cachedir);
 asprintf(&datafilesbindir,"%s/"ENDIANTUPLE,datafilesdir);
 asprintf(&newdatascript,"%s/new_data",pkglibdir);
 asprintf(&paramdatafile,"%s/param_data",datafilesdir);
 if (stat(datafilesdir,&infod))
 {if (stat(cachedir,&infod)) {fprintf(stderr,"**ERROR** cache directory %s does not exist\n",cachedir); exit(-1);}
  if (!S_ISDIR(infod.st_mode)) {fprintf(stderr,"**ERROR** cache directory %s exists but is not a directory\n",cachedir); exit(-1);}
	if (access(cachedir,(R_OK|W_OK|X_OK))) {fprintf(stderr,"**ERROR** cache directory %s yields insufficient permissions\n",cachedir); exit(-1);}
	int os=ownership(infod); umask(ownership2mask[os]);
	if (mkdir(datafilesdir,infod.st_mode)) {fprintf(stderr,"**ERROR** failed to create data cache folder %s\n",datafilesdir); exit(-1);}
	if (os&1) {chown(datafilesdir,datauid,infod.st_gid);}
	stat(datafilesdir,&infod);
  }
 else
 {if (!S_ISDIR(infod.st_mode)) {fprintf(stderr,"**ERROR** data cache folder %s exists but is not a directory\n",datafilesdir); exit(-1);}
  if (access(datafilesdir,(R_OK|W_OK|X_OK))) {fprintf(stderr,"**ERROR** data cache folder %s yields insufficient permissions\n",datafilesdir); exit(-1);}}
 if (stat(datafilesbindir,&infodb))
 {if (mkdir(datafilesbindir,infod.st_mode)) {fprintf(stderr,"**ERROR** failed to create data bin cache folder %s\n",datafilesbindir); exit(-1);}}
 else
 {if (!S_ISDIR(infodb.st_mode)) {fprintf(stderr,"**ERROR** data bin cache folder %s exists but is not a directory\n",datafilesbindir); exit(-1);}
  if (access(datafilesbindir,(R_OK|W_OK|X_OK))) {fprintf(stderr,"**ERROR** data bin cache folder %s yields insufficient permissions\n",datafilesbindir); exit(-1);}}
 umask(MASK);
 if (stat(paramdatafile,&infof))
 {int fd=-1; datamode= infod.st_mode & ~(S_ISVTX);
 	fd=open(paramdatafile,(O_WRONLY|O_CREAT|O_TRUNC),datamode);
	if (fd==-1) {fprintf(stderr,"**ERROR** failed to create data file %s\n",paramdatafile); exit(-1);}
	if (pkgdatafilesdir) {int df=-1; char *pkgparamdatafile=NULL;
   asprintf(&pkgparamdatafile,"%s/param_data",pkgdatafilesdir); df=open(pkgparamdatafile,O_RDONLY);
	 if (-1<df) {char PARAM[PARAMSIZE]; ssize_t nread=0;
	  while (0<(nread=read(df,PARAM,PARAMSIZE)))
		{if ((write(fd,PARAM,(size_t)(nread)))!=nread)
		 {fprintf(stderr,"**ERROR** failed to copy data file %s to %s\n",pkgparamdatafile,paramdatafile); exit(-1);}}
    if (nread<0) {fprintf(stderr,"**ERROR** failed to read data file %s\n",pkgparamdatafile); exit(-1);}
    close(df);}
   else if (errno!=ENOENT) {fprintf(stderr,"**ERROR** failed to open data file %s\n",pkgparamdatafile); exit(-1);}
	 free(pkgparamdatafile);
	 }
	close(fd);
	int os=ownership(infod);
	if (os&1) {chown(paramdatafile,datauid,infod.st_gid);}
	stat(paramdatafile,&infof);
  }
 else
 {if (!S_ISREG(infof.st_mode)) {fprintf(stderr,"**ERROR** %s exists but is not a regular file\n",paramdatafile); exit(-1);}
  if (access(paramdatafile,(R_OK|W_OK))) {fprintf(stderr,"**ERROR** %s yields insufficient permissions\n",paramdatafile); exit(-1);}}
 datamode= infof.st_mode & ~MASK;
 datagid=infof.st_gid;
 if (pkgdatafilesbindir==NULL) {asprintf(&pkgdatafilesbindir,"%s",datafilesbindir); pkgdatamode=datamode;}
 asprintf(&txtdatafiletemplate,"%s/PXXXXXXXXxxxxxxxxXXXXXXXXxxxxxxxx",
  ((pkgdatafilesdir) && (strlen(datafilesdir) < strlen(pkgdatafilesdir)))?pkgdatafilesdir:datafilesdir);
 asprintf(&bindatafiletemplate,"%s/PXXXXXXXXxxxxxxxxXXXXXXXXxxxxxxxx",
 	(strlen(datafilesbindir) < strlen(pkgdatafilesbindir))?pkgdatafilesbindir:datafilesbindir);
 }

static void postpare_main()
{free(txtdatafiletemplate); free(bindatafiletemplate);
 free(paramdatafile); free(newdatascript);
 free(datafilesbindir); free(datafilesdir); free(cachedir);
 free(pkgdatafilesbindir);
 if (pkgdatafilesdir) free(pkgdatafilesdir);
 /* free(pkgcachedir); free(pkglibdir); free(pkgdatadir); */
 /*free(invocationname);*/}

int main(int argc,char **argv)
{char INSTRING[1024]="2w3s1p32,3bp16d1,4p8\0",CSTR[1024]="\0",LSTR[1024]="\0";
 char TYPE[128]; int i=1; llint NT,UB=(((llint) 1)<<45);
 int NO_CM=FALSE,info=0,ROOTNO=FALSE,SLOPPY=0,QD_CHECK;
 NO_QT=FALSE; VERBOSE=VERBOSE_DEFAULT; GLOBAL=TRUE; HECKE=FALSE; TWIST=FALSE; AP_SAVE=0;
 CM_CASE=FALSE; GET=malloc(1024); COND0=1; fp3=0; fp2=0; MAX_TABLE=1<<27;
 MODDEG=FALSE; ANAL_RANK=FALSE; ZEROCHECK=FALSE; RERUN=FALSE; QD_CHECK=TRUE;
#define PROCNAME_PREPOSTPARE(A) { prepare_main(*argv); A; postpare_main(); }
#if defined(ISOC99_FENV) || defined(FPUCONTROLH) || defined(x86)
 fpu_53bits();
#endif
 snprintf(TYPE,sizeof(TYPE),"%s",FLAVOUR); MD_SPEED=2.0;
 while(i<argc)
 {if (!strcmp(argv[i],"-quiet")) {VERBOSE=0; i++;}
  else if ((!strcmp(argv[i],"-sympow")) || (!strcmp(argv[i],"-sp")))
  {strcpy(INSTRING,argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-terse")) {VERBOSE=1; i++;}
  else if (!strcmp(argv[i],"-verbose")) {VERBOSE=2; i++;}
  else if (!strcmp(argv[i],"-help")) help_message();
  else if (!strcmp(argv[i],"-version")) { printf("%s %s\n",VERSION,TYPE); exit(0);}
  else if (!strcmp(argv[i],"-dump-endiantuple")) { printf("%s\n",ENDIANTUPLE); exit(0);}
  else if (!strcmp(argv[i],"-dump-versiontuple")) { printf("%s\n",VERSION); exit(0);}
  else if (!strcmp(argv[i],"-curve")) {strcpy(CSTR,argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-label")) {strcpy(LSTR,argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-info")) {info=i+1; i+=3;}
  else if (!strcmp(argv[i],"-bound")) {UB=atoll(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-maxtable")) {MAX_TABLE=atoll(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-hecke")) {HECKE=TRUE; i++;}
  else if (!strcmp(argv[i],"-mdspeed")) {MD_SPEED=atof(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-nocm")) {NO_CM=TRUE; i++;}
  else if (!strcmp(argv[i],"-noqdcheck")) {QD_CHECK=FALSE; i++;}
  else if (!strcmp(argv[i],"-noqt")) {NO_QT=TRUE; i++;}
  else if (!strcmp(argv[i],"-local")) {GLOBAL=FALSE; i++;}
  else if (!strcmp(argv[i],"-moddeg")) {MODDEG=TRUE; i++;}
  else if (!strcmp(argv[i],"-analrank")) {ANAL_RANK=TRUE; i++;}
  else if (!strcmp(argv[i],"-sloppy")) {SLOPPY=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-rootno")) {ROOTNO=atoi(argv[i+1]); i+=2;}
  else if (!strcmp(argv[i],"-rewarp")) { PROCNAME_PREPOSTPARE(rewarp_params()) return 0;}
#ifdef NEW_DATA
  else if (!strcmp(argv[i],"-shell1")) { PROCNAME_PREPOSTPARE(new_sympow_s1(argv[i+1])) return 0;}
  else if (!strcmp(argv[i],"-pari")) { PROCNAME_PREPOSTPARE(new_sympow_pari(argv[i+1])) return 0;}
  else if (!strcmp(argv[i],"-shell2")) { PROCNAME_PREPOSTPARE(new_sympow_s2(argv[i+1])) return 0;}
  else if (!strcmp(argv[i],"-txt2bin"))
  { PROCNAME_PREPOSTPARE(txt2bin(atoi(argv[i+1]),argv[i+2],stdin,datamode)) return 0;}
  else if (!strcmp(argv[i],"-new_data")) { PROCNAME_PREPOSTPARE(new_data(argv[i+1])) return 0;}
#else
  else if ((!strcmp(argv[i],"-shell1")) || (!strcmp(argv[i],"-pari")) ||
	   (!strcmp(argv[i],"-shell2")) || (!strcmp(argv[i],"-txt2bin")) ||
	   (!strcmp(argv[i],"-new_data")))
    errorit("new_data not possible --- try re-configuring SYMPOW");
#endif
  else
	 {printf("sympow: unrecognised command\nTry sympow -help for more information\n"); exit(-1);}}
 if (VERBOSE) printf("sympow %s %s  (c) Mark Watkins -",VERSION,TYPE);
 if (QD_CHECK) QD_check();
 if (VERBOSE) {printf("-- see README and COPYING for details\n");}
 prepare_main(*argv);
 init_primes(); if (!CSTR[0]) getline0(stdin,CSTR,256); curve_init(CSTR,LSTR);
 if (HECKE && !CM_CASE) errorit("Curve does not have CM");
 if (NO_CM && CM_CASE) errorit("Curve has CM");
 if (VERBOSE && CM_CASE)
   printf("Curve has complex multiplication %i\n",CM_CASE);
 if (info) {localinfos(argv[info],argv[info+1]); return 0;}
 if (ROOTNO) {printf("Root number for sp:%i is %i\n",
		     ROOTNO,global_rootno(ROOTNO)); return 0;}
 if (MODDEG) {HECKE=FALSE; prepare_moddeg(INSTRING);}
 if (ANAL_RANK) {HECKE=FALSE; prep_analrank(UB,SLOPPY); return 0;}
 if (VERBOSE && HECKE) printf("Using Hecke powers\n");
 NT=preparation(CURR_MAX,INSTRING,UB);
 if (VERBOSE>=2) printf("NUM_SUMS is %i\n",NUM_SUMS);
 if (NUM_SUMS>CURR_MAX) errorit("Too many sums!");
 if (NUM_SUMS==0) printf("Nothing to do!\n");
 fflush(stdout); if (GLOBAL && NUM_SUMS) go(NT,NT);
 if (ZEROCHECK)
 {NUM_SUMS=0; WHICH=0; free_data(); RERUN=TRUE;
  NT=process_string(INSTRING,UB); go(NT,NT);}
 postpare_main();
#undef PROCNAME_PREPOSTPARE
 return 0;}
