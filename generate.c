#include "sympow.h"
#define DEBUG (FALSE || GLOBAL_DEBUG)

char Mtxt1[16],Stxt1[16],Mtxt2[16],Stxt2[16],Mbin1[16],Mbin2[16];
int HECKE=FALSE,sp=0,dv=0,CM=FALSE,HILO=FALSE;
int N,mx,prec; char *F1,*F2;

static void procit(char *STR)
{char TOK[8]=" "; char *INP;
 INP=strtok(STR,TOK);
 if (!strcmp(INP,"-hecke")) {HECKE=TRUE; INP=strtok(NULL,TOK);}
 if (!strcmp(INP,"-cm")) {CM=TRUE; INP=strtok(NULL,TOK);}
 if (strcmp(INP,"-sp")) errorit("No symmetric power specified");
 INP=strtok(NULL,TOK); sp=atoi(INP);
 if ((sp<=0) || (sp>25)) errorit("Symmetric power not valid");
 if ((HECKE) && (sp==1)) errorit("Hecke first power is redundant");
 if (CM && ((sp&3) | HECKE)) errorit("CM specification invalid");
 if ((sp&1) || (HECKE))
 {INP=strtok(NULL,TOK);
  if (strcmp(INP,"-dv")) errorit("No derivative specified");
  INP=strtok(NULL,TOK); dv=atoi(INP);
  if ((dv<0) || (dv>9)) errorit("Derivative not valid");}
 if (HECKE) {Mtxt1[0]=Stxt1[0]=Mtxt2[0]=Stxt2[0]=Mbin1[0]=Mbin2[0]='H';}
 else {Mtxt1[0]=Stxt1[0]=Mtxt2[0]=Stxt2[0]=Mbin1[0]=Mbin2[0]='P';}
 Mtxt1[1]=Stxt1[1]=Mtxt2[1]=Stxt2[1]=Mbin1[1]=Mbin2[1]='0'+sp/10;
 Mtxt1[2]=Stxt1[2]=Mtxt2[2]=Stxt2[2]=Mbin1[2]=Mbin2[2]='0'+sp%10;
 if ((sp&1) || (HECKE))
 {if (dv==0) Mtxt1[3]=Stxt1[3]=Mtxt2[3]=Stxt2[3]=Mbin1[3]=Mbin2[3]='E';
  if (dv==1) Mtxt1[3]=Stxt1[3]=Mtxt2[3]=Stxt2[3]=Mbin1[3]=Mbin2[3]='O';
  if (dv>1) Mtxt1[3]=Stxt1[3]=Mtxt2[3]=Stxt2[3]=Mbin1[3]=Mbin2[3]='0'+dv;}
 else
 {if (CM) {Mtxt1[3]=Stxt1[3]=Mbin1[3]='h'; Mtxt2[3]=Stxt2[3]=Mbin2[3]='l';}
  else {Mtxt1[3]=Stxt1[3]=Mbin1[3]='H'; Mtxt2[3]=Stxt2[3]=Mbin2[3]='L';}
  HILO=TRUE;}
 Mtxt1[4]=Mtxt2[4]=Mbin1[4]=Mbin2[4]='M'; Stxt1[4]=Stxt2[4]='S';
 Mtxt1[5]=Stxt1[5]=Mtxt2[5]=Stxt2[5]=Mbin1[5]=Mbin2[5]='.';
 Mtxt1[6]=Mtxt2[6]=Stxt1[6]=Stxt2[6]='t'; Mbin1[6]=Mbin2[6]='b';
 Mtxt1[7]=Mtxt2[7]=Stxt1[7]=Stxt2[7]='x'; Mbin1[7]=Mbin2[7]='i';
 Mtxt1[8]=Mtxt2[8]=Stxt1[8]=Stxt2[8]='t'; Mbin1[8]=Mbin2[8]='n';
 Mtxt1[9]=Stxt1[9]=Mtxt2[9]=Stxt2[9]=Mbin1[9]=Mbin2[9]=0;}

static void pari_params()
{int i; char F3[1024];
 int NA[14]={0,1000,600,400,320,300,300,210,200,200,200,150,150,125};
 F1=malloc(1024); F2=malloc(1024); prec=250; if (sp<14) N=NA[sp]; else N=200;
 if (HECKE)
 {N=1000; mx=dv;
  sprintf(F1,"F(k)=J(k-%i,X)*sinv(k,X)^%i/%i!",sp,dv+1,sp-1); return;}
 mx=sp/2; if (dv>mx) mx=dv; if (sp>15) prec+=(sp-15)*10;
 if (!HILO)
 {strcpy(F1,"F(k)=");
  for (i=1;i<=(1+sp)/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i-1); strcat(F1,F3);}
  sprintf(F3,"sinv(k,X)^%i",dv+1); strcat(F1,F3);}
 else if (sp&3)
 {strcpy(F1,"F(k)=if(k%2==0,");
  for (i=2;i<=1+sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i-1); strcat(F1,F3);}
  strcat(F1,"J(k/2-1,X/2)*sinv(k,X),sqrt(Pi)/2*");
  for (i=1;i<=1+sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i-1); strcat(F1,F3);}
  strcat(F1,"1/J((k-1)/2,X/2)*two1ms(k,X)*sinv(k,X))");
  strcpy(F2,"F(k)=if(k%2==1,");
  for (i=1;i<=sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i); strcat(F2,F3);}
  strcat(F2,"J((k-1)/2,X/2)*sinv(k,X),sqrt(Pi)/2*");
  for (i=1;i<=sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i); strcat(F2,F3);}
  strcat(F2,"J(k-1,X)/J(k/2-1,X/2)*two1ms(k,X)*sinv(k,X))");}
 else
 {strcpy(F1,"F(k)=if(k%2==1,");
  for (i=2;i<=1+sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i-1); strcat(F1,F3);}
  strcat(F1,"J((k-1)/2,X/2)*sinv(k,X)/sqrt(Pi),1/2*");
  for (i=1;i<=1+sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i-1); strcat(F1,F3);}
  strcat(F1,"1/J(k/2-1,X/2)*two1ms(k,X)*sinv(k,X))");
  strcpy(F2,"F(k)=if(k%2==0,");
  for (i=1;i<=sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i); strcat(F2,F3);}
  strcat(F2,"J(k/2,X/2)*sinv(k,X)/sqrt(Pi),");
  for (i=1;i<=sp/2;i++)
  {sprintf(F3,"J(k-%i,X)/%i!*",i,i); strcat(F2,F3);}
  strcat(F2,"J(k,X)/J((k-1)/2,X/2)*two1ms(k,X)*sinv(k,X))");
  if (CM) strcat(F1,"*(X-k+1)*(X-k)"); if (CM) strcat(F2,"*(X-k)*(X-k-1)");}}

int assure_line(char *STR)
{int i,j,k,l=strlen(STR); if (DEBUG) printf("assure_line %s\n",STR);
 if ((STR[0]!='P') && (STR[0]!='H') && (STR[0]!='A') 
     && (STR[0]!='M') && (STR[0]!='m')) return 0;
 if ((STR[1]>'2') || (STR[1]<'0')) return 0;
 if (!ISA_NUMBER(STR[2])) return 0;
 if ((!ISA_NUMBER(STR[3])) && (STR[3]!='E') && (STR[3]!='O')
     && (STR[3]!='h') && (STR[3]!='l') && (STR[3]!='H') && (STR[3]!='L'))
   return 0; if (STR[4]!=',') return 0;
 if (STR[l-1]!='\n') return 0;
 for (i=5;STR[i]!=',';i++) if ((STR[i]!='.') && !ISA_NUMBER(STR[i])) return 0;
 for (j=i+1;STR[j]!=',';j++)
   if ((STR[j]!='-') && !ISA_NUMBER(STR[j])) return 0;
 for (k=j+1;STR[k]!='\n';k++) if (!ISA_NUMBER(STR[k])) return 0; return 1;}

void new_sympow_s1(char *A)
{procit(A);
 printf("echo 'Removing any old data files'\n"); printf("cd datafiles\n");
 printf("%s -f %s %s %s\n",RM,Mtxt1,Stxt1,Mbin1);
 if (HILO) printf("%s -f %s %s %s\n",RM,Mtxt2,Stxt2,Mbin2); printf("cd ..\n");}

void new_sympow_pari(char *A)
{int i; procit(A); pari_params();
 printf("N=%i; dv=%i; mx=%i;\n\\p %i\n",N,dv,mx,prec);
 Stxt1[4]='\0'; printf("STR=\"%s\";\n",Stxt1); Stxt1[4]='S';
 printf("\\r standard1.gp\n%s\n\\r standard2.gp\n",F1);
 printf("\\l datafiles/%s\n\\r standard3.gp\n",Mtxt1);
 printf("\\l datafiles/%s\n",Stxt1);
 if (!HECKE) for (i=0;i<(1+sp)/2;i++) printf("coeffs(%i);\n",i);
 else printf("coeffs(0);\n");
 if (HECKE)
 {for (i=1;i<=dv;i++)
   printf("print(\"TACKS %i\");QD(polcoeff(polcoeff(P,%i,L),0,X));\n",i-1,i);
  return;}
 if (!HILO)
 {for (i=1+sp/2;i<=dv;i++)
   printf
    ("print(\"TACKS %i\");QD(polcoeff(polcoeff(P,%i,L),0,X));\n",i-1-sp/2,i);
  printf("\\q\n"); return;}
 if (sp&3) printf("coeffE(%i);\n",sp/2); else printf("coeffO(%i);\n",sp/2);
 Stxt2[4]='\0'; printf("STR=\"%s\";\n",Stxt2); Stxt2[4]='S';
 printf("\\r standard1.gp\n%s\n\\r standard2.gp\n",F2);
 printf("\\l datafiles/%s\n\\r standard3.gp\n",Mtxt2);
 printf("\\l datafiles/%s\n",Stxt2);
 for (i=0;i<(1+sp)/2;i++) printf("coeffs(%i);\n",i);
 if (sp&3) printf("coeffO(%i);\n",sp/2); else printf("coeffE(%i);\n",sp/2);
 printf("\\q\n");}

static void trimit(char *A)
{printf("%s -v '^\?' %s | %s 's/ E/e/' > .tempfile.123\\\n",GREP,A,SED);
 printf(" && echo 'END' >> .tempfile.123 && mv .tempfile.123 %s\n",A);}

void new_sympow_s2(char *A)
{procit(A);
 printf("echo 'Trimming the data files'\n"); printf("cd datafiles\n");
 trimit(Mtxt1); trimit(Stxt1); if (HILO) {trimit(Mtxt2); trimit(Stxt2);}
 printf("echo 'Turning the meshes into binaries'\n");
 printf("NUM=`%s -c AT %s`\n",GREP,Mtxt1);
 printf("../sympow -txt2bin $NUM %s < %s\n",Mbin1,Mtxt1);
 if (HILO)
 {printf("NUM=`%s -c AT %s`\n",GREP,Mtxt2);
  printf("../sympow -txt2bin $NUM %s < %s\n",Mbin2,Mtxt2);}
 printf("cd ..\n");}

void rewarp_params()
{FILE *F; char PARAM[1024][64]; int j,i=0; char a0,a1,a2,a3;
 printf("Rewarping the param_data file\n");
 F=fopen("datafiles/param_data","r");
 while (1)
 {if (!getline0(F,PARAM[i],64)) break;
  if (!assure_line(PARAM[i]))
  {printf("Found bad param_data: %s\n",PARAM[i]); continue;}
  a0=PARAM[i][0]; a1=PARAM[i][1]; a2=PARAM[i][2]; a3=PARAM[i][3];
  for (j=0;j<i;j++)
  {if ((PARAM[j][0]==a0)  && (PARAM[j][1]==a1) 
       && (PARAM[j][2]==a2) && (PARAM[j][3]==a3))
    {strcpy(PARAM[j],PARAM[i]); i--; break;}} i++;}
 fclose(F); printf("Left with %i entries in param_data\n",i);
 F=fopen("datafiles/param_data","w");
 for (j=0;j<i;j++) fprintf(F,"%s",PARAM[j]); fclose(F);}

static void read_file_mesh(int mesh_count,double *P,FILE *R)
{int mesh_number,LPT_number,QD_num; char INPUT[64],BLAH[64];
 for(mesh_number=0;mesh_number<mesh_count;mesh_number++)
 {for(LPT_number=0;LPT_number<=35;LPT_number++)
  {fscanf(R,"%s %s",INPUT,BLAH);
   for(QD_num=0;QD_num<4;QD_num++)
   {fscanf(R,"%s",INPUT);
    P[mesh_number*36*4+LPT_number*4+QD_num]=atof(INPUT);}}}}

void txt2bin(int num,char *B,FILE *R)
{double *P; FILE *F; P=malloc(36*4*num*sizeof(double));
 read_file_mesh(num,P,R); F=fopen(B,"w");
 fwrite(P,sizeof(double),36*4*num,F); fclose(F); free(P);}

#include <unistd.h>

void new_data(char *S)
{char PATH[128]="new_data",ARGS[128]="",W[128]="";
 int i=0,h=FALSE,c=FALSE,sp=0,dv=-1;
 if (ISA_NUMBER(S[i+1])) {sp=10*(S[i]-'0')+(S[i+1]-'0'); i=2;}
 else {sp=(S[i]-'0'); i=1;}
 if ((sp&1)) {ASSERT(S[i]=='d'); dv=(S[i+1]-'0'); i+=2;}
 else {if (S[i]=='d') {dv=(S[i+1]-'0'); ASSERT(S[i+2]=='h'); i+=2;}}
 if (S[i]=='h') {h=TRUE; strcpy(W,"Hecke");}
 if (S[i]=='c') {c=TRUE; strcpy(W,"CM");}
 if ((!h) && ((sp&1)==0) && (dv!=-1)) errorit("Derivative with even sympow?!");
 if ((h) && (dv==-1)) errorit("Hecke powers must have a derivative");
 printf("Make data for %s symmetric power %i",W,sp);
 if ((h) || (sp&1)) printf(" and derivative %i\n",dv); else printf("\n");
 if ((h) && (sp==1)) errorit("First Hecke power is redundant");
 if ((c) && ((sp&3)!=0)) errorit("CM power must have sp be 0 mod 4");
 if (sp>25) errorit("Symmetric power is likely too large");
 if (h) sprintf(ARGS,"-hecke -sp %i -dv %i",sp,dv);
 else if (c) sprintf(ARGS,"-cm -sp %i",sp);
 else if (sp&1) sprintf(ARGS,"-sp %i -dv %i",sp,dv);
 else sprintf(ARGS,"-sp %i",sp);
 execlp(SH,SH,PATH,SH,GP,ARGS,NULL);}
