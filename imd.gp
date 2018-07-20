\p 250
N=600; mx=1; dv=0; L; X; SP=mx+2; default(seriesprecision,SP+1); H=vector(SP);
{for(i=1,SP,H[i]=vector(N);
 for(j=1,N,if(i==1,if(j==1,H[1][1]=1,H[i][j]=H[i][j-1]+1.0/j),
           if(j==1,H[i][1]=H[i-1][1],H[i][j]=H[i][j-1]+H[i-1][j]/j*1.0))));}
Hf(n,k)=if(k==0,0,H[n][k])
ZETA(n)=if(n==1,Euler,zeta(n))
J(k,v)=if(k<0,(v-k-1)*J(k+1,v),\
                1.0*(-1)^k/k!/v*\
                exp(sum(n=1,SP,(-1)^n*(ZETA(n)/n)*v^n))*\
                (1+sum(n=1,SP,Hf(n,k)*v^n)))
sinv(k,v)=if(k==0,1/v,-1/k-sum(l=1,SP,v^l/k^(l+1)))
two1ms(k,v)=2^(1+k)*sum(l=0,SP,log(1/2)^l/l!*v^l)

{DERIV(M)= sum(i=0,poldegree(M,L),
                (deriv(polcoeff(M,i,L),X)+(i+1)*polcoeff(M,i+1,L)/X)*L^i)}
        
{ev(T,v,C)=U=0; for(i=0,mx,U+=(truncate(polcoeff(T,i,L)+O(X^C))*L^i));
 subst(subst(U,L,log(v)),X,v);}

P=sum(k=0,N,sum(l=1,mx+1,-(-1)^l/(l-1)!*polcoeff(F(k),-l,X)*L^(l-1))*X^k);
K=40; Ds=vector(K); Ds[1]=DERIV(P); for(i=2,K,Ds[i]=DERIV(Ds[i-1]));

ieee(x)=(round(x<<53)*1.0)>>53
{IEEE(x)=if(x==0,return(0));if(length(x)<2,return(0));
         y=ceil(log(abs(x))/log(2));ieee(x/2^y)*2^y;}
{QDc(x)=local(A); A=[IEEE(x),0]; A[2]=IEEE(x-A[1]); A=precision(A,18);
 print(A[1]); print(A[2]); print(0); print(0);}

{doit(x,C)=print("AT ",x); QDc(ev(P,x,C));
 for(d=1,35,print("DERIV ",d); QDc(ev(Ds[d]/d!,x,C)));}
ch(x,C)=ev(P,x,C)-ev(P,x,C+20)

setup(a,b)=for(i=0,31,doit(2^a+2^a*i/32,b));
print("About to find TOO_BIG"); \\ ends up in /dev/null
MM=100; while(log(10^(-34)+abs(ev(P,MM,N)))>-74,MM*=1.1);
write1("datafiles/param_data",STR,",",precision(IEEE(MM),18));
print("Now working backwards..."); \\ ends up in /dev/null
l=-10; while(log(10^(-34)+abs(ch(2^l,20)))<-74,l+=1); l-=1; s=l;
write1("datafiles/param_data",",",l-5); n=20;
print("Starting to write mesh files"); \\ ends up in /dev/null

dSTR=concat("datafiles/",STR); STRm=concat(dSTR,"M.txt");
default(logfile,STRm); default(log,1);
while(2^l*31/32<MM,\
      while(log(10^(-34)+abs(ch(2*2^l,n)))>-74,n+=5); setup(l,n);l+=1);
write("datafiles/param_data",",",l-s);
coeffs(j)=for(i=0,19,print("PS ",i); QDc(polcoeff(polcoeff(P,j,L),i,X)));
coO(j)=for(i=0,9,print("PSL ",2*i+1);QDc(polcoeff(polcoeff(P,j,L),2*i+1,X)));
coE(j)=for(i=0,9,print("PSL ",2*i);QDc(polcoeff(polcoeff(P,j,L),2*i,X)));

STRs=concat(dSTR,"S.txt"); default(logfile,STRs); default(log,1);
coeffs(0); if(doEV,coE(1),coO(1));
default(log,0);
