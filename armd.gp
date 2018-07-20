
STR="M02H"; doEV=1;
F(k)=if(k%2==0,J(k-2,X)/1!*J(k/2-1,X/2)*sinv(k,X),\
    sqrt(Pi)/2*J(k-1,X)/0!*J(k-2,X)/1!*1/J((k-1)/2,X/2)*two1ms(k,X)*sinv(k,X))
\r imd.gp
STR="M02L"; doEV=0;
F(k)=if(k%2==1,J(k-1,X)/1!*J((k-1)/2,X/2)*sinv(k,X),\
    sqrt(Pi)/2*J(k-1,X)/1!*J(k-1,X)/J(k/2-1,X/2)*two1ms(k,X)*sinv(k,X))
\r imd.gp
STR="m01E"; dv=0; sh=1;
\r ihmd.gp
STR="m02E"; dv=0; sh=2;
\r ihmd.gp
STR="A01O"; dv=1;
\r ianalrank.gp
STR="A012"; dv=2;
\r ianalrank.gp
STR="A013"; dv=3;
\r ianalrank.gp
STR="A014"; dv=4;
\r ianalrank.gp
STR="A015"; dv=5;
\r ianalrank.gp
STR="A016"; dv=6;
\r ianalrank.gp
STR="A017"; dv=7;
\r ianalrank.gp
STR="A018"; dv=8;
\r ianalrank.gp
STR="A019"; dv=9;
\r ianalrank.gp

\q
