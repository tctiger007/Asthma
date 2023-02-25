/* IL6*/
filename rf "/home/wwang750/Biofluc_testRVeffectiveness_cytokine/mrm_IL6_2.csv";
proc import datafile=rf
        out=IL6
        dbms=csv
        replace;
        run;
        
data IL6;
set IL6 ;
LINEAR = TIME;
logIL6 = log(IL6);
run;

PROC FORMAT;
VALUE GROUP 0='Healthy' 1='Asthmatic';
/* VALUE Time 1='Day 1' 2='Day 3' 3='Day 6' 4='Day 8' 5='Day 10'; */

        
* spaghetti, regression, and spline plots;
TITLE 'Observed Data, All Subjects';
PROC SGPLOT NOAUTOLEGEND DATA=IL6;
	* observed trends;
	SERIES X=TIME Y=logIL6 / GROUP=ID LINEATTRS=(THICKNESS=1);
	*overall linear ;
	REG X=TIME Y=logIL6 / NOMARKERS LINEATTRS=(COLOR=BLUE PATTERN=1 THICKNESS=3);
	*overall spline;
	PBSPLINE X=TIME Y=logIL6 / NOMARKERS LINEATTRS=(COLOR=RED THICKNESS=3);
RUN;

/* PART 1 */
PROC MEANS DATA=IL6 NWAY MEAN STD;
	CLASS TIME;
	VAR logIL6;
RUN;

	/* ***************************************************** */
/* some random effects models with autocorrelated errors */
/* ***************************************************** */

PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT /SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID R RCORR; 

PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT /SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=AR(1) R RCORR; 
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT /SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=ARH(1) R RCORR; 
		
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT /SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=ARMA(1,1) R RCORR; 
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT /SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(2) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT /SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(3) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT /SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(4) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID R RCORR; 

PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=AR(1) R RCORR; 
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=ARH(1) R RCORR; 
		
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=ARMA(1,1) R RCORR; 
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(2) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(3) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(4) R RCORR;	
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME TIME*TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID R RCORR; 

PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME TIME*TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=AR(1) R RCORR; 
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME TIME*TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=ARH(1) R RCORR; 
		
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME TIME*TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=ARMA(1,1) R RCORR; 
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME TIME*TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(2) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME TIME*TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(3) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	RANDOM INTERCEPT TIME TIME*TIME/SUB=ID TYPE=UN G GCORR;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(4) R RCORR;

/*  */
/* ******************************** */
/* some covariance structure models */
/* ******************************** */
/*  */
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=CS R RCORR; 
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=AR(1) R RCORR;

PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=ARH(1) R RCORR;

PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=ARMA(1,1) R RCORR;

PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(2) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(3) R RCORR;
	
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=TOEP(4) R RCORR;
	
	
/* Unstructured, saturated model	 */
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION;
	REPEATED LINEAR /SUB=ID TYPE=UN R RCORR;
	
	
	
	/* 	Final model	 */
PROC MIXED METHOD=ML COVTEST;
	CLASS ID LINEAR;
MODEL logIL6=TIME  GROUP GROUP*TIME  /SOLUTION outp=pred;
	REPEATED LINEAR /SUB=ID TYPE=UN R RCORR;/* 	 */
/*  */
PROC EXPORT
DATA =work.pred
DBMS=xlsx
outfile="/home/wwang750/Biofluc_testRVeffectiveness_cytokine/Output/IL6.xlsx" REPLACE; 
RUN;