/**********************************************************************************************************************
* Author: Alex Keil
* Project: INWORKS (International Nuclear Workers Study)
* Description: Creating survival curves from g-formula simulation output using Cox models and baseline hazard estimators
* Keywords: radiation, occupational, g-formula
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/


%MACRO mksurvdata (ds=nc, invar=inage, outvar=age, suff=, intervention=,breakties=1);
  %IF &dobootstrap=1 %THEN %DO;
     %LET dsname= cidata_&outvar._&ds&suff._&bootiter;
  %END;
  %IF &dobootstrap=0 %THEN %DO;
     %LET dsname= cidata_&outvar._&ds&suff;
  %END;

 *PHREG cannot handle negative times, so convert sas dates to number of days from the minimum start time (will be converted back at end);
 %LET mintime = -900000; *cause an obvious problem in the graph if assignment of this does not work;
 TITLE "gformula.&dsname";
 DATA ___sur;
  SET &DS;
  %IF "&intervention" NE "" %THEN IF intervention=&intervention;;
  %LET mintime = 0; * turning this off for age time scale;
 DATA ___sur;
  SET ___sur;
  &invar =  &invar - &mintime;
  &outvar = &outvar - &mintime;
  d_multi = d_&CAUSE1 + 2*d_&CAUSE2;
  *random breaking of ties for mortality outcomes only (introduced 3/11/2017);
  %IF &breakties=1 %THEN %DO;
  CALL STREAMINIT(12323100);
  IF d_multi>0 THEN DO;
   ct=1;
   __tmp = &outvar;
   DO WHILE (__tmp <= &invar OR ct=1) ;
    __tmp = &outvar - (RAND('uniform'))*0.0001;
    ct=ct+1; 
   END;
   &outvar=__tmp;
   DROP ct;
  END;
  %END;
 RUN;
 *causes of death;
 PROC PHREG DATA = ___sur NOPRINT;
  MODEL (&invar &outvar)*d_multi(0) =  / TIES=EFRON;
  OUTPUT OUT=surv_multi_&outvar._&ds(KEEP=d_multi &invar &outvar surv_y atrisk) 
    SURVIVAL=surv_y ATRISK=atrisk / METHOD=BRESLOW ;
 RUN;
 /* programming note: this restricts the survival dataset to the observations where there are jumps*/
 PROC SORT DATA = surv_multi_&outvar._&ds(WHERE=(d_multi>0)); BY &outvar;RUN;
 
 DATA &dsname;
  SET surv_multi_&outvar._&ds(KEEP=&outvar surv_y d_multi atrisk);
    &outvar = &outvar + &mintime;
 RUN;
 PROC SORT DATA = &dsname; BY &outvar DESCENDING surv_y;
 
 DATA &dsname
 (KEEP = &outvar ci_: atrisk) 
 ;
  SET &dsname;
  BY &outvar DESCENDING surv_y;
  RETAIN lastsurv_y 1 ci_allcause ci_&CAUSE1  ci_&CAUSE2 0;
   IF FIRST.&outvar THEN DO; 
    IF surv_y<=.z THEN surv_y=lastsurv_y;
    ci_allcause=1-surv_y;
  END;
   ci_&CAUSE1  + MAX((lastsurv_y-surv_y)*(d_multi=1), 0);
   ci_&CAUSE2  + MAX((lastsurv_y-surv_y)*(d_multi=2), 0);
   *ci_lw  = MAX(0, 1-surv_l, ci_lw);
   altci_allcause = SUM(ci_&CAUSE1, ci_&CAUSE2);
   IF surv_y<=.z THEN surv_y = 1-ci_allcause;
   IF LAST.&OUTVAR THEN OUTPUT;
   lastsurv_y=surv_y;
  RUN;
  
  PROC SORT DATA = &dsname; BY DESCENDING &outvar;
  
  /* programming note: this takes the last observation for a year/5year/month period and assigns survival/ci to that rounded time */
  DATA &dsname(KEEP = time ci_: atrisk mark);
    SET &dsname;
    BY DESCENDING &outvar; 
    RETAIN nexttime 0;
   *OUTPUT test; 
    IF LAST.&OUTVAR THEN DO;
      IF  FLOOR(&outvar / 5) NE FLOOR(nexttime / 5) THEN mark=5;
      ELSE IF  FLOOR(&outvar ) NE FLOOR(nexttime ) THEN mark=1;
      ELSE IF  FLOOR(&outvar *12) NE FLOOR(nexttime *12) THEN mark=0.083;
      ELSE mark = 0;
      time = ROUND(ROUND(&outvar*12)/12, .001); 
      OUTPUT;
    END;
    nexttime = &outvar;
 RUN;
 PROC SORT DATA = &dsname(WHERE=(mark>0)); BY time ci_allcause; RUN;
 * remove ties for compact output dataset;
 DATA &dsname;
  SET &dsname;
  BY time ci_allcause;
  IF last.time THEN OUTPUT;
 RUN;
 
 PROC PRINT DATA= &dsname(WHERE=(mark>1));
 RUN;
 *export as CSV;
 DATA __xpt; SET &dsname(RENAME=(time=&outvar));
 PROC EXPORT DATA = __xpt(DROP=mark) DBMS = CSV OUTFILE = "&outdir./&dsname..csv" REPLACE; RUN;
 PROC SQL; DROP TABLE __xpt; QUIT;
%MEND;

*step 0 - data preparation;
%MACRO reduceobs(invar=inage, outvar=outage, newinvar=agest, newoutvar=ageend);
  PROC SQL;
   CREATE TABLE obs AS SELECT 
     wkerid, cohort, d_&CAUSE1, d_&CAUSE2, 1 AS sampwt,
     MIN(&invar) AS &newinvar, MAX(&outvar) AS &newoutvar
    FROM an(WHERE=(inage>=16)) /*added 2/22*/
    GROUP BY cohort, wkerid
    HAVING &outvar = &newoutvar
    ORDER BY cohort, wkerid;
  QUIT;
%MEND;


%MACRO procobs();
  %IF &dobootstrap=0 %THEN %DO;
   %REDUCEobs(invar=inage, outvar=outage, newinvar=agest, newoutvar=ageend);
   DATA obsfr;
    SET obs (WHERE=(cohort IN(5)));
   DATA obsuk;
    SET obs (WHERE=(cohort IN(13)));
   DATA obsus;
    SET obs (WHERE=(cohort IN(14)));
   RUN;
   
   *observed data;
   %mksurvdata(DS=obsuk, invar=agest, outvar=ageend, suff=, breakties=1);
   %mksurvdata(DS=obsfr, invar=agest, outvar=ageend, suff=, breakties=1);
   %mksurvdata(DS=obsus, invar=agest, outvar=ageend, suff=, breakties=1);
   %mksurvdata(DS=obs,   invar=agest, outvar=ageend, suff=, breakties=1);
   
   PROC DATASETS LIBRARY=WORK; DELETE obs obsuk obsus obsfr;QUIT;
  %END;
%MEND;
  
%MACRO reducegf_old(invar=inage, outvar=outage, newinvar=agest, newoutvar=ageend);
  *make a smaller dataset, assuming gf uses full person time data as output;
  PROC SQL;
   CREATE TABLE gf AS SELECT 
     intervention, gfid, cohort, d_&CAUSE1, d_&CAUSE2, censadmin, 1 AS sampwt,
     MIN(&invar) AS &newinvar, MAX(&outvar) AS &newoutvar
    FROM gformula2
    GROUP BY gfid /* note gfid must be unique for each simulated individual - no duplicates across interventions */
    HAVING &outvar = &newoutvar
    ORDER BY intervention, cohort, gfid;
  QUIT;
%MEND;

%MACRO reducegf(invar=inage, outvar=outage, newinvar=agest, newoutvar=ageend);
  *make a smaller dataset assuming gf only outputs a single observation per intervention/pseudo-person combination;
  PROC SQL;
   CREATE TABLE gf AS SELECT 
     intervention, gfid, cohort, d_&CAUSE1, d_&CAUSE2, censadmin, 1 AS sampwt,
     &invar AS &newinvar, &outvar AS &newoutvar
    FROM gformula2
    ORDER BY intervention, cohort, gfid;
  QUIT;
%MEND;


%PROCOBS();
%REDUCEgf(invar=inage, outvar=outage, newinvar=agest, newoutvar=ageend);

PROC PRINT DATA = gf (obs=10);
 TITLE "Data used for ci calculation";
RUN;

*g formula;
DATA gffr;
 SET gf (WHERE=(cohort IN(5)));
DATA gfuk;
 SET gf (WHERE=(cohort IN(13)));
DATA gfus;
 SET gf (WHERE=(cohort IN(14)));
RUN;

* create survival curve data (risk by age);
%MACRO allcurves(invar,outvar);
   %LET __inti = 1;
   %DO %WHILE("%SCAN(&INTS, &__inti, %STR(,))"^="");
     %LET Jint = %QUOTE(%SCAN(&INTS, &__inti, %STR(,)));
	 %LET int = %sysfunc(tranwrd(%sysfunc(tranwrd(&JINT,.,p)),-,m));
      %mksurvdata(DS=gfuk, invar=&invar, outvar=&outvar, suff=&int, intervention=&jint);
      %mksurvdata(DS=gffr, invar=&invar, outvar=&outvar, suff=&int, intervention=&jint);
      %mksurvdata(DS=gfus, invar=&invar, outvar=&outvar, suff=&int, intervention=&jint);
      %mksurvdata(DS=gf,   invar=&invar, outvar=&outvar, suff=&int, intervention=&jint);
  %LET __inti = %EVAL(&__inti + 1);
  %END;
%MEND;

%allcurves(invar=agest, outvar=ageend);


RUN;QUIT;RUN;

