/**********************************************************************************************************************
* Author: Alex Keil
* Project: INWORKS (International Nuclear Workers cohort)
* Description: G-formula macro for INWORKS data (model, simulate data under interventions)
* Keywords: radiation, occupational, g-formula
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/


PROC FORMAT;
 VALUE cohort
  14="US"
  13="UK"
   5="France"
  .="Overall";
 VALUE intform
  .1  = 'b) <0.1 mSv'
  .5  = 'c) <0.5 mSv'
  .75 = 'd) <0.75 mSv'
  1   = 'e) <1 mSv'
  2   = 'f) <2 mSv'
  3   = 'g) <3 mSv'
  5   = 'h) <5 mSv'
  10   = 'i) <10 mSv'
  -1   = 'j) =1 mSv'
  -2   = 'k) =2 mSv'
  -3   = 'l) =3 mSv'
  -10   = 'm) =10 mSv'
  -20   = 'n) =20 mSv'
  -999  = 'o) Never exposed'
  ;

********* BEGIN PROGRAMMING STATEMENTS ************************;
DATA an;
  SET anpy;
  BY wkerid indate;
  %TFTRANSFORM();
  %TVtimeTRANSFORM();
  %TVworkTRANSFORM();
  %TVdoseTRANSFORM();
  %TVlagworkTRANSFORM();
  %TVlagdoseTRANSFORM();
  female = sex - 1;
PROC SORT DATA = an;
  BY cohort wkerid inage;
DATA anfirst;
  SET an;
  BY cohort wkerid inage;
  IF first.wkerid THEN OUTPUT;
RUN;



/************************************************************************;
*step 1 - modeling;

* modeling effects of dose: be flexible with coding to allow different elements of flexibility
* think about attained age: bins of age at exposure: 9 bins -> plan to this in simulation
*  age at exposure X time since exposure (bin those (<30, 30-45, 45+ age; published?))
*
* From Richardson et al (2018)
*  an estimate of the coefficient of primary interest, Bj, was adjusted to account for the 
*  effects of country, attained age (in 5-year intervals), sex, year of birth (in 10-year 
*  intervals), socioeconomic status (in five categories, based on job title, for French, 
*  US, and UK workers employed by the Atomic Energy Authority and Atomic Weapons Establishment;
*  other UK workers were classified as nonmanual or manual skilled workers, based on employment 
*  category), duration of employment or radiation work (in 10-year intervals), and exposure to 
*  neutrons (whether a worker had a positive recorded neutron dose, and if so, whether their 
*  recorded neutron dose equivalent ever exceeded 10% of their total external radiation dose 
*  equivalent
************************************************************************/
/* define predictors */

*leaving employment;
%LET lpreds =  usa uk
			   frn*female usa*female uk*female
			   frn*ses2 frn*ses3 frn*ses4 frn*ses9
               uk*ses2 uk*ses3 uk*ses4 uk*ses9
               usa*ses2 usa*ses3 usa*ses4 usa*ses9 
			   frn*outyear_awcen frn*oyrawcenb_sp1  frn*oyrawcenb_sp2 frn*oyrawcenb_sp3 frn*oyrawcenb_sp4 frn*oyrawcenb_sp5 
			   uk*outyear_awcen uk*oyrawcenb_sp1  uk*oyrawcenb_sp2 uk*oyrawcenb_sp3 uk*oyrawcenb_sp4 uk*oyrawcenb_sp5  
			   usa*outyear_awcen usa*oyrawcenb_sp1 usa*oyrawcenb_sp2 usa*oyrawcenb_sp3 usa*oyrawcenb_sp4 usa*oyrawcenb_sp5               
               frn*outagecen frn*outageawb_sp1 frn*outageawb_sp2 frn*outageawb_sp3
               uk*outagecen uk*outageawb_sp1 uk*outageawb_sp2 uk*outageawb_sp3
               usa*outagecen usa*outageawb_sp1 usa*outageawb_sp2 usa*outageawb_sp3
			   frn*atwork_lag5 frn*atwork_lag10 frn*atwork_lag15
			   usa*atwork_lag5 usa*atwork_lag10 usa*atwork_lag15
			   uk*atwork_lag5 uk*atwork_lag10 uk*atwork_lag15
               frn*dosel_1_5 
               frn*ds_cu_lag6 frn*ds_cu_lag6wb_sp1 
               uk*dosel_1_5 
               uk*ds_cu_lag6 uk*ds_cu_lag6wb_sp1 
               usa*dosel_1_5 
               usa*ds_cu_lag6 usa*ds_cu_lag6wb_sp1 
               ;




*death from CAUSE1;
%LET y1preds = usa uk
			   frn*female usa*female uk*female
			   frn*ses2 frn*ses3 frn*ses4 frn*ses9
               usa*ses2 usa*ses3 usa*ses4 usa*ses9 
               uk*ses2 uk*ses3 uk*ses4 uk*ses9 
               frn*outagecen frn*outagedcanb_sp1 frn*outagedcanb_sp2 frn*outagedcanb_sp3 frn*outagedcanb_sp4 frn*outagedcanb_sp5
               uk*outagecen uk*outagedcanb_sp1 uk*outagedcanb_sp2 uk*outagedcanb_sp3 uk*outagedcanb_sp4 uk*outagedcanb_sp5
               usa*outagecen usa*outagedcanb_sp1 usa*outagedcanb_sp2 usa*outagedcanb_sp3 usa*outagedcanb_sp4 usa*outagedcanb_sp5
               frn*yobcenb_sp1 frn*yobcenb_sp2 frn*yobcenb_sp3 frn*yobcenb_sp4 frn*yobcenb_sp5
               usa*yobcenb_sp1 usa*yobcenb_sp2 usa*yobcenb_sp3 usa*yobcenb_sp4 usa*yobcenb_sp5
               uk*yobcenb_sp1 uk*yobcenb_sp2 uk*yobcenb_sp3 uk*yobcenb_sp4 uk*yobcenb_sp5                              
               dose_agetimecat1 dose_agetimecat2 dose_agetimecat3
			   dose_agetimecat4 dose_agetimecat5 dose_agetimecat6
			   dose_agetimecat7 dose_agetimecat8 dose_agetimecat9
               ;


*death from CAUSE2;
%LET y2preds = usa uk
			   frn*female usa*female uk*female
			   frn*ses2 frn*ses3 frn*ses4 frn*ses9
               usa*ses2 usa*ses3 usa*ses4 usa*ses9 
               uk*ses2 uk*ses3 uk*ses4 uk*ses9
               frn*outagecen frn*outagedcanb_sp1 frn*outagedcanb_sp2 frn*outagedcanb_sp3 frn*outagedcanb_sp4 frn*outagedcanb_sp5
               uk*outagecen uk*outagedcanb_sp1 uk*outagedcanb_sp2 uk*outagedcanb_sp3 uk*outagedcanb_sp4 uk*outagedcanb_sp5
               usa*outagecen usa*outagedcanb_sp1 usa*outagedcanb_sp2 usa*outagedcanb_sp3 usa*outagedcanb_sp4 usa*outagedcanb_sp5
               frn*yobcenb_sp1 frn*yobcenb_sp2 frn*yobcenb_sp3 frn*yobcenb_sp4 frn*yobcenb_sp5
               usa*yobcenb_sp1 usa*yobcenb_sp2 usa*yobcenb_sp3 usa*yobcenb_sp4 usa*yobcenb_sp5
               uk*yobcenb_sp1 uk*yobcenb_sp2 uk*yobcenb_sp3 uk*yobcenb_sp4 uk*yobcenb_sp5               
               dose_agetimecat1 dose_agetimecat2 dose_agetimecat3
			   dose_agetimecat4 dose_agetimecat5 dose_agetimecat6
			   dose_agetimecat7 dose_agetimecat8 dose_agetimecat9
               ;


/*************************************** modeling ************************************************/
%MACRO models(modeldata=an);
  %PUT MODEL: employment status;
  PROC GENMOD DATA =  &modeldata DESCENDING;
   WHERE facop = 1 AND (wyr>0 OR leftwork=1); *admin censoring of work status  by cohort;
   MODEL leftwork = &lPREDS
    / D=B LINK=LOGIT;
   ODS OUTPUT ParameterEstimates=L_model (KEEP = parameter estimate  WHERE=(parameter ^= "Scale"));
  RUN;
  
  %PUT MODEL: lung cancer deaths;
  PROC GENMOD DATA =  &modeldata DESCENDING;
   MODEL d_&CAUSE1 = &y1preds
    / D=B LINK=LOGIT;
   ODS OUTPUT ParameterEstimates=Y1_model (KEEP = parameter estimate  WHERE=(parameter ^= "Scale"));
  RUN;
    
  %PUT MODEL: deaths attributed to non-lung cancer causes ;
  PROC GENMOD DATA =  &modeldata DESCENDING;
   WHERE  d_&CAUSE1 = 0; *conditional on first outcome not happening;
   MODEL d_&CAUSE2 = &y2preds
    / D=B LINK=LOGIT;
   ODS OUTPUT ParameterEstimates=Y2_model (KEEP = parameter estimate WHERE=(parameter ^= "Scale"));
  RUN;
%MEND;

%MODELS(modeldata=an);

%MACRO addcohort(model);
  DATA &model; SET &model; DO COHORT = 5,13,14; OUTPUT; END;RUN;
  PROC SORT DATA=&model; BY cohort;RUN;
%MEND;

%ADDCOHORT(l_model);
%ADDCOHORT(y1_model);
%ADDCOHORT(y2_model);

/************************************************************************************************************************
****** end modeling
************************************************************************************************************************/

*max exposure by cohort and decade (natural caps on exposure);
PROC MEANS DATA = an NOPRINT;
 *BY cohort;
 CLASS cohort indecade;
 VAR logdosert;
 WHERE pyrs > 0.99;
 OUTPUT OUT=maxexp(DROP=_: WHERE=(cohort > .Z AND indecade > .z)) MAX=max_logdosert MIN=min_logdosert;
RUN;

DATA maxexp;
 SET maxexp;
 RETAIN lastmin lastmax 0;
   IF min_logdosert <= .z THEN  min_logdosert = lastmin;
   IF max_logdosert <= .z THEN  max_logdosert = lastmax;
 OUTPUT;
 lastmin = min_logdosert;
 lastmax = max_logdosert;

%AssignMinMaxDose(maxexp);

PROC PRINT DATA = maxexp;
 TITLE "min/max observed log(dose rate) | pyr>0.75";
RUN;

************************************************************************;
*step 2 - MC sample of baseline, macro automation of covariates;
************************************************************************;
PROC SURVEYSELECT DATA = anfirst OUT=mcsample SAMPSIZE=&NMC METHOD=URS OUTHITS;
 TITLE "Initial MC sample of baseline observations";
RUN;


%MACRO count();
%GLOBAL  lmod y1mod y2mod
         nl ny1 ny2;
   /* employment status */
   %LET i=1; %LET lmod=_l1; ;%DO %UNTIL(%QSCAN(&lpreds, %EVAL(&i), " ")=);
      %LET lmod=&lmod + _l%EVAL(&i+1) * %QSCAN(&lpreds, %EVAL(&i), " "); %LET i = %EVAL(&i+1);
    %END;%LET nl = &i;
   
   /* outcome 1 */
   %LET i=1; %LET y1mod=_da1; %DO %UNTIL(%QSCAN(&y1preds, %EVAL(&i), " ")=);
      %LET y1mod=&y1mod + _da%EVAL(&i+1) * %QSCAN(&y1preds, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
   %END;%LET ny1 = &i;
   
   /* outcome(s) 2 */
   %LET i=1; %LET y2mod=_db1; %DO %UNTIL(%QSCAN(&y2preds, %EVAL(&i), " ")=);
      %LET y2mod=&y2mod + _db%EVAL(&i+1) * %QSCAN(&y2preds, %EVAL(&i), " ");  %LET i = %EVAL(&i+1);
   %END;%LET ny2 = &i;
%MEND;

%COUNT;
*employment;
PROC TRANSPOSE DATA = L_model OUT=c_l(DROP=_NAME_); BY cohort; ID parameter; ;RUN;
*mortality;
PROC TRANSPOSE DATA = Y1_model OUT=c_d1(DROP=_NAME_); BY cohort; ID parameter; ;RUN;
PROC TRANSPOSE DATA = Y2_model OUT=c_d2(DROP=_NAME_); BY cohort; ID parameter; ;RUN;


DATA st(KEEP=cohort); SET c_l;
DATA c_l; SET c_l(DROP=cohort);
 ARRAY coefs[*] _NUMERIC_; ARRAY _l[&nl];
 DO j = 1 TO DIM(coefs); _l[j] = coefs[j];END;
 KEEP _l:;
DATA c_d1; SET c_d1(DROP=cohort);
 ARRAY coefs[*] _NUMERIC_;ARRAY _da[&ny1];
 DO j = 1 TO DIM(coefs);_da[j] = coefs[j];END;
 KEEP _da:;
DATA c_d2; SET c_d2(DROP=cohort);
 ARRAY coefs[*] _NUMERIC_;ARRAY _db[&ny2];
 DO j = 1 TO DIM(coefs);_db[j] = coefs[j];END;
 KEEP _db:;
OPTIONS MERGENOBY=NOWARN;
DATA c_l; MERGE st c_l;
DATA c_d1; MERGE st c_d1;
DATA c_d2; MERGE st c_d2; 
RUN;
OPTIONS MERGENOBY=WARN;

************************************************************************;
*step 3 - MC sample of time-varying quantities;
************************************************************************;
%MACRO GFORMULAsim();
  DATA gformula2;
   LENGTH newid 8;
   MERGE mcsample
    c_aA c_a c_l c_d1 c_d2;
    BY cohort;
   newid = _N_;
   OUTPUT gformula2;
  RUN;
  
  DATA gformula2
    ;
   LENGTH cohort gfid newid intervention inage outage d_&CAUSE1 d_&CAUSE2 d_death ltfu censadmin leftwork dose dose_cum atwork anydose wyr wyr_cum agelastwork 8;
   SET gformula2(
                 DROP = d_: 
                      cancer solid leuk can_nole sol_nolu oral eso stom colon rectum liver panc perit larynx lung breast 
                      uterus ovary prost kidney bladder brain thyr conn gallbldr livernog bone pleura allskin testis breastf pleuraex
  					outage outyear totfu
  					numberhits
    );
    DROP _:;
   FORMAT intervention intform.;
  
    /* create containers */
    %VSinits(outcomes=d_&CAUSE1 d_&CAUSE2 ltfu);
    %TVinits(pref=,size=26, agest=16, ageend=85, setzero=0);
    %TVinits(pref=_old,size=26, agest=16, ageend=85, setzero=0); /* note this does not overwrite age specific exposures */
    * store original variables for re-use for an individual across interventions;
    %TRANSFER(varlist= inage indate inyear dose_cum wyr_cum, prefA=_old, prefB=);
    %cumTRANSFER(prefA=_old, prefB=, size=26, agest=16, ageend=85);
    %dosewindowTRANSFER(prefA=_old, prefB=, agepts = &ageendpoints, lagpts= &lagendpoints);
    %TFTRANSFORM(); *transformations on time fixed variables (could be cut?);
    
    /* main loop - carrying out each intervention in each pseudo-individual */
    intnum = 0;
    DO intervention = &INTS;
      gfid = 10*(intnum+1)*(&nmc+1)+newid+1; *are these unique?;
    
      *initializing time-varying variables;
      *employment;
      atwork=1;
      leftwork=0;
      *years of exposure;
      wyr = 1;
  	  pyrs = 1;
      *bring in id specific baseline levels of covariates;
      %TRANSFER(varlist= inage indate inyear dose_cum wyr_cum, prefA=, prefB=_old);
      %CUMTRANSFER(prefA=, prefB=_old, size=26, agest=16, ageend=85);
      %dosewindowTRANSFER(prefA=, prefB=_old, agepts = &ageendpoints, lagpts= &lagendpoints);
  
      *utility;
      firstobs = 1;
      done = 0;
      * time on study;
      intos = 0;
      __startage = inage;
      __startdate = indate;
      *mortality and follow up;
      *MAIN time loop over an individual;
      DO WHILE(done=0);
         *external time varying quantities;
         outage = inage+1;
         outyear = inyear+1;
         outtos = intos + 1;
         /* is facility still in operation ?*/
         facop = 1;
         IF mainfacy IN (11,12) AND inyear >= 2004.9972678 THEN facop = 0;*us;
         ELSE IF mainfacy IN (13) AND inyear >= 2003.9972678 THEN facop = 0;*us;
         ELSE IF mainfacy IN (21, 23, 24, 25, 221, 222) AND inyear >= 2002 THEN facop = 0; * uk;
         ELSE IF mainfacy IN (454, 460) AND inyear >= 1998.9972678 THEN facop = 0; *france;
         ELSE IF mainfacy IN (456) AND inyear >= 1991.9972678 THEN facop = 0;*france;
         ELSE IF mainfacy IN (457) AND inyear >= 1988.1420765027 THEN facop = 0;*france;
         ELSE IF mainfacy IN (458) AND inyear >= 1991.4958904 THEN facop = 0;*france;
         ELSE IF inyear >= 2004.9972678 THEN facop = 0; /* in case there are facilities that do not show up in the test data */       
         
         * all variable transformations on lagged/external time varying variables here;
         %TVtimeTRANSFORM(); *decade created here;
  	     %SetMinMaxDose(clist=&cohortlist, dlist=&decadelist, maxlist=&maxratelist, minlist=&minratelist);
         %TVlagworkTRANSFORM();
         %TVlagdoseTRANSFORM();
      
         ****** internal (modeled) time varying quantities;
          **************;
          /*  employment */
          **************;
         IF firstobs THEN DO;
           firstobs = 0; *do nothing here - every individual employed in first observation;
           wyr = 1;
           atwork = 1;
         END;
         ELSE DO;
          IF facop = 0 THEN atwork = 0;
          IF atwork THEN DO;
           *if already at work (this is effectively lagged one year), simulate probability of leaving work;
           lol = &LMOD;
           %BOUND(lol, upperbound=700, lowerbound=-700);
           atwork = 1 - RAND('BERNOULLI', 1/(1+exp(-lol)));
           wyr = atwork;
           IF atwork=0 THEN DO; leftwork=1; END;
           ELSE leftwork = 0;
           agelastwork=outage; * will only be valid in last observations (but that is ok when this step only outputs last observation);
          END; * if at work;
          ELSE DO;
           *if not already at work, do nothing;
           atwork = 0;
           wyr = 0;
           leftwork = 0;
          END; *if not working;
         END; *if not firstobs;;
         /*%TVlagworkTRANSFORM*/;
         %TVworkTRANSFORM;
          **************;
          /*  exposure */
          **************;
         IF atwork>0 THEN DO;
           *some studies have unexposed - this is the log odds of any exposure;
       
           /* BEGIN INTERVENTIONS */
           %deterministicINTERVENTIONS(); *negative values of intervention = deterministic intervention (conditional on being at work);
           * transform dose if still noted to be exposed;
           IF anydose THEN dose = EXP(logdosert);
           ELSE IF anydose=0 THEN DO; 
             logdosert=.z; 
             dose = 0;
           END;
           **************;
           /*  end exposure */
           **************;
         END;*atwork;
         **************;
         /*  end employment */
         **************;
         IF atwork=0 THEN DO;
          anydose=0;
          logdosert=.z; 
          dose=0;
          wyr=0;
         END;
         wyr_cum = SUM(wyr_cum_lag1, wyr);
         dose_cum = SUM(ds_cu_lag1, dose);

         %CUMdoseCREATE(inline=no,prefa=, prefb=, size=26); * here lag 1;
         %CUMworkCREATE(inline=no,prefa=, prefb=, size=26); * lagged employment indicators;
         %dosewindowCREATE(pref=,size=26, agest=16, ageend=85, agepts = &ageendpoints, lagpts= &lagendpoints, id = gfid);    
         /**/
         %TVdoseTRANSFORM;
         **************;
         /*  outcomes */
         **************;
    
         *mortality, main cause;
         lod1 = &y1mod;
         %BOUND(lod1, upperbound=700, lowerbound=-700);
         pc1 = 1/(1+exp(-lod1));
         d_&CAUSE1 = RAND('BERNOULLI', pc1);
      
         *mortality, complementary cause;
         IF d_&CAUSE1 = 1 THEN d_&CAUSE2 = 0;
         ELSE DO;
           lod2 = &y2mod;
           %BOUND(lod2, upperbound=700, lowerbound=-700);
           pc2 = 1/(1+exp(-lod2));
           d_&CAUSE2 = RAND('BERNOULLI', pc2);
         END;
         
         d_death = d_&CAUSE1+d_&CAUSE2;
         *censoring by loss to follow-up;
         ltfu=0;
         *exit if reach outage 90 or censoring date;
         IF outage >= eofage AND d_death=0 THEN censadmin=1; ELSE censadmin=0;
         IF outage >= 90 OR d_death>0 OR censadmin=1 OR ltfu=1 THEN done=1;
         *do not include those younger than 16 years of age (for stability - this is later reduced to age 20);
         
         *IF inage >= 16 THEN OUTPUT gformula2;
         IF done THEN DO;
             inage = MAX(16, __startage);
             __finalinc = MAX(0, __startage-inage);
             indate = __startdate + __finalinc;
             IF outage >= 16 THEN OUTPUT gformula2;
             
             KEEP
              intervention gfid
             cohort frn usa uk female ses2 ses3 ses4 ses9
             d_&CAUSE1 d_&CAUSE2 d_death censadmin ltfu
             yob inyear outyear inage  outage wyr logdosert wyr_cum dose_cum agelastwork
             dose_agetimecat1 dose_agetimecat2 dose_agetimecat3 
             dose_agetimecat4 dose_agetimecat5 dose_agetimecat6 
             dose_agetimecat7 dose_agetimecat8 dose_agetimecat9
             dosel_1_5 dosel_6_10 dose_lag1
             ;
         END;
         *lagged exposure and employment variables;
         %LAGS(pref=,size=26); * lag by one year: sets lag1_new = current, lag2_new=lag1_current, etc.;
         *incrementing time by one year;
         inage=inage+1; 
         inyear = inyear+1;
         intos = intos+1;
      END;*time loop;
      intnum = intnum+1;
   END; *intervention loop;
  RUN;
  
  
%MEND;

%GFORMULAsim();



RUN;QUIT;RUN;
