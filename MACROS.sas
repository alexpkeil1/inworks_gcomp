/**********************************************************************************************************************
* Author: Alex Keil
* Project: INWORKS (International Nuclear Workers cohort)
* Description: SAS macros that get used in other programs
* Keywords: radiation, occupational, g-formula
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/


%MACRO VSinits(outcomes=d_cancer d_solid d_leuk d_can_nole d_sol_nolu d_oral d_eso d_stom d_colon d_rectum d_liver d_panc d_perit d_larynx d_lung d_breast d_uterus 
                      d_ovary d_prost d_kidney d_bladder d_brain d_thyr d_conn d_gallbldr d_livernog d_bone d_pleura d_allskin d_testis d_breastf d_pleuraex);
  ARRAY d_outcomes[*] &outcomes;
  DO j = 1 TO DIM(d_outcomes);
    d_outcomes[j] = 0;
  END;
%MEND;

%MACRO TVinits(pref=,size=26, agest=16, ageend=85, setzero=1);
  ARRAY &pref.wyr_lag[&size] ;
  ARRAY &pref.wyr_cum_lag[&size] ;
  ARRAY &pref.dose_lag[&size] ;
  ARRAY &pref.ds_cu_lag[&size] ;
  ARRAY &pref.dose_age[&agest : &ageend] &pref.dose_age&agest-&pref.dose_age&ageend;
  *age at exposure/time since exposure 2 dim array;
  ARRAY &pref.dose_agetimecat[9];
   * initialize if setzero = 1;
  %IF (&SETZERO=1) %THEN %DO;
    DO i = 1 TO &size;
      &pref.wyr_lag[i] = 0;
      &pref.wyr_cum_lag[i] = 0;
      &pref.dose_lag[i] = 0;
      &pref.ds_cu_lag[i] = 0;
    END;
	DO i = &AGEST TO &AGEEND;
      &pref.dose_age[i] = 0; /*todo: this should come over from the observed data!!! */
	END;
	DO i = 1 TO 9;
      &pref.dose_agetimecat[i] = 0;
	END;
  %END;
%MEND;

 %MACRO CUMdoseCREATE(inline=,prefa=, prefb=, size=26);
    %IF ((&inline =)) %THEN %DO;
      anydose = (dose>0);
      logdosert = .z; * 3/22;
      IF dose > 0 & wyr > 0 THEN logdosert = LOG(dose/wyr);
      IF dose > 0 & wyr = 0 THEN DO;
        * check for errors in dose;
        PUT "check for ERR0R: wyr=0 " dose= wkerid= inage= inyear= outyear= end_emp_= start_em=;
      END;
      * cumulative variables;
      dose_cum = dose_cum + dose;
    %END;
	*todo: separate lagged exposures from those that rely on current dose;
    * cumulative exposures;
    %IF &prefa ^= &prefb %THEN %DO;
      %IF %EVAL(&size >=1)  %THEN  &prefa.ds_cu_lag1  = &prefB.ds_cu_lag[1];;
      %IF %EVAL(&size >=5)  %THEN  &prefa.ds_cu_lag5  = &prefB.ds_cu_lag[5];;
      %IF %EVAL(&size >=6)  %THEN  &prefa.ds_cu_lag6  = &prefB.ds_cu_lag[6];;
      %IF %EVAL(&size >=10) %THEN &prefa.ds_cu_lag10 = &prefB.ds_cu_lag[10];;
      %IF %EVAL(&size >=15) %THEN &prefa.ds_cu_lag15 = &prefB.ds_cu_lag[15];;
      %IF %EVAL(&size >=20) %THEN &prefa.ds_cu_lag20 = &prefB.ds_cu_lag[20];;
      %IF %EVAL(&size >=25) %THEN &prefa.ds_cu_lag25 = &prefB.ds_cu_lag[25];;
    %END;
 %MEND;

 %MACRO CUMworkCREATE(inline=,prefa=, prefb=, size=26);

   %IF %EVAL(&size >=1)  %THEN   atwork_lag1 = (wyr_lag1>0);;
   %IF %EVAL(&size >=5)  %THEN   atwork_lag5 = (wyr_lag5>0);;
   %IF %EVAL(&size >=10) %THEN  atwork_lag10 = (wyr_lag10>0);;
   %IF %EVAL(&size >=15) %THEN  atwork_lag15 = (wyr_lag15>0);;
   %IF %EVAL(&size >=20) %THEN  atwork_lag20 = (wyr_lag20>0);;
   %IF %EVAL(&size >=25) %THEN  atwork_lag25 = (wyr_lag25>0);;
   * cumulative employment variables;
     %IF ((&inline =)) %THEN %DO;
      * cumulative variables;
      wyr_cum = wyr_cum + wyr;
    %END;
	*todo: separate lagged exposures from those that rely on current dose;
    * cumulative working;
    %IF &prefa ^= &prefb %THEN %DO;
      %IF %EVAL(&size >=1)  %THEN  &prefa.wyr_cum_lag1  = &prefB.wyr_cum_lag[1];;
      %IF %EVAL(&size >=5)  %THEN  &prefa.wyr_cum_lag5  = &prefB.wyr_cum_lag[5];;
      %IF %EVAL(&size >=10) %THEN &prefa.wyr_cum_lag10 = &prefB.wyr_cum_lag[10];;
      %IF %EVAL(&size >=15) %THEN &prefa.wyr_cum_lag15 = &prefB.wyr_cum_lag[15];;
      %IF %EVAL(&size >=20) %THEN &prefa.wyr_cum_lag20 = &prefB.wyr_cum_lag[20];;
      %IF %EVAL(&size >=25) %THEN &prefa.wyr_cum_lag25 = &prefB.wyr_cum_lag[25];;
    %END;
%MEND;

 %MACRO dosewindowCREATE(pref, size=26, agest=16, ageend=85, 
    /* */
    agepts=0 29 44 100, /* grant lists 30, 30-45, 45+ */
	/* ages in (min,max] using integer years*/
    /* todo: fix dose accumulation*/
    /* todo: set work coefs to zero in sim*/
    lagpts=9 19 29 100,   /* grant lists 10-20, 20-30, 30+ */
	/* lags in (min,max] using integer years*/
    id=wkerid);
    * time since exposure windows;
    %IF %EVAL(&size >=5)  %THEN &pref.dosel_1_5 =   SUM(OF &pref.dose_lag1 - &pref.dose_lag5);;
    %IF %EVAL(&size >=5)  %THEN &pref.dosel_2_5 =   SUM(OF &pref.dose_lag2 - &pref.dose_lag5);;
    %IF %EVAL(&size >=10) %THEN &pref.dosel_6_10 =  SUM(OF &pref.dose_lag6 - &pref.dose_lag10);;
    %IF %EVAL(&size >=15) %THEN &pref.dosel_11_15 = SUM(OF &pref.dose_lag11 - &pref.dose_lag15);;
    %IF %EVAL(&size >=20) %THEN &pref.dosel_16_20 = SUM(OF &pref.dose_lag16 - &pref.dose_lag20);;
    %IF %EVAL(&size >=20) %THEN &pref.dosel_10_20 = SUM(OF &pref.dose_lag10 - &pref.dose_lag20);;
    %IF %EVAL(&size >=25) %THEN &pref.dosel_21_25 = SUM(OF &pref.dose_lag21 - &pref.dose_lag25);;
    * age at exposure windows (100% of annual dose assigned to maximum age in record, rather than split across ages);
    RETAIN holdid holddose holdage 0;
	agecurr = MIN(&ageend, MAX(&agest, FLOOR(inage)));
	IF (FLOOR(inage) > holdage OR holdid ^= &id) THEN holddose = dose; ELSE holddose = holddose + dose;
    &pref.dose_age[agecurr] =  holddose;
    holdage = FLOOR(inage);
    holdid = &id;
    
	%LET __J = 1;
	%LET ageidx = 1;
	%DO %WHILE(%SCAN(&agepts, &ageidx+1) NE );
	* define age + lag constraints on a time window relative to the current age (e.g. 55);
       * e.g. potentially accumulate exposure from age (0,29] = [1,29] if ageidx = 1;	
      %LET lagidx = 1;
      %DO %WHILE(%SCAN(&lagpts, &lagidx+1) NE );
		  &pref.dose_agetimecat[&__J] = 0;
		* possible lag window for dose accumulation;
	    _lagmin =  %SCAN(&lagpts, %EVAL(&lagidx+0));  * e.g. 29 years prior if lagidx=3;
	    _lagmax =  %SCAN(&lagpts, %EVAL(&lagidx+1));  * e.g. 100 years prior if lagidx = 3;
		* possible age window for dose accumulation, given the current lag;
	    _agelagmin = MAX(&agest,  agecurr - _lagmax + 1); * e.g. current age 58.9 - 99 = max(-41.9 years, 16 years) if lagidx=3;
	    _agelagmax = MIN(&ageend, agecurr - _lagmin - 0); * e.g. current age 58.9 - 29 = 29.9 years if lagidx=3;
	    * now check whether that age range is within the bounds of the current window under consideration;
		IF _agelagmin <= _agelagmax THEN DO;
 		  * age at min lag is higher than upper age of window;
 		  IF _agelagmax > %SCAN(&agepts, &ageidx+1) THEN _agelagmax = %SCAN(&agepts, &ageidx+1); * e.g. shrink upper limit from 29.9 years to 29 years - this actually will pull all dose up until 30th birthday;
 		  * age at max lag is younger than lowest age of window;
		  IF _agelagmin <= %SCAN(&agepts, &ageidx) THEN _agelagmin = %SCAN(&agepts, &ageidx)+1; * e.g. change lower limit from -41.9 years to 1 years (really 16);
		  * now accumulate dose within the lag/age window that fits for the current observation;
		  DO _thisage = _agelagmin TO _agelagmax; * accumulate exposure from age 1 to 29 = [1,30);
		   * note this loop will do nothing (short circuit) if _agelagmax< _agelagmin;
            &pref.dose_agetimecat[&__J] = &pref.dose_agetimecat[&__J] + &pref.dose_age[_thisage]; * slight accumulation err0r here from bringing in annual dose for partial year (empirically, this does not matter much);
		  END;
		END;
		*OUTPUT; *only for testing!;
	    %LET __J = %EVAL(&__J + 1);
        %LET lagidx = %EVAL(&lagidx + 1);
	  %END;
      %LET ageidx = %EVAL(&ageidx + 1);
	%END;
    &pref.dose_agetimecatcheck = SUM(OF &pref.dose_agetimecat1-&pref.dose_agetimecat9);
 %MEND;

 %MACRO lags(pref=,size=26);
    DO i = &size TO 2 BY -1;
      &pref.wyr_lag[i]    = &pref.wyr_lag[i-1];
      &pref.wyr_cum_lag[i] = &pref.wyr_cum_lag[i-1];
      &pref.dose_lag[i]       = &pref.dose_lag[i-1];
      &pref.ds_cu_lag[i]    = &pref.ds_cu_lag[i-1];
    END;
    &pref.wyr_lag[1] = wyr;
    &pref.wyr_cum_lag[1] = wyr_cum;
    &pref.dose_lag[1] = dose;
    &pref.ds_cu_lag[1] = dose_cum;
%MEND;

%MACRO cumTRANSFER(prefA=, prefB=, size=26, agest=16, ageend=85);
    * setting &prefA variables to values of &prefB (cumulative);
    DO __i = 1 TO &size;
     &prefA.dose_lag[__i]  = &prefB.dose_lag[__i];
     &prefA.ds_cu_lag[__i]  = &prefB.ds_cu_lag[__i];
     &prefA.wyr_cum_lag[__i]  = &prefB.wyr_cum_lag[__i];
     &prefA.wyr_lag[__i]  = &prefB.wyr_lag[__i];
    END;
    DO i = &AGEST TO &AGEEND;
      &prefA.dose_age[i] =  &prefB.dose_age[i];
	END;
 %MEND;

 
%MACRO dosewindowTRANSFER(prefA=, prefB=, 
    agepts=0 35 45 100, 
    lagpts=0 5 15 100
    );
    * setting &prefA variables to values of &prefB (age/time dose windows);
	%LET __J = 1;
	%LET ageidx = 1;
	%DO %WHILE(%SCAN(&agepts, &ageidx+1) NE );
	/*%DO ageidx = 1 %TO 3;*/
      %LET lagidx = 1;
      %DO %WHILE(%SCAN(&lagpts, &lagidx+1) NE );
		&prefa.dose_agetimecat[&__J] = &prefb.dose_agetimecat[&__J];
	    %LET __J = %EVAL(&__J + 1);
        %LET lagidx = %EVAL(&lagidx + 1);
	  %END;
      %LET ageidx = %EVAL(&ageidx + 1);
	%END;
%MEND;

%MACRO transfer(varlist=, prefA=, prefB=);
     * setting &prefA variables to values of &prefB (&varlist);
   %LET __j = 1;
    %DO %WHILE (%SCAN(&varlist, %EVAL(&__j), " ")^=);
      %LET var = %SCAN(&varlist, %EVAL(&__j), " ");
      &prefA.&var = &prefb.&var;
      %LET __j=%EVAL(&__j+1);
	%END;
%MEND;

%MACRO lags(pref=,size=26);
    DO i = &size TO 2 BY -1;
      &pref.wyr_lag[i]    = &pref.wyr_lag[i-1];
      &pref.wyr_cum_lag[i] = &pref.wyr_cum_lag[i-1];
      &pref.dose_lag[i]       = &pref.dose_lag[i-1];
      &pref.ds_cu_lag[i]    = &pref.ds_cu_lag[i-1];
    END;
    &pref.wyr_lag[1] = wyr;
    &pref.wyr_cum_lag[1] = wyr_cum;
    &pref.dose_lag[1] = dose;
    &pref.ds_cu_lag[1] = dose_cum;
%MEND;

%MACRO TFTRANSFORM();
  *transformations of time fixed quantities;
    ******** centered variable(s) ********;
    yobcen = yob - 1900;
    yobcat1 = (1800 <= yob < 1920);
    yobcat2 = (1920 <= yob < 1930);
    yobcat3 = (1930 <= yob < 1940);
    yobcat4 = (1940 <= yob < 1950);
    yobcat5 = (1950 <= yob        );
    **** country specific year of birth term splines;
    IF cohort=5 THEN DO;
      * yobcen knots: 1.306507 19.64918 25.46617 30.47397 34.70722 45.16206 62.2237 ;
      yobcenb_sp1 = ((yobcen>1.30650684931497) * (yobcen-1.30650684931497)/((62.223698630137-1.30650684931497)**(2/3)))**3
          + ((45.1620605810466-1.30650684931497)* ((yobcen>62.223698630137) * (yobcen-62.223698630137)/((62.223698630137-1.30650684931497)**(2/3)))**3 - 
             (62.223698630137-1.30650684931497) * ((yobcen>45.1620605810466) * (yobcen-45.1620605810466)/((62.223698630137-1.30650684931497)**(2/3)))**3)/(62.223698630137-45.1620605810466);
      yobcenb_sp2 = ((yobcen>19.6491750986302) * (yobcen-19.6491750986302)/((62.223698630137-1.30650684931497)**(2/3)))**3
          + ((45.1620605810466-19.6491750986302)* ((yobcen>62.223698630137) * (yobcen-62.223698630137)/((62.223698630137-1.30650684931497)**(2/3)))**3 - 
             (62.223698630137-19.6491750986302) * ((yobcen>45.1620605810466) * (yobcen-45.1620605810466)/((62.223698630137-1.30650684931497)**(2/3)))**3)/(62.223698630137-45.1620605810466);
      yobcenb_sp3 = ((yobcen>25.4661663726028) * (yobcen-25.4661663726028)/((62.223698630137-1.30650684931497)**(2/3)))**3
          + ((45.1620605810466-25.4661663726028)* ((yobcen>62.223698630137) * (yobcen-62.223698630137)/((62.223698630137-1.30650684931497)**(2/3)))**3 - 
             (62.223698630137-25.4661663726028) * ((yobcen>45.1620605810466) * (yobcen-45.1620605810466)/((62.223698630137-1.30650684931497)**(2/3)))**3)/(62.223698630137-45.1620605810466);
      yobcenb_sp4 = ((yobcen>30.4739726027398) * (yobcen-30.4739726027398)/((62.223698630137-1.30650684931497)**(2/3)))**3
          + ((45.1620605810466-30.4739726027398)* ((yobcen>62.223698630137) * (yobcen-62.223698630137)/((62.223698630137-1.30650684931497)**(2/3)))**3 - 
             (62.223698630137-30.4739726027398) * ((yobcen>45.1620605810466) * (yobcen-45.1620605810466)/((62.223698630137-1.30650684931497)**(2/3)))**3)/(62.223698630137-45.1620605810466);
      yobcenb_sp5 = ((yobcen>34.7072216410959) * (yobcen-34.7072216410959)/((62.223698630137-1.30650684931497)**(2/3)))**3
          + ((45.1620605810466-34.7072216410959)* ((yobcen>62.223698630137) * (yobcen-62.223698630137)/((62.223698630137-1.30650684931497)**(2/3)))**3 - 
             (62.223698630137-34.7072216410959) * ((yobcen>45.1620605810466) * (yobcen-45.1620605810466)/((62.223698630137-1.30650684931497)**(2/3)))**3)/(62.223698630137-45.1620605810466);
    END;
    ELSE IF cohort=13 THEN DO;
      * yobcen knots: -4.396409 8.531778 14.51033 20.1612 24.85186 30.21061 45.27392 ;
      yobcenb_sp1 = ((yobcen>-4.3964086009431) * (yobcen--4.3964086009431)/((45.2739243206826--4.3964086009431)**(2/3)))**3
          + ((30.2106138520549--4.3964086009431)* ((yobcen>45.2739243206826) * (yobcen-45.2739243206826)/((45.2739243206826--4.3964086009431)**(2/3)))**3 - 
             (45.2739243206826--4.3964086009431) * ((yobcen>30.2106138520549) * (yobcen-30.2106138520549)/((45.2739243206826--4.3964086009431)**(2/3)))**3)/(45.2739243206826-30.2106138520549);
      yobcenb_sp2 = ((yobcen>8.5317782349726) * (yobcen-8.5317782349726)/((45.2739243206826--4.3964086009431)**(2/3)))**3
          + ((30.2106138520549-8.5317782349726)* ((yobcen>45.2739243206826) * (yobcen-45.2739243206826)/((45.2739243206826--4.3964086009431)**(2/3)))**3 - 
             (45.2739243206826-8.5317782349726) * ((yobcen>30.2106138520549) * (yobcen-30.2106138520549)/((45.2739243206826--4.3964086009431)**(2/3)))**3)/(45.2739243206826-30.2106138520549);
      yobcenb_sp3 = ((yobcen>14.5103291616438) * (yobcen-14.5103291616438)/((45.2739243206826--4.3964086009431)**(2/3)))**3
          + ((30.2106138520549-14.5103291616438)* ((yobcen>45.2739243206826) * (yobcen-45.2739243206826)/((45.2739243206826--4.3964086009431)**(2/3)))**3 - 
             (45.2739243206826-14.5103291616438) * ((yobcen>30.2106138520549) * (yobcen-30.2106138520549)/((45.2739243206826--4.3964086009431)**(2/3)))**3)/(45.2739243206826-30.2106138520549);
      yobcenb_sp4 = ((yobcen>20.1612021857924) * (yobcen-20.1612021857924)/((45.2739243206826--4.3964086009431)**(2/3)))**3
          + ((30.2106138520549-20.1612021857924)* ((yobcen>45.2739243206826) * (yobcen-45.2739243206826)/((45.2739243206826--4.3964086009431)**(2/3)))**3 - 
             (45.2739243206826-20.1612021857924) * ((yobcen>30.2106138520549) * (yobcen-30.2106138520549)/((45.2739243206826--4.3964086009431)**(2/3)))**3)/(45.2739243206826-30.2106138520549);
      yobcenb_sp5 = ((yobcen>24.8518637650272) * (yobcen-24.8518637650272)/((45.2739243206826--4.3964086009431)**(2/3)))**3
          + ((30.2106138520549-24.8518637650272)* ((yobcen>45.2739243206826) * (yobcen-45.2739243206826)/((45.2739243206826--4.3964086009431)**(2/3)))**3 - 
             (45.2739243206826-24.8518637650272) * ((yobcen>30.2106138520549) * (yobcen-30.2106138520549)/((45.2739243206826--4.3964086009431)**(2/3)))**3)/(45.2739243206826-30.2106138520549);
    END;
    ELSE IF cohort=14 THEN DO;
      * yobcen knots: -5.223425 7.616071 14.84966 18.69452 23.82171 29.31956 47.26288 ;
      yobcenb_sp1 = ((yobcen>-5.22342465753422) * (yobcen--5.22342465753422)/((47.2628767123287--5.22342465753422)**(2/3)))**3
          + ((29.3195554904109--5.22342465753422)* ((yobcen>47.2628767123287) * (yobcen-47.2628767123287)/((47.2628767123287--5.22342465753422)**(2/3)))**3 - 
             (47.2628767123287--5.22342465753422) * ((yobcen>29.3195554904109) * (yobcen-29.3195554904109)/((47.2628767123287--5.22342465753422)**(2/3)))**3)/(47.2628767123287-29.3195554904109);
      yobcenb_sp2 = ((yobcen>7.61607050958894) * (yobcen-7.61607050958894)/((47.2628767123287--5.22342465753422)**(2/3)))**3
          + ((29.3195554904109-7.61607050958894)* ((yobcen>47.2628767123287) * (yobcen-47.2628767123287)/((47.2628767123287--5.22342465753422)**(2/3)))**3 - 
             (47.2628767123287-7.61607050958894) * ((yobcen>29.3195554904109) * (yobcen-29.3195554904109)/((47.2628767123287--5.22342465753422)**(2/3)))**3)/(47.2628767123287-29.3195554904109);
      yobcenb_sp3 = ((yobcen>14.8496555369864) * (yobcen-14.8496555369864)/((47.2628767123287--5.22342465753422)**(2/3)))**3
          + ((29.3195554904109-14.8496555369864)* ((yobcen>47.2628767123287) * (yobcen-47.2628767123287)/((47.2628767123287--5.22342465753422)**(2/3)))**3 - 
             (47.2628767123287-14.8496555369864) * ((yobcen>29.3195554904109) * (yobcen-29.3195554904109)/((47.2628767123287--5.22342465753422)**(2/3)))**3)/(47.2628767123287-29.3195554904109);
      yobcenb_sp4 = ((yobcen>18.6945205479452) * (yobcen-18.6945205479452)/((47.2628767123287--5.22342465753422)**(2/3)))**3
          + ((29.3195554904109-18.6945205479452)* ((yobcen>47.2628767123287) * (yobcen-47.2628767123287)/((47.2628767123287--5.22342465753422)**(2/3)))**3 - 
             (47.2628767123287-18.6945205479452) * ((yobcen>29.3195554904109) * (yobcen-29.3195554904109)/((47.2628767123287--5.22342465753422)**(2/3)))**3)/(47.2628767123287-29.3195554904109);
      yobcenb_sp5 = ((yobcen>23.8217056356164) * (yobcen-23.8217056356164)/((47.2628767123287--5.22342465753422)**(2/3)))**3
          + ((29.3195554904109-23.8217056356164)* ((yobcen>47.2628767123287) * (yobcen-47.2628767123287)/((47.2628767123287--5.22342465753422)**(2/3)))**3 - 
             (47.2628767123287-23.8217056356164) * ((yobcen>29.3195554904109) * (yobcen-29.3195554904109)/((47.2628767123287--5.22342465753422)**(2/3)))**3)/(47.2628767123287-29.3195554904109);
    END;   
%MEND;

 
%MACRO TVtimeTRANSFORM();
  *transformations of external time varying quantities (e.g. age, date);
    ******** centered variable(s) ********;
    outyear_awcen = outyear-1960;
    indecade = FLOOR(inyear/10)*10;
    outagecen = outage-45;
    outyear_cen = outyear-1960;
    ybya = yobcen*(outage-60);

    **** country specific calendar time splines;
    IF cohort=5 THEN DO;
      * outyear_awcen knots: 12 22 30 36 41 44.13666 44.55191 ;
      oyrawcenb_sp1 = ((outyear_awcen>12) * (outyear_awcen-12)/((44.5519125683061-12)**(2/3)))**3
          + ((44.1366581502732-12)* ((outyear_awcen>44.5519125683061) * (outyear_awcen-44.5519125683061)/((44.5519125683061-12)**(2/3)))**3 - 
             (44.5519125683061-12) * ((outyear_awcen>44.1366581502732) * (outyear_awcen-44.1366581502732)/((44.5519125683061-12)**(2/3)))**3)/(44.5519125683061-44.1366581502732);
      oyrawcenb_sp2 = ((outyear_awcen>22) * (outyear_awcen-22)/((44.5519125683061-12)**(2/3)))**3
          + ((44.1366581502732-22)* ((outyear_awcen>44.5519125683061) * (outyear_awcen-44.5519125683061)/((44.5519125683061-12)**(2/3)))**3 - 
             (44.5519125683061-22) * ((outyear_awcen>44.1366581502732) * (outyear_awcen-44.1366581502732)/((44.5519125683061-12)**(2/3)))**3)/(44.5519125683061-44.1366581502732);
      oyrawcenb_sp3 = ((outyear_awcen>30) * (outyear_awcen-30)/((44.5519125683061-12)**(2/3)))**3
          + ((44.1366581502732-30)* ((outyear_awcen>44.5519125683061) * (outyear_awcen-44.5519125683061)/((44.5519125683061-12)**(2/3)))**3 - 
             (44.5519125683061-30) * ((outyear_awcen>44.1366581502732) * (outyear_awcen-44.1366581502732)/((44.5519125683061-12)**(2/3)))**3)/(44.5519125683061-44.1366581502732);
      oyrawcenb_sp4 = ((outyear_awcen>36) * (outyear_awcen-36)/((44.5519125683061-12)**(2/3)))**3
          + ((44.1366581502732-36)* ((outyear_awcen>44.5519125683061) * (outyear_awcen-44.5519125683061)/((44.5519125683061-12)**(2/3)))**3 - 
             (44.5519125683061-36) * ((outyear_awcen>44.1366581502732) * (outyear_awcen-44.1366581502732)/((44.5519125683061-12)**(2/3)))**3)/(44.5519125683061-44.1366581502732);
      oyrawcenb_sp5 = ((outyear_awcen>41) * (outyear_awcen-41)/((44.5519125683061-12)**(2/3)))**3
          + ((44.1366581502732-41)* ((outyear_awcen>44.5519125683061) * (outyear_awcen-44.5519125683061)/((44.5519125683061-12)**(2/3)))**3 - 
             (44.5519125683061-41) * ((outyear_awcen>44.1366581502732) * (outyear_awcen-44.1366581502732)/((44.5519125683061-12)**(2/3)))**3)/(44.5519125683061-44.1366581502732);
    END;
    ELSE IF cohort=13 THEN DO;
      * outyear_awcen knots: -1 12 19 26 31 36 41.74103 ;
      oyrawcenb_sp1 = ((outyear_awcen>-1) * (outyear_awcen--1)/((41.7410273972602--1)**(2/3)))**3
          + ((36--1)* ((outyear_awcen>41.7410273972602) * (outyear_awcen-41.7410273972602)/((41.7410273972602--1)**(2/3)))**3 - 
             (41.7410273972602--1) * ((outyear_awcen>36) * (outyear_awcen-36)/((41.7410273972602--1)**(2/3)))**3)/(41.7410273972602-36);
      oyrawcenb_sp2 = ((outyear_awcen>12) * (outyear_awcen-12)/((41.7410273972602--1)**(2/3)))**3
          + ((36-12)* ((outyear_awcen>41.7410273972602) * (outyear_awcen-41.7410273972602)/((41.7410273972602--1)**(2/3)))**3 - 
             (41.7410273972602-12) * ((outyear_awcen>36) * (outyear_awcen-36)/((41.7410273972602--1)**(2/3)))**3)/(41.7410273972602-36);
      oyrawcenb_sp3 = ((outyear_awcen>19) * (outyear_awcen-19)/((41.7410273972602--1)**(2/3)))**3
          + ((36-19)* ((outyear_awcen>41.7410273972602) * (outyear_awcen-41.7410273972602)/((41.7410273972602--1)**(2/3)))**3 - 
             (41.7410273972602-19) * ((outyear_awcen>36) * (outyear_awcen-36)/((41.7410273972602--1)**(2/3)))**3)/(41.7410273972602-36);
      oyrawcenb_sp4 = ((outyear_awcen>26) * (outyear_awcen-26)/((41.7410273972602--1)**(2/3)))**3
          + ((36-26)* ((outyear_awcen>41.7410273972602) * (outyear_awcen-41.7410273972602)/((41.7410273972602--1)**(2/3)))**3 - 
             (41.7410273972602-26) * ((outyear_awcen>36) * (outyear_awcen-36)/((41.7410273972602--1)**(2/3)))**3)/(41.7410273972602-36);
      oyrawcenb_sp5 = ((outyear_awcen>31) * (outyear_awcen-31)/((41.7410273972602--1)**(2/3)))**3
          + ((36-31)* ((outyear_awcen>41.7410273972602) * (outyear_awcen-41.7410273972602)/((41.7410273972602--1)**(2/3)))**3 - 
             (41.7410273972602-31) * ((outyear_awcen>36) * (outyear_awcen-36)/((41.7410273972602--1)**(2/3)))**3)/(41.7410273972602-36);
    END;
    ELSE IF cohort=14 THEN DO;
      * outyear_awcen knots: -8 2.681144 10 18 25 31 37 ;
      oyrawcenb_sp1 = ((outyear_awcen>-8) * (outyear_awcen--8)/((37--8)**(2/3)))**3
          + ((31--8)* ((outyear_awcen>37) * (outyear_awcen-37)/((37--8)**(2/3)))**3 - 
             (37--8) * ((outyear_awcen>31) * (outyear_awcen-31)/((37--8)**(2/3)))**3)/(37-31);
      oyrawcenb_sp2 = ((outyear_awcen>2.68114390684929) * (outyear_awcen-2.68114390684929)/((37--8)**(2/3)))**3
          + ((31-2.68114390684929)* ((outyear_awcen>37) * (outyear_awcen-37)/((37--8)**(2/3)))**3 - 
             (37-2.68114390684929) * ((outyear_awcen>31) * (outyear_awcen-31)/((37--8)**(2/3)))**3)/(37-31);
      oyrawcenb_sp3 = ((outyear_awcen>10) * (outyear_awcen-10)/((37--8)**(2/3)))**3
          + ((31-10)* ((outyear_awcen>37) * (outyear_awcen-37)/((37--8)**(2/3)))**3 - 
             (37-10) * ((outyear_awcen>31) * (outyear_awcen-31)/((37--8)**(2/3)))**3)/(37-31);
      oyrawcenb_sp4 = ((outyear_awcen>18) * (outyear_awcen-18)/((37--8)**(2/3)))**3
          + ((31-18)* ((outyear_awcen>37) * (outyear_awcen-37)/((37--8)**(2/3)))**3 - 
             (37-18) * ((outyear_awcen>31) * (outyear_awcen-31)/((37--8)**(2/3)))**3)/(37-31);
      oyrawcenb_sp5 = ((outyear_awcen>25) * (outyear_awcen-25)/((37--8)**(2/3)))**3
          + ((31-25)* ((outyear_awcen>37) * (outyear_awcen-37)/((37--8)**(2/3)))**3 - 
             (37-25) * ((outyear_awcen>31) * (outyear_awcen-31)/((37--8)**(2/3)))**3)/(37-31);
    END;
    
    ******** splines on attained age ********;

    **** country specific age splines: age at work;
    IF cohort=5 THEN DO;
      * French age spline: age at work;
      * outagecen knots: -17.50602 -2.446366 7.356778 12.98439 17.81356 ;
      outageawb_sp1 = ((outagecen>-17.506018414552) * (outagecen--17.506018414552)/((17.8135642637922--17.506018414552)**(2/3)))**3
          + ((12.9843904858148--17.506018414552)* ((outagecen>17.8135642637922) * (outagecen-17.8135642637922)/((17.8135642637922--17.506018414552)**(2/3)))**3 - 
             (17.8135642637922--17.506018414552) * ((outagecen>12.9843904858148) * (outagecen-12.9843904858148)/((17.8135642637922--17.506018414552)**(2/3)))**3)/(17.8135642637922-12.9843904858148);    outageawb_sp2 = ((outagecen>-2.44636593307882) * (outagecen--2.44636593307882)/((17.8135642637922--17.506018414552)**(2/3)))**3
          + ((12.9843904858148--2.44636593307882)* ((outagecen>17.8135642637922) * (outagecen-17.8135642637922)/((17.8135642637922--17.506018414552)**(2/3)))**3 - 
             (17.8135642637922--2.44636593307882) * ((outagecen>12.9843904858148) * (outagecen-12.9843904858148)/((17.8135642637922--17.506018414552)**(2/3)))**3)/(17.8135642637922-12.9843904858148);    outageawb_sp3 = ((outagecen>7.35677820196123) * (outagecen-7.35677820196123)/((17.8135642637922--17.506018414552)**(2/3)))**3
          + ((12.9843904858148-7.35677820196123)* ((outagecen>17.8135642637922) * (outagecen-17.8135642637922)/((17.8135642637922--17.506018414552)**(2/3)))**3 - 
             (17.8135642637922-7.35677820196123) * ((outagecen>12.9843904858148) * (outagecen-12.9843904858148)/((17.8135642637922--17.506018414552)**(2/3)))**3)/(17.8135642637922-12.9843904858148);
    END;
    ELSE IF cohort=13 THEN DO;
      * UK age spline: age at work;
      * outagecen knots: -21.54858 -14.33835 -5.10137 7.502055 20.93432 ;
      outageawb_sp1 = ((outagecen>-21.5485814806497) * (outagecen--21.5485814806497)/((20.9343214312449--21.5485814806497)**(2/3)))**3
          + ((7.50205479452055--21.5485814806497)* ((outagecen>20.9343214312449) * (outagecen-20.9343214312449)/((20.9343214312449--21.5485814806497)**(2/3)))**3 - 
             (20.9343214312449--21.5485814806497) * ((outagecen>7.50205479452055) * (outagecen-7.50205479452055)/((20.9343214312449--21.5485814806497)**(2/3)))**3)/(20.9343214312449-7.50205479452055);    outageawb_sp2 = ((outagecen>-14.338354292986) * (outagecen--14.338354292986)/((20.9343214312449--21.5485814806497)**(2/3)))**3
          + ((7.50205479452055--14.338354292986)* ((outagecen>20.9343214312449) * (outagecen-20.9343214312449)/((20.9343214312449--21.5485814806497)**(2/3)))**3 - 
             (20.9343214312449--14.338354292986) * ((outagecen>7.50205479452055) * (outagecen-7.50205479452055)/((20.9343214312449--21.5485814806497)**(2/3)))**3)/(20.9343214312449-7.50205479452055);    outageawb_sp3 = ((outagecen>-5.1013698630137) * (outagecen--5.1013698630137)/((20.9343214312449--21.5485814806497)**(2/3)))**3
          + ((7.50205479452055--5.1013698630137)* ((outagecen>20.9343214312449) * (outagecen-20.9343214312449)/((20.9343214312449--21.5485814806497)**(2/3)))**3 - 
             (20.9343214312449--5.1013698630137) * ((outagecen>7.50205479452055) * (outagecen-7.50205479452055)/((20.9343214312449--21.5485814806497)**(2/3)))**3)/(20.9343214312449-7.50205479452055);
    END;
    ELSE IF cohort=14 THEN DO;
      * US age spline: age at work;
      * outagecen knots: -19.97972 -10.15592 -1.054652 10.75507 21.70701 ;
      outageawb_sp1 = ((outagecen>-19.9797245302792) * (outagecen--19.9797245302792)/((21.7070139980537--19.9797245302792)**(2/3)))**3
          + ((10.7550714873868--19.9797245302792)* ((outagecen>21.7070139980537) * (outagecen-21.7070139980537)/((21.7070139980537--19.9797245302792)**(2/3)))**3 - 
             (21.7070139980537--19.9797245302792) * ((outagecen>10.7550714873868) * (outagecen-10.7550714873868)/((21.7070139980537--19.9797245302792)**(2/3)))**3)/(21.7070139980537-10.7550714873868);    outageawb_sp2 = ((outagecen>-10.1559188562018) * (outagecen--10.1559188562018)/((21.7070139980537--19.9797245302792)**(2/3)))**3
          + ((10.7550714873868--10.1559188562018)* ((outagecen>21.7070139980537) * (outagecen-21.7070139980537)/((21.7070139980537--19.9797245302792)**(2/3)))**3 - 
             (21.7070139980537--10.1559188562018) * ((outagecen>10.7550714873868) * (outagecen-10.7550714873868)/((21.7070139980537--19.9797245302792)**(2/3)))**3)/(21.7070139980537-10.7550714873868);    outageawb_sp3 = ((outagecen>-1.05465229433341) * (outagecen--1.05465229433341)/((21.7070139980537--19.9797245302792)**(2/3)))**3
          + ((10.7550714873868--1.05465229433341)* ((outagecen>21.7070139980537) * (outagecen-21.7070139980537)/((21.7070139980537--19.9797245302792)**(2/3)))**3 - 
             (21.7070139980537--1.05465229433341) * ((outagecen>10.7550714873868) * (outagecen-10.7550714873868)/((21.7070139980537--19.9797245302792)**(2/3)))**3)/(21.7070139980537-10.7550714873868);
    END;
    
    **** country specific age splines (all person time);
=    IF cohort=5 THEN DO;
      * French age spline: mortality;
      * outagecen knots: -13.47386 2.744547 12.40614 17.97617 23.08442 30.25848 43.10336 ;
      outagedcanb_sp1 = ((outagecen>-13.4738625645632) * (outagecen--13.4738625645632)/((43.1033561643836--13.4738625645632)**(2/3)))**3
          + ((30.2584796537241--13.4738625645632)* ((outagecen>43.1033561643836) * (outagecen-43.1033561643836)/((43.1033561643836--13.4738625645632)**(2/3)))**3 - 
             (43.1033561643836--13.4738625645632) * ((outagecen>30.2584796537241) * (outagecen-30.2584796537241)/((43.1033561643836--13.4738625645632)**(2/3)))**3)/(43.1033561643836-30.2584796537241);
      outagedcanb_sp2 = ((outagecen>2.74454694342391) * (outagecen-2.74454694342391)/((43.1033561643836--13.4738625645632)**(2/3)))**3
          + ((30.2584796537241-2.74454694342391)* ((outagecen>43.1033561643836) * (outagecen-43.1033561643836)/((43.1033561643836--13.4738625645632)**(2/3)))**3 - 
             (43.1033561643836-2.74454694342391) * ((outagecen>30.2584796537241) * (outagecen-30.2584796537241)/((43.1033561643836--13.4738625645632)**(2/3)))**3)/(43.1033561643836-30.2584796537241);
      outagedcanb_sp3 = ((outagecen>12.4061360465828) * (outagecen-12.4061360465828)/((43.1033561643836--13.4738625645632)**(2/3)))**3
          + ((30.2584796537241-12.4061360465828)* ((outagecen>43.1033561643836) * (outagecen-43.1033561643836)/((43.1033561643836--13.4738625645632)**(2/3)))**3 - 
             (43.1033561643836-12.4061360465828) * ((outagecen>30.2584796537241) * (outagecen-30.2584796537241)/((43.1033561643836--13.4738625645632)**(2/3)))**3)/(43.1033561643836-30.2584796537241);
      outagedcanb_sp4 = ((outagecen>17.9761696234748) * (outagecen-17.9761696234748)/((43.1033561643836--13.4738625645632)**(2/3)))**3
          + ((30.2584796537241-17.9761696234748)* ((outagecen>43.1033561643836) * (outagecen-43.1033561643836)/((43.1033561643836--13.4738625645632)**(2/3)))**3 - 
             (43.1033561643836-17.9761696234748) * ((outagecen>30.2584796537241) * (outagecen-30.2584796537241)/((43.1033561643836--13.4738625645632)**(2/3)))**3)/(43.1033561643836-30.2584796537241);
      outagedcanb_sp5 = ((outagecen>23.0844223726027) * (outagecen-23.0844223726027)/((43.1033561643836--13.4738625645632)**(2/3)))**3
          + ((30.2584796537241-23.0844223726027)* ((outagecen>43.1033561643836) * (outagecen-43.1033561643836)/((43.1033561643836--13.4738625645632)**(2/3)))**3 - 
             (43.1033561643836-23.0844223726027) * ((outagecen>30.2584796537241) * (outagecen-30.2584796537241)/((43.1033561643836--13.4738625645632)**(2/3)))**3)/(43.1033561643836-30.2584796537241);
    END;
    ELSE IF cohort=13 THEN DO;
      * UK age spline: mortality;
      * outagecen knots: -7.867397 11.6663 19.46981 25.39945 29.35892 34.36786 42.80897 ;
      outagedcanb_sp1 = ((outagecen>-7.86739726027397) * (outagecen--7.86739726027397)/((42.8089726027397--7.86739726027397)**(2/3)))**3
          + ((34.3678580493151--7.86739726027397)* ((outagecen>42.8089726027397) * (outagecen-42.8089726027397)/((42.8089726027397--7.86739726027397)**(2/3)))**3 - 
             (42.8089726027397--7.86739726027397) * ((outagecen>34.3678580493151) * (outagecen-34.3678580493151)/((42.8089726027397--7.86739726027397)**(2/3)))**3)/(42.8089726027397-34.3678580493151);
      outagedcanb_sp2 = ((outagecen>11.6662977972603) * (outagecen-11.6662977972603)/((42.8089726027397--7.86739726027397)**(2/3)))**3
          + ((34.3678580493151-11.6662977972603)* ((outagecen>42.8089726027397) * (outagecen-42.8089726027397)/((42.8089726027397--7.86739726027397)**(2/3)))**3 - 
             (42.8089726027397-11.6662977972603) * ((outagecen>34.3678580493151) * (outagecen-34.3678580493151)/((42.8089726027397--7.86739726027397)**(2/3)))**3)/(42.8089726027397-34.3678580493151);
      outagedcanb_sp3 = ((outagecen>19.4698092575342) * (outagecen-19.4698092575342)/((42.8089726027397--7.86739726027397)**(2/3)))**3
          + ((34.3678580493151-19.4698092575342)* ((outagecen>42.8089726027397) * (outagecen-42.8089726027397)/((42.8089726027397--7.86739726027397)**(2/3)))**3 - 
             (42.8089726027397-19.4698092575342) * ((outagecen>34.3678580493151) * (outagecen-34.3678580493151)/((42.8089726027397--7.86739726027397)**(2/3)))**3)/(42.8089726027397-34.3678580493151);
      outagedcanb_sp4 = ((outagecen>25.3994498091174) * (outagecen-25.3994498091174)/((42.8089726027397--7.86739726027397)**(2/3)))**3
          + ((34.3678580493151-25.3994498091174)* ((outagecen>42.8089726027397) * (outagecen-42.8089726027397)/((42.8089726027397--7.86739726027397)**(2/3)))**3 - 
             (42.8089726027397-25.3994498091174) * ((outagecen>34.3678580493151) * (outagecen-34.3678580493151)/((42.8089726027397--7.86739726027397)**(2/3)))**3)/(42.8089726027397-34.3678580493151);
      outagedcanb_sp5 = ((outagecen>29.3589186359683) * (outagecen-29.3589186359683)/((42.8089726027397--7.86739726027397)**(2/3)))**3
          + ((34.3678580493151-29.3589186359683)* ((outagecen>42.8089726027397) * (outagecen-42.8089726027397)/((42.8089726027397--7.86739726027397)**(2/3)))**3 - 
             (42.8089726027397-29.3589186359683) * ((outagecen>34.3678580493151) * (outagecen-34.3678580493151)/((42.8089726027397--7.86739726027397)**(2/3)))**3)/(42.8089726027397-34.3678580493151);
    END;
    ELSE IF cohort=14 THEN DO;
      * US age spline: mortality;
      * outagecen knots: -2.878989 13.20514 21.0795 25.72385 31.12084 36.29168 47.20753 ;
      outagedcanb_sp1 = ((outagecen>-2.8789888838985) * (outagecen--2.8789888838985)/((47.2075342465753--2.8789888838985)**(2/3)))**3
          + ((36.2916769939292--2.8789888838985)* ((outagecen>47.2075342465753) * (outagecen-47.2075342465753)/((47.2075342465753--2.8789888838985)**(2/3)))**3 - 
             (47.2075342465753--2.8789888838985) * ((outagecen>36.2916769939292) * (outagecen-36.2916769939292)/((47.2075342465753--2.8789888838985)**(2/3)))**3)/(47.2075342465753-36.2916769939292);
      outagedcanb_sp2 = ((outagecen>13.2051366152631) * (outagecen-13.2051366152631)/((47.2075342465753--2.8789888838985)**(2/3)))**3
          + ((36.2916769939292-13.2051366152631)* ((outagecen>47.2075342465753) * (outagecen-47.2075342465753)/((47.2075342465753--2.8789888838985)**(2/3)))**3 - 
             (47.2075342465753-13.2051366152631) * ((outagecen>36.2916769939292) * (outagecen-36.2916769939292)/((47.2075342465753--2.8789888838985)**(2/3)))**3)/(47.2075342465753-36.2916769939292);
      outagedcanb_sp3 = ((outagecen>21.0794977471143) * (outagecen-21.0794977471143)/((47.2075342465753--2.8789888838985)**(2/3)))**3
          + ((36.2916769939292-21.0794977471143)* ((outagecen>47.2075342465753) * (outagecen-47.2075342465753)/((47.2075342465753--2.8789888838985)**(2/3)))**3 - 
             (47.2075342465753-21.0794977471143) * ((outagecen>36.2916769939292) * (outagecen-36.2916769939292)/((47.2075342465753--2.8789888838985)**(2/3)))**3)/(47.2075342465753-36.2916769939292);
      outagedcanb_sp4 = ((outagecen>25.7238528332959) * (outagecen-25.7238528332959)/((47.2075342465753--2.8789888838985)**(2/3)))**3
          + ((36.2916769939292-25.7238528332959)* ((outagecen>47.2075342465753) * (outagecen-47.2075342465753)/((47.2075342465753--2.8789888838985)**(2/3)))**3 - 
             (47.2075342465753-25.7238528332959) * ((outagecen>36.2916769939292) * (outagecen-36.2916769939292)/((47.2075342465753--2.8789888838985)**(2/3)))**3)/(47.2075342465753-36.2916769939292);
      outagedcanb_sp5 = ((outagecen>31.1208439945205) * (outagecen-31.1208439945205)/((47.2075342465753--2.8789888838985)**(2/3)))**3
          + ((36.2916769939292-31.1208439945205)* ((outagecen>47.2075342465753) * (outagecen-47.2075342465753)/((47.2075342465753--2.8789888838985)**(2/3)))**3 - 
             (47.2075342465753-31.1208439945205) * ((outagecen>36.2916769939292) * (outagecen-36.2916769939292)/((47.2075342465753--2.8789888838985)**(2/3)))**3)/(47.2075342465753-36.2916769939292);
    END;        
%MEND;

%MACRO TVlagworkTRANSFORM();
  *transformations of time varying lagged employment variables;          
    * wyr_cum_lag1 knots: 1.612154 14.90274 35.5889 ;	
    wyr_cum_lag1d_sp1 = ((wyr_cum_lag1>1.61215435287072) * (wyr_cum_lag1-1.61215435287072)/((35.588904109589-1.61215435287072)**(2/3)))**3
        + ((14.9027397260274-1.61215435287072)* ((wyr_cum_lag1>35.588904109589) * (wyr_cum_lag1-35.588904109589)/((35.588904109589-1.61215435287072)**(2/3)))**3 - 
           (35.588904109589-1.61215435287072) * ((wyr_cum_lag1>14.9027397260274) * (wyr_cum_lag1-14.9027397260274)/((35.588904109589-1.61215435287072)**(2/3)))**3)/(35.588904109589-14.9027397260274);
    * wyr_cum_lag1 knots: 1.292059 14.22298 32.52767 ;
    wyr_cum_lag1aw_sp1 = ((wyr_cum_lag1>1.29205891159518) * (wyr_cum_lag1-1.29205891159518)/((32.5276712328767-1.29205891159518)**(2/3)))**3
        + ((14.2229807620331-1.29205891159518)* ((wyr_cum_lag1>32.5276712328767) * (wyr_cum_lag1-32.5276712328767)/((32.5276712328767-1.29205891159518)**(2/3)))**3 - 
           (32.5276712328767-1.29205891159518) * ((wyr_cum_lag1>14.2229807620331) * (wyr_cum_lag1-14.2229807620331)/((32.5276712328767-1.29205891159518)**(2/3)))**3)/(32.5276712328767-14.2229807620331);
%MEND;

%MACRO TVworkTRANSFORM();
  *transformations of time varying current employment variables (not actually needed);          
%MEND;

%MACRO TVlagdoseTRANSFORM();
  *transformations of time varying lagged exposure variables;   
  
    l1dose_lag1 = log(dose_lag1 +      0.01);
    l1dose_lag2 = log(dose_lag2 +      0.01);
    l1dose_lag3 = log(dose_lag3 +      0.01);
    l1dose_lag4 = log(dose_lag4 +      0.01);
    l1dose_lag5 = log(dose_lag5 +      0.01);
    l1dose_lag6 = log(dose_lag6 +      0.01);
    l1dose_lag7 = log(dose_lag7 +      0.01);
    l1dose_lag8 = log(dose_lag8 +      0.01);
    l1dose_lag9 = log(dose_lag9 +      0.01);
    l1dose_lag10 = log(dose_lag10 +    0.01);
    l1dose_lag15 = log(dose_lag15 +    0.01);
    l1dose_lag20 = log(dose_lag20 +    0.01);
    l1dose_lag25 = log(dose_lag25 +    0.01);  
    cl1dosel_1_5 = log(dosel_1_5 +     0.01);
    cl1dosel_2_5 = log(dosel_2_5 +     0.01);
    cl1dosel_6_10 = log(dosel_6_10 +   0.01);
    cl1dosel_11_15 = log(dosel_11_15 + 0.01);
    cl1dose_lag5 = log(ds_cu_lag5 +    0.01);
    cl1dose_lag10 = log(ds_cu_lag10 +  0.01);
    cl1dose_lag15 = log(ds_cu_lag15 +  0.01);
    
    **** country specific dose splines: 1-5, 6 year lagged cumulative exposure;
    /* note: these are modified splines that only use basis functions corresponding to unique knots of a 7df spline */
    /* the default quantile based approach resulted in multiple knots at zero, so these were excluded */

    **** country specific dose splines: 1-5, 6 year lagged cumulative exposure;
    IF cohort=5 THEN DO;
      * French dose spline: 1-5, 6 year lagged cumulative exposure;
      * ds_cu_lag6 knots: 6.94023 24.7091 149.9005 ;
      ds_cu_lag6wb_sp1 = ((ds_cu_lag6>6.94022969) * (ds_cu_lag6-6.94022969)/((149.900525-6.94022969)**(2/3)))**3
          + ((24.709100585-6.94022969)* ((ds_cu_lag6>149.900525) * (ds_cu_lag6-149.900525)/((149.900525-6.94022969)**(2/3)))**3 - 
             (149.900525-6.94022969) * ((ds_cu_lag6>24.709100585) * (ds_cu_lag6-24.709100585)/((149.900525-6.94022969)**(2/3)))**3)/(149.900525-24.709100585);
    END;
    ELSE IF cohort=13 THEN DO;
      * UK dose spline: 1-5, 6 year lagged cumulative exposure;
      * ds_cu_lag6 knots: 2.3764 15.34475 198.2462 ;
      ds_cu_lag6wb_sp1 = ((ds_cu_lag6>2.37640003) * (ds_cu_lag6-2.37640003)/((198.24625-2.37640003)**(2/3)))**3
          + ((15.34474824-2.37640003)* ((ds_cu_lag6>198.24625) * (ds_cu_lag6-198.24625)/((198.24625-2.37640003)**(2/3)))**3 - 
             (198.24625-2.37640003) * ((ds_cu_lag6>15.34474824) * (ds_cu_lag6-15.34474824)/((198.24625-2.37640003)**(2/3)))**3)/(198.24625-15.34474824);
    END;
    ELSE IF cohort=14 THEN DO;
      * US dose spline: 1-5, 6 year lagged cumulative exposure;
      * ds_cu_lag6 knots: 4.165197 19.17678 142.2469 ;
      ds_cu_lag6wb_sp1 = ((ds_cu_lag6>4.165197328) * (ds_cu_lag6-4.165197328)/((142.24695-4.165197328)**(2/3)))**3
          + ((19.176777688-4.165197328)* ((ds_cu_lag6>142.24695) * (ds_cu_lag6-142.24695)/((142.24695-4.165197328)**(2/3)))**3 - 
             (142.24695-4.165197328) * ((ds_cu_lag6>19.176777688) * (ds_cu_lag6-19.176777688)/((142.24695-4.165197328)**(2/3)))**3)/(142.24695-19.176777688);
    END;   
%MEND;

%MACRO TVdoseTRANSFORM();
  *transformations of time varying current exposure variables (not actually used);          
%MEND;

%MACRO deterministicINTERVENTIONS();
* INTERVENTIONS with a negative value are determined as deterministic interventions 
 where all employed person years are set to negative the value of the intervention.;
       IF intervention < 0 THEN DO;
         IF intervention=-999 THEN DO;
           * never exposed;
           logdosert = .N;
	       anydose = 0; *need this to avoid recalculating;
           dose = 0;
         END;
         ELSE DO;
	       anydose = 1; 
           loglimit = LOG(-intervention);
           logdosert = loglimit;
         END;         
      END;
%MEND;


%MACRO bound(var, upperbound=700, lowerbound=-700);
       /* prevent numerical under/overflow */
       IF .z < &var < &lowerbound THEN &var=&lowerbound;
       ELSE IF &var > &upperbound  THEN &var=&upperbound;
%MEND;


%MACRO AssignMinMaxDose(ds);
  /* create set of macros that assigns min and max ln-exposure rate by cohort/decade */
  /* requires input dataset with variables: cohort indecade max_logdosert min_logdosert */
  %GLOBAL cohortlist;
  %GLOBAL decadelist;
  %GLOBAL maxratelist;
  %GLOBAL minratelist;
  %LET cohort=;
  %LET decade=;
  %LET maxrate=;
  %LET minrate=;
  %LET combos=;
  DATA &ds;
   SET &ds;
   CALL SYMPUT("combos", _N_); 
  RUN;
  %DO __ci = 1 %TO &combos;
    DATA tmp;
     SET &ds(OBS=&__ci FIRSTOBS=&__ci);
     CALL SYMPUT("cohort", TRIM(PUT(cohort, BEST9.)));
     CALL SYMPUT("decade", TRIM(PUT(indecade, BEST9.)));
     CALL SYMPUT("maxrate", TRIM(PUT(max_logdosert, BEST9.)));
     CALL SYMPUT("minrate", TRIM(PUT(min_logdosert, BEST9.)));
    RUN;
    %LET cohortlist=&cohortlist &cohort;
    %LET decadelist=&decadelist &decade;
    %LET maxratelist=&maxratelist &maxrate;
    %LET minratelist=&minratelist &minrate;
  %END; 
  %PUT &combos;
  %PUT &cohortlist &decadelist &maxratelist &minratelist;
%MEND;

%MACRO SetMinMaxDose(clist=&cohortlist, dlist=&decadelist, maxlist=&maxratelist, minlist=&minratelist);
  /* create code that uses min/max exposure rate by decade/cohort */
  %LET maxscan = 1;
  IF 0 THEN DO; END;
  %DO %WHILE(%QSCAN(&clist, %EVAL(&maxscan), " ")^=);
    ELSE IF cohort = %QSCAN(&cohortlist, %EVAL(&maxscan), " ") AND indecade =  %QSCAN(&dlist, %EVAL(&maxscan), " ") THEN DO;
      max_logdosert = %QSCAN(&maxlist, %EVAL(&maxscan), " ");
      min_logdosert = %QSCAN(&minlist, %EVAL(&maxscan), " ");
    END;
    %LET maxscan = %EVAL(&MAXSCAN+1);
  %END;
%MEND;