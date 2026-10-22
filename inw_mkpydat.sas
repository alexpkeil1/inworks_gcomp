/**********************************************************************************************************************
* Author: Alex Keil
* Project: INWORKS (International Nuclear Workers Study)
* Description: Make person-time data
* Keywords: radiation, occupational, g-formula
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/


PROC SQL NOPRINT; SELECT COUNT(1) INTO :N FROM &indata;QUIT;; * get sample size;


%MACRO pickdata();
* choose between bootstrap sample or original data ;
  %IF &dobootstrap = 0 %THEN %DO;
    DATA an;
     SET &INDATA (obs=&n);
    RUN;
  %END;
  
  %IF &dobootstrap = 1 %THEN %DO;
    PROC SURVEYSELECT DATA=&INDATA OUT=an METHOD=URS OUTHITS SAMPSIZE=&N;RUN;
    DATA an;
      SET an;
	  wkerid_old = wkerid;
      wkerid = _N_;
  %END;
%MEND;

%PICKDATA();

%MACRO mkpydata();
* create person-period data from person-level data ;

DATA anpy ;
  SET an;
  /* containers */
  ARRAY hp[*] hp1932-hp2006;
  ARRAY outcomes[*] cancer solid leuk can_nole sol_nolu oral eso stom colon rectum liver panc perit larynx lung breast 
                    uterus ovary prost kidney bladder brain thyr conn gallbldr livernog bone pleura allskin testis breastf pleuraex;
  /* time fixed indicators */
  usa = (cohort=14);
  uk = (cohort=13);
  frn = (cohort=5);
  intos = 0;
  outtos = 0;
  
  ses1 = (ses= 1 );
  ses2 = (ses= 2 );
  ses3 = (ses= 3 );
  ses4 = (ses= 4 );
  ses9 = (ses= 9 );
  
  agelastwork = YRDIF(dob,end_emp_, 'ACT/ACT');
  /* create containers */
  wyr_cum=0;
  dose_cum=0;
  %vsINITS(outcomes=d_cancer d_solid d_leuk d_can_nole d_sol_nolu d_oral d_eso d_stom d_colon d_rectum d_liver d_panc d_perit d_larynx d_lung d_breast d_uterus 
                      d_ovary d_prost d_kidney d_bladder d_brain d_thyr d_conn d_gallbldr d_livernog d_bone d_pleura d_allskin d_testis d_breastf d_pleuraex
                      atrisk d_death d_noncancer);
  * create time varying exposure/work variables (note that these are not actually calculated here);
  %tvINITS(pref=,size=26, setzero=1);

  /* time varying */
  DO _inyearmin = 1932 to 2020;
    /* time variables */
    _indatemin = MDY(1,1,_inyearmin); 
    _outdatemax = MIN(dlo, MDY(12,31,_inyearmin));
    IF dlo < _indatemin 
      THEN LEAVE;
    IF bfu < _indatemin OR bfu > _outdatemax THEN DO;
      indate = _indatemin;
    END;
    ELSE DO;
      indate = MAX(bfu, _indatemin);
 	 atrisk = 1;
    END;
    IF dlo > (_outdatemax) THEN DO;
      outdate = _outdatemax;
    END;
    ELSE DO;
      outdate = MIN(dlo, _outdatemax);
    END;
    
    * decimal year of birth;
    yob = 1900+YRDIF(MDY(1,1,1900), dob, 'ACT/ACT');

	* age at end of follow up;
	eofage = YRDIF(dob, MDY(8, 1, 2006), 'ACT/ACT');

    *age;
    inage = YRDIF(dob,indate, 'ACT/ACT');
    outage = YRDIF(dob,outdate, 'ACT/ACT');
    *calendar year;
    _yrindays = FLOOR(indate-MDY(1,1,YEAR(indate)));
    _yroutdays = FLOOR(outdate-MDY(1,1,YEAR(indate))+1);
    _totyrdays = MDY(12,31,YEAR(indate))-MDY(1,1,YEAR(indate))+1;
    inyear = YEAR(indate) + _yrindays/_totyrdays;
    outyear = YEAR(indate) + _yroutdays/_totyrdays;
 
    totfu = YRDIF(bfu, outdate+1, 'ACT/ACT');
    /* at risk time */
    pyrs = atrisk*(outyear - inyear);
    *time on study;
    outtos = intos+pyrs;
       
    /* employment time */
    wyr =  MAX(0, MIN(
                   MIN(1, YRDIF(indate, end_emp_, 'ACT/ACT')),
                   MAX(0, YRDIF(start_em, outdate, 'ACT/ACT'))
                   ));               
    timesincehire = YRDIF(start_em, outdate, 'ACT/ACT');
    * last date of employment is 12/31/2004, last DLO is 2006-08-27;
    IF _inyearmin <= 2004 THEN dose = MAX(0,hp[_inyearmin-1931]); *hp1932-hp2006;
    ELSE dose = 0;
    /*admin censoring dates? for workers by country*/
    facop = 1;
    IF mainfacy IN (11,12) THEN DO;
      IF inyear >= 2004.9972678 THEN facop = 0;
    END;
    ELSE IF mainfacy IN (13) THEN DO;
      IF inyear >= 2003.9972678 THEN facop = 0;
    END;
    ELSE IF mainfacy IN (21, 23, 24, 25, 221, 222) THEN DO;
      IF inyear >= 2002 THEN facop = 0;
    END;
    ELSE IF mainfacy IN (454, 460) THEN DO;
      IF inyear >= 1998.9972678 THEN facop = 0;
    END;
    ELSE IF mainfacy IN (456) THEN DO;
      IF inyear >= 1991.9972678 THEN facop = 0;
    END;
    ELSE IF mainfacy IN (457) THEN DO;
      IF inyear >= 1988.1420765027 THEN facop = 0;
    END;
    ELSE IF mainfacy IN (458) THEN DO;
      IF inyear >= 1991.4958904 THEN facop = 0;
    END;
    ELSE IF inyear >= 2004.9972678 THEN facop = 0; /* in case there are other facilities */
    

    /* cumulative variables */
	* LOG DOSE RATE CREATED HERE: log(dose/wyr);
      %CUMdoseCREATE(inline=,prefa=, prefb=, size=26);
      %CUMworkCREATE(inline=,prefa=, prefb=, size=26);
    %dosewindowCREATE(pref=, size=26, agest=16, ageend=85, agepts = &ageendpoints, lagpts= &lagendpoints, id=wkerid);
    /* vital status */
    IF dlo = outdate THEN DO;
      d_death = vs;
 	  DO j = 1 TO DIM(outcomes);
        d_outcomes[j] = outcomes[j];
 	  END;
    END;
    * competing events for cancer incidence (non cancer related deaths);
    d_nonlung = (1-d_lung)*d_death;
    d_noncancer =  (1-d_cancer)*d_death;
    d_nonsolid =  (1-d_solid)*d_death;
    * fix for single day risk periods;
    IF inage = outage THEN DO;
      outage = outage + 1/367;
      outyear = outyear + 1/367;
      outtos = outtos + 1/367;
    END;
    /* output */
    IF pyrs>0 THEN OUTPUT anpy;
      /* lags */
      intos=outtos;
      %LAGS(size=26);
   FORMAT indate outdate _indatemin _outdatemax MMDDYY10.;
   DROP _: hp: i j;
  END;
RUN;


%MEND;

%MKPYDATA();


RUN;QUIT;RUN;

