/**********************************************************************************************************************
* Author: Alex Keil
* Project: INWORKS (International Nuclear Workers Study)
* Description: A central program that runs all other programs to create and clean analytic datasets for use in g-formula analysis
* Keywords: radiation, occupational, g-formula
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/


************ input file locations / names ********** ;
* data directories;
  %LET prgdir = <redacted>\prg;
  %LET rawdir = <redacted>\rawdata;
  %LET intdir = <redacted>\dataprocessed;
  %LET outdir = <redacted>\output;

* names of raw input files (XPT format) (these should be in the "rawdir" directory);
  %LET raww = demo.xpt;
  %LET rawx = exposures.xpt;
  %LET rawy = outcomes.xpt;


************ global MACRO variables ********** ;
* complementary outcomes - cause1 is outcome of interest, cause2 is composite of all competing risks;
  %LET cause1 = cancer;    * incident event coded as d_&cause1, CI coded as ci_&cause1;
  %LET cause2 = noncancer; * incident event coded as d_&cause2, CI coded as ci_&cause2;


* number of MC samples to take in g-formula simulation;
  %LET nmc_main  = 80000; 
  %LET nmc_boot  = 80000; 

* number of bootstrap iterations;
  %LET bootstart = 1;
  %LET bootiters = 200;

* age/time-since exposure window endpoints;
  %LET ageendpoints = 0 30 45 100; * lowest is minimum age at which exposure accumulates;
  %LET lagendpoints = 10 20 30 100;* lowest is minimum lag at which exposure accumulates;

* annual exposure limits (in mSv) to compare (will also include natural course);
* negative values are deterministic values of assigned exposure (if employed);
* -999 is setting exposure to zero;
  %LET INTS = %STR(-999, -1.0, -2.0, -5.0, -20.0);
  
* person-level data from inworks (merged raw input files);
  %LET INDATA = work.inworks_merged;


********* end global MACRO variaables ********* ;

************ global function definitions ********** ;
PROC FCMP OUTLIB=work.misc.addyrs;
 FUNCTION addyrs(startdate, years);
  floor = INTNX('YEAR', startdate, years, "sameday");
  days = (years-FLOOR(years))*(INTNX('YEAR', floor, 1, 'sameday')-floor);
  enddate = floor+days;
 RETURN(enddate);
 ENDSUB;
RUN;
OPTIONS CMPLIB=work.misc;
********* end global function definitions ********* ;

LIBNAME raw "&rawdir";
LIBNAME dat "&intdir";
LIBNAME gformula "&outdir";

* ALL INPUT DATA HERE;
LIBNAME in1 XPORT "&rawdir\&raww";
LIBNAME in2 XPORT "&rawdir\&rawx";
LIBNAME in3 XPORT "&rawdir\&rawy";


* main result;
%LET dobootstrap=0;
%LET nmc  = &nmc_main; 

%INCLUDE "&prgdir/MACROS.sas";
%INCLUDE "&prgdir/inw_readraw.sas";
%INCLUDE "&prgdir/inw_mkpydat.sas";
%INCLUDE "&prgdir/inw_gformula.sas";
%INCLUDE "&prgdir/inw_showmodels.sas";

%LET bootiter=;
%INCLUDE "&prgdir/inw_survcurve.sas";

* bootstrap confidence intervals;
%MACRO bootstrap();
OPTIONS NONOTES NOMPRINT NOSYMBOLGEN NOSOURCE; ODS LISTING CLOSE;
  %IF &bootiters NE 0 %THEN %DO;
    %LET nmc  = &nmc_boot; 
    %LET dobootstrap = 1;
    %DO bootiter = &bootstart %TO &bootiters;
      %INCLUDE "&prgdir/inw_mkpydat.sas";
      %INCLUDE "&prgdir/inw_gformula.sas";
      %INCLUDE "&prgdir/inw_survcurve.sas";
    %END;
  %END;  
OPTIONS NOTES SOURCE; ODS LISTING;
%MEND;

%BOOTSTRAP();

RUN;QUIT;RUN;
/*DM ODSRESULTS 'clear;' CONTINUE; *clear ODS generated datasets;*/

PROC PRINTTO ;
RUN; QUIT; RUN;
