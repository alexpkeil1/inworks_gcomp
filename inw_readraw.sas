/**********************************************************************************************************************
* Author: Alex Keil
* Project: INWORKS (International Nuclear Workers Study)
* Description: Read in raw data files and create person-level data
* Keywords: radiation, occupational, g-formula
* Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
**********************************************************************************************************************/


OPTIONS nofmterr; 

DATA cohort; 
  SET &raww;
RUN;


DATA exposure; 
  SET &rawx; 
RUN;



DATA deaths;
  SET &rawy;
RUN;
 
PROC SORT DATA=cohort;  BY WKERID;
PROC SORT DATA=exposure;BY WKERID;
PROC SORT DATA=deaths;  BY WKERID; 
RUN;

DATA &indata;
 MERGE cohort(IN=ina) exposure deaths; 
 BY WKERID;
 IF ina;
RUN;QUIT;RUN;

