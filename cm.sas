
* logrank test for two or more groups;
proc import datafile="/folders/myfolders/Project1/all_tert_cm.csv" 
out=mydata dbms=csv replace;
datarow=2;
getnames=yes;
guessingrows=100;
run;



* Global test; 
proc lifetest data=mydata plots=survival;*(atrisk=0 to 2500 by 500);  
time Overall_Survival__Months_*censor(1); 
strata TERT_mutation/test=(logrank TARONE PETO MODPETO  
FLEMING(0,1)    ); 
run;


proc phreg data = mydata plots=survival;
  class TERT_mutation;
  model Overall_Survival__Months_*censor(1)=TERT_mutation;
run;

