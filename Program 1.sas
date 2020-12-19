








* logrank test for two or more groups;
proc import datafile="/folders/myfolders/Project1/all_tert_buc.csv" 
out=mydata_buc dbms=csv replace;
datarow=2;
getnames=yes;
guessingrows=100;
run;



* Global test; 
proc lifetest data=mydata_buc plots=survival;* intervals=(0 to 24 by 1);* (atrisk=0 to 24 by 1);*(atrisk=0 to 2500 by 500);  
time Overall_Survival__Months_*censor(1); 
strata TERT_mutation/test=(logrank TARONE PETO MODPETO  
FLEMING(0,1)    ); 
run;


proc phreg data = mydata_buc plots=survival;
  class TERT_mutation;
  model Overall_Survival__Months_*censor(1)=TERT_mutation;
run;




ods text="TESTING full model phreg stepwise";
proc phreg data = mydata_buc;
  class TERT_mutation Sex Smoking_History Sample_Class Sample_Type Specimen_Type Specimen_Preservation_Type Matched_Status;
  model Overall_Survival__Months_*censor(1)= Sex Smoking_History TERT_mutation Tumor_Purity Fraction_Genome_Altered Mutation_Count DNA_Input Sample_Class Sample_Type Specimen_Type Specimen_Preservation_Type Matched_Status/ selection=stepwise slentry=0.99
                           slstay=0.99 details;
  *sex: test z3 = z4 = z3xz4 = 0;
  *ods output  FitStatistics = sexaic_t6;
  *ods output TestStmts = sexstat_t6;
run;

ods text="TESTING small model phreg stepwise";
proc phreg data = mydata_buc;
  class TERT_mutation;
  model Overall_Survival__Months_*censor(1)= TERT_mutation / selection=stepwise slentry=0.4
                           slstay=0.25 details;
  *sex: test z3 = z4 = z3xz4 = 0;
  *ods output  FitStatistics = sexaic_t6;
  *ods output TestStmts = sexstat_t6;
run;




proc phreg data = mydata_buc plots=survival;
  class TERT_mutation Sex Smoking_History Sample_Class Sample_Type Specimen_Type Specimen_Preservation_Type Matched_Status;
  model Overall_Survival__Months_*censor(1)= Sex Smoking_History TERT_mutation Tumor_Purity Fraction_Genome_Altered Mutation_Count DNA_Input Sample_Class Sample_Type Specimen_Type Specimen_Preservation_Type Matched_Status;
run;














* A Cox regression model is fit to the data with these eight covariates
* residuals are computed;
* the default method of estimating a survivor function is Breslow (1972) estimator, i.e..method=ch;
proc phreg data = mydata_buc;
  class TERT_mutation Sex Smoking_History Sample_Class Sample_Type Specimen_Type Specimen_Preservation_Type Matched_Status;
  model Overall_Survival__Months_*censor(1)= Sex Smoking_History TERT_mutation Tumor_Purity Fraction_Genome_Altered Mutation_Count DNA_Input Sample_Class Sample_Type Specimen_Type Specimen_Preservation_Type Matched_Status;	
  output out = plot1_1 LOGSURV = logsurv1 /method = ch; /*-logsurv is the cox-snell residual*/
run;
data plot1_1;
set plot1_1;
snell = -logsurv1;
cons = 1;
run;


proc phreg data = plot1_1;
model snell*censor(1) = cons;
output out = plot1_2 logsurv = logsurv2 /method = ch;
run;

data plot1_2;
set plot1_2;
cumhaz = - logsurv2;
run;

proc sort data = plot1_2;
by snell;
run;


proc sgplot data = plot1_2;
step y=cumhaz x=snell /MARKERFILLATTRS=(color="red");
lineparm x=0 y=0 slope=1; /** intercept, slope **/
label cumhaz = "Estimated Cumulative Hazard Rates";
label snell = "Residual";
run;





* A Cox regression model is fit to the data with these eight covariates
* residuals are computed;
* the default method of estimating a survivor function is Breslow (1972) estimator, i.e..method=ch;
proc phreg data = mydata_buc;
  class TERT_mutation Sample_Type;
  model Overall_Survival__Months_*censor(1)= TERT_mutation Fraction_Genome_Altered Mutation_Count Sample_Type;	
  output out = plot1_1 LOGSURV = logsurv1 /method = ch; /*-logsurv is the cox-snell residual*/
run;
data plot1_1;
set plot1_1;
snell = -logsurv1;
cons = 1;
run;


proc phreg data = plot1_1;
model snell*censor(1) = cons;
output out = plot1_2 logsurv = logsurv2 /method = ch;
run;

data plot1_2;
set plot1_2;
cumhaz = - logsurv2;
run;

proc sort data = plot1_2;
by snell;
run;


proc sgplot data = plot1_2;
step y=cumhaz x=snell /MARKERFILLATTRS=(color="red");
lineparm x=0 y=0 slope=1; /** intercept, slope **/
label cumhaz = "Estimated Cumulative Hazard Rates";
label snell = "Residual";
run;




* Martingale residual;
* leave out patient age;
proc phreg data = mydata_buc;
  class TERT_mutation Sample_Type;
  model Overall_Survival__Months_*censor(1)= TERT_mutation Mutation_Count Sample_Type;	

output out = plot2_1 RESMART = mgale ;
run;
ods listing gpath='/folders/myfolders/Lab8/';
ods graphics / imagename="p2" imagefmt=png;
proc loess data=plot2_1;
model mgale = Fraction_Genome_Altered / smooth=0.6 direct;
run;



* Martingale residual;
* leave out patient age;
proc phreg data = mydata_buc;
  class TERT_mutation Sample_Type;
  model Overall_Survival__Months_*censor(1)= TERT_mutation Fraction_Genome_Altered Sample_Type;	

output out = plot2_1 RESMART = mgale ;
run;
ods listing gpath='/folders/myfolders/Lab8/';
ods graphics / imagename="p2" imagefmt=png;
proc loess data=plot2_1;
model mgale = Mutation_Count / smooth=0.6 direct;
run;






data mydata_buc;
set mydata_buc;
cons = 1;
run;



* The baseline cumulative hazards are estimated using
* Breslow’s estimator for each stratu;
proc phreg data = mydata_buc ;
  class TERT_mutation;
  model Overall_Survival__Months_*censor(1)=  cons/rl;	
strata TERT_mutation;
output out = base logsurv = ls /method = ch;
run;
data base;
set base;
logH = log (-ls);
if TERT_mutation= 0 then logH1 = logH;
else if TERT_mutation= 1 then logH2 = logH;
proc sort;by Overall_Survival__Months_ TERT_mutation;
proc print;var Overall_Survival__Months_ TERT_mutation logH logH1 logH2;
run;
ods listing gpath='/folders/myfolders/Lab8/';
ods graphics / imagename="p3" imagefmt=png;
proc sgplot data = base;
where logH ne .;
series x=Overall_Survival__Months_ y=logH /group=TERT_mutation ;
run;






proc sort data = base;
by Overall_Survival__Months_;
run;
data base;
set base;
retain temp1 temp2 temp3;
if logH1 ~= . then temp1 = logH1;
if logH2 ~= . then temp2 = logH2;
diff2v1 = temp2 - temp1;
proc print;
var TERT_mutation Overall_Survival__Months_  diff2v1;
run;
ods listing gpath='/folders/myfolders/Lab8/';
ods graphics / imagename="p4" imagefmt=png;
proc sgplot data =base;
where Overall_Survival__Months_ <=24;
step x=Overall_Survival__Months_ y=diff2v1 ;
lineparm x=0 y=0 slope=0; /** intercept, slope **/
run;