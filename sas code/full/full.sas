
%macro nlm(data);
proc nlmixed data=&data qpoints=10;
parms beta1=0.1 beta2=0.1 beta3=0.1 beta4=0.1 beta5=0.1 beta6=0.1 beta7=1 beta8=0.1 beta9=0.1 
beta10=0.1 beta11=0.1 beta12=0.1 sigma1=1.0 sigma2=1.2 ro_sigma=0.5 sigmae=1;
bounds sigma1 sigma2 sigmae>0;
temp0=X1*beta1+X2*beta2+X3*beta3+X4*beta4+X5*beta5+X6*beta6+u1;
temp1=R1*temp0;
temp2=log(1+exp(temp0));
temp3=R1/2*log(sigmae);
temp4=R1*((log(Y1)-X1*beta7-X2*beta8-X3*beta9-X4*beta10-X5*beta11-X6*beta12-u2)**2)/2/sigmae;
loglik=temp1-temp2-temp3-temp4;
model Y1~ general(loglik);
random u1 u2 ~normal([0,0],[sigma1,ro_sigma,sigma2]) subject=in;
ods output ParameterEstimates=two.test1 ;
run;
%mend();
%macro estimate();
libname two 'E:\two-part\code\sas code\full';
%do jj=1%to 1;
proc import out=two.two_part2&jj 
datafile="E:\two-part\code\sas code\full\data\two_part2&jj..csv" 
dbms=csv
replace;
run;
%nlm(two.two_part2&jj) 
data two.one;
set _null_;
run;
data two.one;
set two.one two.test1;
KEEP parameter Estimate;
run;
proc transpose data=two.one out=two.one&jj;
id parameter;
run;
%end;
%mend();
%estimate()
%macro toge();
data two.estfullall;  
set _null_;
run;
%do ii=1%to 1;
proc transpose data=two.one&ii out=two.onet&ii;
run;
data two.estfullall;   
set two.estfullall two.onet&ii;   
run;
%end;
%mend();
%toge()

%macro estimate1();
proc export data=two.estfullall
outfile="E:\two-part\code\sas code\full\Results_save1\estfullall.csv"
dbms=csv
label
replace;
run;
%mend();
%estimate1();

