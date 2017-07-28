
%macro nlm(data);
proc nlmixed data=&data qpoints=10;
parms beta1=0.1 beta2=0.1 beta7=0.1  beta12=0.1 sigma1=1.0 sigma2=1.2 ro_sigma=0.5 sigmae=1;
bounds sigma1 sigma2 sigmae>0;
temp0=X1*beta1+X2*beta2+u1;
temp1=R1*temp0;
temp2=log(1+exp(temp0));
temp3=R1/2*log(sigmae);
temp4=R1*((log(Y1)-X1*beta7-X6*beta12-u2)**2)/2/sigmae;
loglik=temp1-temp2-temp3-temp4;
model Y1~ general(loglik);
random u1 u2 ~normal([0,0],[sigma1,ro_sigma,sigma2]) subject=in;
ods output ParameterEstimates=two.test1 ;
run;
%mend();
%macro estimate();
libname two 'E:\two-part\code\sas code\oracle';
%do jj=1%to 1;
proc import out=two.two_part2&jj 
datafile="E:\two-part\code\sas code\oracle\data\two_part2&jj..csv" 
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
data two.estoracleall;  
set _null_;
run;
%do ii=1%to 1;
proc transpose data=two.one&ii out=two.onet&ii;
run;
data two.estoracleall;   
set two.estoracleall two.onet&ii;   
run;
%end;
%mend();
%toge()

%macro estimate1();
proc export data=two.estoracleall
outfile="E:\two-part\code\sas code\oracle\Results_save1\estoracleall.csv"
dbms=csv
label
replace;
run;
%mend();
%estimate1();

