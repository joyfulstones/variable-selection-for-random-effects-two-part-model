
%macro nlm0(data);
proc nlmixed data=&data qpoints=10 hess;
parms gamma1=0.1 gamma2=0.1 gamma3=0.1 gamma4=0.1 gamma5=0.1 gamma6=0.1 gamma7=1 gamma8=0.1 gamma9=0.1 
gamma10=0.1 gamma11=0.1 gamma12=0.1 sigma1=1.0 sigma2=1.2 ro_sigma=0.5 sigmae=1;
bounds sigma1 sigma2 sigmae>0;
omega1=(exp(2*100/4*gamma1*gamma1)-1)/(exp(2*100/4*gamma1*gamma1)+1)*one;
omega2=(exp(2*100/4*gamma2*gamma2)-1)/(exp(2*100/4*gamma2*gamma2)+1)*one;
omega3=(exp(2*100/4*gamma3*gamma3)-1)/(exp(2*100/4*gamma3*gamma3)+1)*one;
omega4=(exp(2*100/4*gamma4*gamma4)-1)/(exp(2*100/4*gamma4*gamma4)+1)*one;
omega5=(exp(2*100/4*gamma5*gamma5)-1)/(exp(2*100/4*gamma5*gamma5)+1)*one;
omega6=(exp(2*100/4*gamma6*gamma6)-1)/(exp(2*100/4*gamma6*gamma6)+1)*one;
omega7=(exp(2*100/4*gamma7*gamma7)-1)/(exp(2*100/4*gamma7*gamma7)+1)*one;
omega8=(exp(2*100/4*gamma8*gamma8)-1)/(exp(2*100/4*gamma8*gamma8)+1)*one;
omega9=(exp(2*100/4*gamma9*gamma9)-1)/(exp(2*100/4*gamma9*gamma9)+1)*one;
omega10=(exp(2*100/4*gamma10*gamma10)-1)/(exp(2*100/4*gamma10*gamma10)+1)*one;
omega11=(exp(2*100/4*gamma11*gamma11)-1)/(exp(2*100/4*gamma11*gamma11)+1)*one;
omega12=(exp(2*100/4*gamma12*gamma12)-1)/(exp(2*100/4*gamma12*gamma12)+1)*one;
beta1=omega1*gamma1;
beta2=omega2*gamma2;
beta3=omega3*gamma3;
beta4=omega4*gamma4;
beta5=omega5*gamma5;
beta6=omega6*gamma6;
beta7=omega7*gamma7;
beta8=omega8*gamma8;
beta9=omega9*gamma9;
beta10=omega10*gamma10;
beta11=omega11*gamma11;
beta12=omega12*gamma12;
temp0=X1*beta1+X2*beta2+X3*beta3+X4*beta4+X5*beta5+X6*beta6+u1;
temp1=R1*temp0;
temp2=log(1+exp(temp0));
temp3=R1/2*log(sigmae);
temp4=R1*((log(Y1)-X1*beta7-X2*beta8-X3*beta9-X4*beta10-X5*beta11-X6*beta12-u2)**2)/2/sigmae;
temp5=(omega1+omega2+omega3+omega4+omega5+omega6+omega7+omega8+omega9+omega10+omega11+omega12)/5/100*(0/4*log(500)/2+(4-0)/4*log(100)/2); 
loglik=temp1-temp2-temp3-temp4-temp5;
model Y1~ general(loglik);
random u1 u2 ~normal([0,0],[sigma1,ro_sigma,sigma2]) subject=in;
ods output ParameterEstimates=two.test1 Hessian=two.cov FitStatistics=two.test2;
run;
%mend();
%macro nlm1(data);
proc nlmixed data=&data qpoints=10 hess;
parms gamma1=0.1 gamma2=0.1 gamma3=0.1 gamma4=0.1 gamma5=0.1 gamma6=0.1 gamma7=1 gamma8=0.1 gamma9=0.1 
gamma10=0.1 gamma11=0.1 gamma12=0.1 sigma1=1.0 sigma2=1.2 ro_sigma=0.5 sigmae=1;
bounds sigma1 sigma2 sigmae>0;
omega1=(exp(2*100/4*gamma1*gamma1)-1)/(exp(2*100/4*gamma1*gamma1)+1)*one;
omega2=(exp(2*100/4*gamma2*gamma2)-1)/(exp(2*100/4*gamma2*gamma2)+1)*one;
omega3=(exp(2*100/4*gamma3*gamma3)-1)/(exp(2*100/4*gamma3*gamma3)+1)*one;
omega4=(exp(2*100/4*gamma4*gamma4)-1)/(exp(2*100/4*gamma4*gamma4)+1)*one;
omega5=(exp(2*100/4*gamma5*gamma5)-1)/(exp(2*100/4*gamma5*gamma5)+1)*one;
omega6=(exp(2*100/4*gamma6*gamma6)-1)/(exp(2*100/4*gamma6*gamma6)+1)*one;
omega7=(exp(2*100/4*gamma7*gamma7)-1)/(exp(2*100/4*gamma7*gamma7)+1)*one;
omega8=(exp(2*100/4*gamma8*gamma8)-1)/(exp(2*100/4*gamma8*gamma8)+1)*one;
omega9=(exp(2*100/4*gamma9*gamma9)-1)/(exp(2*100/4*gamma9*gamma9)+1)*one;
omega10=(exp(2*100/4*gamma10*gamma10)-1)/(exp(2*100/4*gamma10*gamma10)+1)*one;
omega11=(exp(2*100/4*gamma11*gamma11)-1)/(exp(2*100/4*gamma11*gamma11)+1)*one;
omega12=(exp(2*100/4*gamma12*gamma12)-1)/(exp(2*100/4*gamma12*gamma12)+1)*one;
beta1=omega1*gamma1;
beta2=omega2*gamma2;
beta3=omega3*gamma3;
beta4=omega4*gamma4;
beta5=omega5*gamma5;
beta6=omega6*gamma6;
beta7=omega7*gamma7;
beta8=omega8*gamma8;
beta9=omega9*gamma9;
beta10=omega10*gamma10;
beta11=omega11*gamma11;
beta12=omega12*gamma12;
temp0=X1*beta1+X2*beta2+X3*beta3+X4*beta4+X5*beta5+X6*beta6+u1;
temp1=R1*temp0;
temp2=log(1+exp(temp0));
temp3=R1/2*log(sigmae);
temp4=R1*((log(Y1)-X1*beta7-X2*beta8-X3*beta9-X4*beta10-X5*beta11-X6*beta12-u2)**2)/2/sigmae;
temp5=(omega1+omega2+omega3+omega4+omega5+omega6+omega7+omega8+omega9+omega10+omega11+omega12)/5/100*(1/4*log(500)/2+(4-1)/4*log(100)/2); 
loglik=temp1-temp2-temp3-temp4-temp5;
model Y1~ general(loglik);
random u1 u2 ~normal([0,0],[sigma1,ro_sigma,sigma2]) subject=in;
ods output ParameterEstimates=two.test1 Hessian=two.cov FitStatistics=two.test2;
run;
%mend();
%macro nlm2(data);
proc nlmixed data=&data qpoints=10 hess;
parms gamma1=0.1 gamma2=0.1 gamma3=0.1 gamma4=0.1 gamma5=0.1 gamma6=0.1 gamma7=1 gamma8=0.1 gamma9=0.1 
gamma10=0.1 gamma11=0.1 gamma12=0.1 sigma1=1.0 sigma2=1.2 ro_sigma=0.5 sigmae=1;
bounds sigma1 sigma2 sigmae>0;
omega1=(exp(2*100/4*gamma1*gamma1)-1)/(exp(2*100/4*gamma1*gamma1)+1)*one;
omega2=(exp(2*100/4*gamma2*gamma2)-1)/(exp(2*100/4*gamma2*gamma2)+1)*one;
omega3=(exp(2*100/4*gamma3*gamma3)-1)/(exp(2*100/4*gamma3*gamma3)+1)*one;
omega4=(exp(2*100/4*gamma4*gamma4)-1)/(exp(2*100/4*gamma4*gamma4)+1)*one;
omega5=(exp(2*100/4*gamma5*gamma5)-1)/(exp(2*100/4*gamma5*gamma5)+1)*one;
omega6=(exp(2*100/4*gamma6*gamma6)-1)/(exp(2*100/4*gamma6*gamma6)+1)*one;
omega7=(exp(2*100/4*gamma7*gamma7)-1)/(exp(2*100/4*gamma7*gamma7)+1)*one;
omega8=(exp(2*100/4*gamma8*gamma8)-1)/(exp(2*100/4*gamma8*gamma8)+1)*one;
omega9=(exp(2*100/4*gamma9*gamma9)-1)/(exp(2*100/4*gamma9*gamma9)+1)*one;
omega10=(exp(2*100/4*gamma10*gamma10)-1)/(exp(2*100/4*gamma10*gamma10)+1)*one;
omega11=(exp(2*100/4*gamma11*gamma11)-1)/(exp(2*100/4*gamma11*gamma11)+1)*one;
omega12=(exp(2*100/4*gamma12*gamma12)-1)/(exp(2*100/4*gamma12*gamma12)+1)*one;
beta1=omega1*gamma1;
beta2=omega2*gamma2;
beta3=omega3*gamma3;
beta4=omega4*gamma4;
beta5=omega5*gamma5;
beta6=omega6*gamma6;
beta7=omega7*gamma7;
beta8=omega8*gamma8;
beta9=omega9*gamma9;
beta10=omega10*gamma10;
beta11=omega11*gamma11;
beta12=omega12*gamma12;
temp0=X1*beta1+X2*beta2+X3*beta3+X4*beta4+X5*beta5+X6*beta6+u1;
temp1=R1*temp0;
temp2=log(1+exp(temp0));
temp3=R1/2*log(sigmae);
temp4=R1*((log(Y1)-X1*beta7-X2*beta8-X3*beta9-X4*beta10-X5*beta11-X6*beta12-u2)**2)/2/sigmae;
temp5=(omega1+omega2+omega3+omega4+omega5+omega6+omega7+omega8+omega9+omega10+omega11+omega12)/5/100*(2/4*log(500)/2+(4-2)/4*log(100)/2); 
loglik=temp1-temp2-temp3-temp4-temp5;
model Y1~ general(loglik);
random u1 u2 ~normal([0,0],[sigma1,ro_sigma,sigma2]) subject=in;
ods output ParameterEstimates=two.test1 Hessian=two.cov FitStatistics=two.test2;
run;
%mend();
%macro nlm3(data);
proc nlmixed data=&data qpoints=10 hess;
parms gamma1=0.1 gamma2=0.1 gamma3=0.1 gamma4=0.1 gamma5=0.1 gamma6=0.1 gamma7=1 gamma8=0.1 gamma9=0.1 
gamma10=0.1 gamma11=0.1 gamma12=0.1 sigma1=1.0 sigma2=1.2 ro_sigma=0.5 sigmae=1;
bounds sigma1 sigma2 sigmae>0;
omega1=(exp(2*100/4*gamma1*gamma1)-1)/(exp(2*100/4*gamma1*gamma1)+1)*one;
omega2=(exp(2*100/4*gamma2*gamma2)-1)/(exp(2*100/4*gamma2*gamma2)+1)*one;
omega3=(exp(2*100/4*gamma3*gamma3)-1)/(exp(2*100/4*gamma3*gamma3)+1)*one;
omega4=(exp(2*100/4*gamma4*gamma4)-1)/(exp(2*100/4*gamma4*gamma4)+1)*one;
omega5=(exp(2*100/4*gamma5*gamma5)-1)/(exp(2*100/4*gamma5*gamma5)+1)*one;
omega6=(exp(2*100/4*gamma6*gamma6)-1)/(exp(2*100/4*gamma6*gamma6)+1)*one;
omega7=(exp(2*100/4*gamma7*gamma7)-1)/(exp(2*100/4*gamma7*gamma7)+1)*one;
omega8=(exp(2*100/4*gamma8*gamma8)-1)/(exp(2*100/4*gamma8*gamma8)+1)*one;
omega9=(exp(2*100/4*gamma9*gamma9)-1)/(exp(2*100/4*gamma9*gamma9)+1)*one;
omega10=(exp(2*100/4*gamma10*gamma10)-1)/(exp(2*100/4*gamma10*gamma10)+1)*one;
omega11=(exp(2*100/4*gamma11*gamma11)-1)/(exp(2*100/4*gamma11*gamma11)+1)*one;
omega12=(exp(2*100/4*gamma12*gamma12)-1)/(exp(2*100/4*gamma12*gamma12)+1)*one;
beta1=omega1*gamma1;
beta2=omega2*gamma2;
beta3=omega3*gamma3;
beta4=omega4*gamma4;
beta5=omega5*gamma5;
beta6=omega6*gamma6;
beta7=omega7*gamma7;
beta8=omega8*gamma8;
beta9=omega9*gamma9;
beta10=omega10*gamma10;
beta11=omega11*gamma11;
beta12=omega12*gamma12;
temp0=X1*beta1+X2*beta2+X3*beta3+X4*beta4+X5*beta5+X6*beta6+u1;
temp1=R1*temp0;
temp2=log(1+exp(temp0));
temp3=R1/2*log(sigmae);
temp4=R1*((log(Y1)-X1*beta7-X2*beta8-X3*beta9-X4*beta10-X5*beta11-X6*beta12-u2)**2)/2/sigmae;
temp5=(omega1+omega2+omega3+omega4+omega5+omega6+omega7+omega8+omega9+omega10+omega11+omega12)/5/100*(3/4*log(500)/2+(4-3)/4*log(100)/2); 
loglik=temp1-temp2-temp3-temp4-temp5;
model Y1~ general(loglik);
random u1 u2 ~normal([0,0],[sigma1,ro_sigma,sigma2]) subject=in;
ods output ParameterEstimates=two.test1 Hessian=two.cov FitStatistics=two.test2;
run;
%mend();
%macro nlm4(data);
proc nlmixed data=&data qpoints=10 hess;
parms gamma1=0.1 gamma2=0.1 gamma3=0.1 gamma4=0.1 gamma5=0.1 gamma6=0.1 gamma7=1 gamma8=0.1 gamma9=0.1 
gamma10=0.1 gamma11=0.1 gamma12=0.1 sigma1=1.0 sigma2=1.2 ro_sigma=0.5 sigmae=1;
bounds sigma1 sigma2 sigmae>0;
omega1=(exp(2*100/4*gamma1*gamma1)-1)/(exp(2*100/4*gamma1*gamma1)+1)*one;
omega2=(exp(2*100/4*gamma2*gamma2)-1)/(exp(2*100/4*gamma2*gamma2)+1)*one;
omega3=(exp(2*100/4*gamma3*gamma3)-1)/(exp(2*100/4*gamma3*gamma3)+1)*one;
omega4=(exp(2*100/4*gamma4*gamma4)-1)/(exp(2*100/4*gamma4*gamma4)+1)*one;
omega5=(exp(2*100/4*gamma5*gamma5)-1)/(exp(2*100/4*gamma5*gamma5)+1)*one;
omega6=(exp(2*100/4*gamma6*gamma6)-1)/(exp(2*100/4*gamma6*gamma6)+1)*one;
omega7=(exp(2*100/4*gamma7*gamma7)-1)/(exp(2*100/4*gamma7*gamma7)+1)*one;
omega8=(exp(2*100/4*gamma8*gamma8)-1)/(exp(2*100/4*gamma8*gamma8)+1)*one;
omega9=(exp(2*100/4*gamma9*gamma9)-1)/(exp(2*100/4*gamma9*gamma9)+1)*one;
omega10=(exp(2*100/4*gamma10*gamma10)-1)/(exp(2*100/4*gamma10*gamma10)+1)*one;
omega11=(exp(2*100/4*gamma11*gamma11)-1)/(exp(2*100/4*gamma11*gamma11)+1)*one;
omega12=(exp(2*100/4*gamma12*gamma12)-1)/(exp(2*100/4*gamma12*gamma12)+1)*one;
beta1=omega1*gamma1;
beta2=omega2*gamma2;
beta3=omega3*gamma3;
beta4=omega4*gamma4;
beta5=omega5*gamma5;
beta6=omega6*gamma6;
beta7=omega7*gamma7;
beta8=omega8*gamma8;
beta9=omega9*gamma9;
beta10=omega10*gamma10;
beta11=omega11*gamma11;
beta12=omega12*gamma12;
temp0=X1*beta1+X2*beta2+X3*beta3+X4*beta4+X5*beta5+X6*beta6+u1;
temp1=R1*temp0;
temp2=log(1+exp(temp0));
temp3=R1/2*log(sigmae);
temp4=R1*((log(Y1)-X1*beta7-X2*beta8-X3*beta9-X4*beta10-X5*beta11-X6*beta12-u2)**2)/2/sigmae;
temp5=(omega1+omega2+omega3+omega4+omega5+omega6+omega7+omega8+omega9+omega10+omega11+omega12)/5/100*(4/4*log(500)/2+(4-4)/4*log(100)/2); 
loglik=temp1-temp2-temp3-temp4-temp5;
model Y1~ general(loglik);
random u1 u2 ~normal([0,0],[sigma1,ro_sigma,sigma2]) subject=in;
ods output ParameterEstimates=two.test1 Hessian=two.cov FitStatistics=two.test2;
run;
%mend();
%macro estimate();
libname two 'E:\two-part\code\sas code\bootstrap';
%do jj=1%to 100;
%do kk=0%to 4;
proc import out=two.two_part21r&jj 
datafile="E:\two-part\code\sas code\bootstrap\data\two_part21r&jj..csv" 
dbms=csv
replace;
run;
%nlm&kk(two.two_part21r&jj) 
data two.one;
set _null_;
run;
data two.one;
set two.one two.test1;
KEEP parameter Estimate;
run;
data two.likelihood&jj.L&kk;
set two.test2;
run;

proc transpose data=two.one out=two.one&jj.L&kk;;
id parameter;
run;

data two.one&jj.L&kk;
set two.one&jj.L&kk;
beta1=(exp(50*gamma1*gamma1)-1)/(exp(50*gamma1*gamma1)+1)*gamma1;
beta2=(exp(50*gamma2*gamma2)-1)/(exp(50*gamma2*gamma2)+1)*gamma2;
beta3=(exp(50*gamma3*gamma3)-1)/(exp(50*gamma3*gamma3)+1)*gamma3;
beta4=(exp(50*gamma4*gamma4)-1)/(exp(50*gamma4*gamma4)+1)*gamma4;
beta5=(exp(50*gamma5*gamma5)-1)/(exp(50*gamma5*gamma5)+1)*gamma5;
beta6=(exp(50*gamma6*gamma6)-1)/(exp(50*gamma6*gamma6)+1)*gamma6;
beta7=(exp(50*gamma7*gamma7)-1)/(exp(50*gamma7*gamma7)+1)*gamma7;
beta8=(exp(50*gamma8*gamma8)-1)/(exp(50*gamma8*gamma8)+1)*gamma8;
beta9=(exp(50*gamma9*gamma9)-1)/(exp(50*gamma9*gamma9)+1)*gamma9;
beta10=(exp(50*gamma10*gamma10)-1)/(exp(50*gamma10*gamma10)+1)*gamma10;
beta11=(exp(50*gamma11*gamma11)-1)/(exp(50*gamma11*gamma11)+1)*gamma11;
beta12=(exp(50*gamma12*gamma12)-1)/(exp(50*gamma12*gamma12)+1)*gamma12;
run;
data two.cov&jj.L&kk;
set  two.cov;
run;
%end;
%end;
%mend();
%estimate()
%macro toge();
%do kk=0%to 4;
data two.estL&kk.allr1;  
set _null_;
run;
data two.hesL&kk.allr1;  
set _null_;
run;
data two.liL&kk.allr1;  
set _null_;
run;
%do ii=1%to 100;


proc transpose data=two.one&ii.L&kk out=two.onet&ii.L&kk;
run;
data two.estL&kk.allr1;   
set two.estL&kk.allr1 two.onet&ii.L&kk;   
run;
data two.hesL&kk.allr1;  
set two.hesL&kk.allr1 two.cov&ii.L&kk;  
run;
data two.liL&kk.allr1;  
set two.liL&kk.allr1 two.likelihood&ii.L&kk;  
run;

%end;
%end;
%mend();
%toge()

%macro estimate1();
%do kk=0%to 4;
proc export data=two.estL&kk.allr1
outfile="E:\two-part\code\sas code\bootstrap\Results_save1\estL&kk.allr1.csv"
dbms=csv
label
replace;
run;

proc export data=two.hesL&kk.allr1
outfile="E:\two-part\code\sas code\bootstrap\Results_save1\hesL&kk.allr1.csv"
dbms=csv
label
replace;
run;

proc export data=two.liL&kk.allr1
outfile="E:\two-part\code\sas code\bootstrap\Results_save1\liL&kk.allr1.csv"
dbms=csv
label
replace;
run;
%end;
%mend();
%estimate1();

