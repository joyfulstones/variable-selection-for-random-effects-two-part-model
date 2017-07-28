function [in,R1,X1,X2,X3,X4,X5,X6,Y1]=data1(beta1,beta2,n,sigma1,sigma2,ro)
sigmanu=[sigma1.^2,ro*sigma1*sigma2;ro*sigma1*sigma2,sigma2.^2];
nu=mvnrnd([0,0],sigmanu,n);
nueplson=zeros(5,1);
sigmaeplson=eye(5);%eplson's variance is 1
eplson=mvnrnd(nueplson,sigmaeplson,n);
R=zeros(n,5);
P=zeros(n,5);
S=zeros(n,5);
T=[-2,-1,0,1,2];
nux=zeros(5,1);
sigmax=zeros(5,5);
for i=1:5
    for j=1:5
        sigmax(i,j)=0.5.^(abs(i-j));
    end
end
Xtemp=mvnrnd(nux,sigmax,n);
X1temp=Xtemp(:,1);
X2temp=Xtemp(:,2);
X3temp=Xtemp(:,3);
X4temp=Xtemp(:,4);
X5temp=Xtemp(:,5);
X=zeros(n,30);
Y1=zeros(5*n,1);
R1=zeros(5*n,1);
X1=zeros(5*n,1);
X2=zeros(5*n,1);
X3=zeros(5*n,1);
X4=zeros(5*n,1);
X5=zeros(5*n,1);
X6=zeros(5*n,1);
in=zeros(5*n,1);
for j=1:5
    X(:,(j-1)*6+1:(j*6))=[repmat(T(j),n,1),X1temp,X2temp,X3temp,X4temp,X5temp];
end

for j=1:5
    temp1=(sum((repmat(beta1,n,1).*X(:,(j-1)*6+1:(j*6)))'))';
    temp2=(sum((repmat(beta2,n,1).*X(:,(j-1)*6+1:(j*6)))'))';
    P(:,j)=exp(nu(:,1)+temp1)./(1+exp(nu(:,1)+temp1));
    U=unifrnd(0,1,n,1);
    R(:,j)=(U<=P(:,j));
    S(:,j)=exp(nu(:,2)+eplson(:,j)+temp2);
end

Y=S;%in fact Y=S.*R,but the likelihood does not change when using S instead of Y.
    % when y=0,1./y is meaningless, Using s instead of y, we can prevent that happen
for i=1:n
    Y1((i-1)*5+1:i*5)=(Y(i,:))';
    R1((i-1)*5+1:i*5)=(R(i,:))';
    in((i-1)*5+1:i*5)=repmat(i,5,1);
    X1((i-1)*5+1:i*5)=(X(i,1:6:end))';
    X2((i-1)*5+1:i*5)=(X(i,2:6:end))';
    X3((i-1)*5+1:i*5)=(X(i,3:6:end))';
    X4((i-1)*5+1:i*5)=(X(i,4:6:end))';
    X5((i-1)*5+1:i*5)=(X(i,5:6:end))';
    X6((i-1)*5+1:i*5)=(X(i,6:6:end))';
end
