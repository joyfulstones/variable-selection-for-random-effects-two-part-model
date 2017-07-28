warning('off');
n=100;
an=n/4;
lambda=linspace(log(n)/2/n,log(5*n)/2/n,5);
Gcv=zeros(1,5);
dg=zeros(16,1);
est=zeros(28,1);
likelihood=zeros(1,1);
hessian=zeros(16,16);
SEE_beta=zeros(1,12);%SEE of beta
SEE_sigma=zeros(1,4);%SEE of sigma
for i=0:4
    temp1=strcat('estl',num2str(i),'all.csv');
    temp2=strcat('lil',num2str(i),'all.csv');
    temp3=strcat('hesl',num2str(i),'all.csv');
    num1=csvread(temp1,1,1);
    est=[est,num1];
    num2=csvread(temp2,1,1);
    likelihoodtemp=num2(1:4:end,:);
    likelihood=[likelihood,likelihoodtemp];
    num3=csvread(temp3,1,2);
    hessian=[hessian,num3];
end
est=est(:,2:6);
likelihood=likelihood(:,2:6)./(-2);
hessian=hessian(:,17:1:end);
est=est.*(abs(est)>1e-4);
for i=1:1
    for j=1:5
        temp5=(((8.*an+32.*an.^2.*est((i-1)*28+1:(i-1)*28+12,j).^2).*exp(2.*an.*est((i-1)*28+1:(i-1)*28+12,j).^2)...
            + (8.*an-32.*an.^2.*est((i-1)*28+1:(i-1)*28+12,j).^2).*exp(4.*an.*est((i-1)*28+1:(i-1)*28+12,j).^2))...
            .\(exp(2.*an.*est((i-1)*28+1:(i-1)*28+12,j).^2)+1).^3).*n.*lambda(j);
        hessian(16*(i-1)+1:16*i,16*(j-1)+1:16*j)=-hessian(16*(i-1)+1:16*i,16*(j-1)+1:16*j)+diag([temp5',0,0,0,0]);
        dg(1:12)=8.*an.*exp(2.*an.*est((i-1)*28+1:(i-1)*28+12,j).^2)./(exp(2.*an.*est((i-1)*28+1:(i-1)*28+12,j).^2)+1).^2;
        dg(1:12)=dg(1:12).*(abs(est((i-1)*28+1:(i-1)*28+12,j))>1e-4);
        sigmatheta=diag(dg);
        elambda=trace(inv(hessian(16*(i-1)+1:16*i,16*(j-1)+1:16*j)-n.*lambda(j).*sigmatheta)*hessian(16*(i-1)+1:16*i,16*(j-1)+1:16*j));
        likelihood(i,j)=likelihood(i,j)+n.*lambda(j).*sum(1-2./(exp((est((i-1)*28+1:(i-1)*28+12,j)).^2.*2.*an)+1));
        Gcv(i,j)=-likelihood(i,j)./n./((1-elambda./n).^2);
    end
end
[mem,pos]=max(Gcv);
est=est(:,pos);
betaest=est(17:28,1);% estimation of beta
sigmaest=est(13:16,1);%estimation of sigma
for i=1:1
    est=zeros(2800,1);
    likelihood=zeros(100,1);
    hessian=zeros(1600,16);
    Gcv=zeros(100,5);
    for k=0:4
        temp1=strcat('estl',num2str(k),'allr',num2str(i),'.csv');
        temp2=strcat('lil',num2str(k),'allr',num2str(i),'.csv');
        temp3=strcat('hesl',num2str(k),'allr',num2str(i),'.csv');
        num1=csvread(temp1,1,1);
        est=[est,num1];
        num2=csvread(temp2,1,1);
        likelihoodtemp=num2(1:4:end,:);
        likelihood=[likelihood,likelihoodtemp];
        num3=csvread(temp3,1,2);
        hessian=[hessian,num3];
    end
    est=est(:,2:6);
    likelihood=likelihood(:,2:6)./(-2);
    hessian=hessian(:,17:1:end);
    est=est.*(abs(est)>1e-4);
    for L=1:100
        for j=1:5
            temp5=(((8.*an+32.*an.^2.*est((L-1)*28+1:(L-1)*28+12,j).^2).*exp(2.*an.*est((L-1)*28+1:(L-1)*28+12,j).^2)...
                + (8.*an-32.*an.^2.*est((L-1)*28+1:(L-1)*28+12,j).^2).*exp(4.*an.*est((L-1)*28+1:(L-1)*28+12,j).^2))...
                .\(exp(2.*an.*est((L-1)*28+1:(L-1)*28+12,j).^2)+1).^3).*n.*lambda(j);
            hessian(16*(L-1)+1:16*L,16*(j-1)+1:16*j)=-hessian(16*(L-1)+1:16*L,16*(j-1)+1:16*j)+diag([temp5',0,0,0,0]);
            dg(1:12)=8.*an.*exp(2.*an.*est((L-1)*28+1:(L-1)*28+12,j).^2)./(exp(2.*an.*est((L-1)*28+1:(L-1)*28+12,j).^2)+1).^2;
            dg(1:12)=dg(1:12).*(abs(est((L-1)*28+1:(L-1)*28+12,j))>1e-4);
            sigmatheta=diag(dg);
            elambda=trace(inv(hessian(16*(L-1)+1:16*L,16*(j-1)+1:16*j)-n.*lambda(j).*sigmatheta)*hessian(16*(L-1)+1:16*L,16*(j-1)+1:16*j));
            likelihood(L,j)=likelihood(L,j)+n.*lambda(j).*sum(1-2./(exp((est((L-1)*28+1:(L-1)*28+12,j)).^2.*2.*an)+1));
            Gcv(L,j)=-likelihood(L,j)./n./((1-elambda./n).^2);
        end
    end
    [mem,pos]=max(Gcv');
    temp6=reshape(repmat(pos,12,1),1200,1);
    temp7=reshape(repmat(1:100,12,1),1200,1);
    temp8=reshape(repmat((17:28)',100,1),1200,1);
    temp9=2800*(temp6-1)+28*(temp7-1)+temp8;
    betaest_r=est(temp9);
    temp10=reshape(repmat(pos,4,1),400,1);
    temp11=reshape(repmat(1:100,4,1),400,1);
    temp12=reshape(repmat((13:16)',100,1),400,1);
    temp13=2800*(temp10-1)+28*(temp11-1)+temp12;
    sigmaest_r=est(temp13);
    betaestr_matrix=(reshape(betaest_r,12,100))';
    sigmaestr_matrix=(reshape(sigmaest_r,4,100))';
    SEE_beta(i,:)=sqrt(var(betaestr_matrix));
    SEE_sigma(i,:)=sqrt(var(sigmaestr_matrix));
end
conf_beta=[betaest-1.96*SEE_beta',betaest+1.96*SEE_beta'];%confidence interval of beta
conf_sigma=[sigmaest-1.96*SEE_sigma',sigmaest+1.96*SEE_sigma'];%confidence interval of sigma





        


