tic;
s=RandStream('mt19937ar','seed','shuffle'); 
RandStream.setGlobalStream(s);
beta1=[0.5,1,0,0,0,0];
beta2=[0.5,0,0,0,0,0.5];
n=100;
sigma1=1;
sigma2=1.2;
ro=0.5;
B=1;
for i=1:B
one=ones(5*n,1);
[in,R1,X1,X2,X3,X4,X5,X6,Y1]=data1(beta1,beta2,n,sigma1,sigma2,ro);
temp=[one,in,R1,X1,X2,X3,X4,X5,X6,Y1];
columns ={'one','in','R1','X1','X2','X3','X4','X5','X6','Y1'};
data=table(one,in,R1,X1,X2,X3,X4,X5,X6,Y1,'VariableNames', columns);
filename=sprintf('two_part2%d.csv',i);
writetable(data,filename);
idd=zeros(5*n,1);
for j=1:100
   id=unidrnd(n,n,1);
   for k=1:n
       idd((k-1)*5+1:1:k*5,:)=5*(id(k)-1)+[1,2,3,4,5];
   end
   temp2=temp(idd,:);
   one=temp2(:,1);
   R1=temp2(:,3);
   X1=temp2(:,4);
   X2=temp2(:,5);
   X3=temp2(:,6);
   X4=temp2(:,7);
   X5=temp2(:,8);
   X6=temp2(:,9);
   Y1=temp2(:,10);
   in=reshape(repmat(1:n,5,1),5*n,1);
   columns ={'one','in','R1','X1','X2','X3','X4','X5','X6','Y1'};
   data=table(one,in,R1,X1,X2,X3,X4,X5,X6,Y1,'VariableNames', columns);
  filename=sprintf('two_part2%dr%d.csv',i,j);
  writetable(data,filename);
end      
end
toc;