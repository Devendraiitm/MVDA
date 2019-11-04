clear variables
clc
% wish to factorize Z=A*P 
% calculating W*
load ncadata.mat
Z=measabs;
%specifying Astruct
Astruct=[1 1 0;1 0 1;0 1 1;1 0 1; 1 1 0;1 0 1;0 1 1];
% specifying rank of p
p=3;
% calling the function fastNCA written by me
[A,P]=fastNCA(Z,Astruct,p);
% doing analysis
i=1:176;
plot(i,P(1,:))
hold on
plot(i,P(2,:))
plot(i,-P(3,:))
xlabel("x");
ylabel("Absorption")
title("Extracted Pure components spectra using fast NCA");
figure
plot(i,pureabs(1,:))
hold on
plot(i,pureabs(2,:))
plot(i,pureabs(3,:))
xlabel("x");
ylabel("Absorption")
title("Pure components specta");
% it is observed P(3,:) entries are all negative making them positive
P(3,:)=-P(3,:);

% calculating correlations
correlation1=[];
for i=1:3
c1=cov(pureabs(i,:)');
c2=cov(P(i,:)');
c3=cov(pureabs(i,:)',P(i,:)');
correlation1=[correlation1;c3(1,2)/sqrt(c1(1,1)*c2(1,1))];
end




