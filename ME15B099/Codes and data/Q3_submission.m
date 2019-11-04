clear variables
clc
load Inorfull.mat
% no. of pure species is known in this case which=3
%setting negative absorbance in DATA matrix =0
Z=DATA;
[nm,nw]=size(DATA);
for i=1:nm
    for j=1:nw
        if Z(i,j)<0
            Z(i,j)=0;
        end
    end
end
% so Z contains only non negative absorption data
% calculating the initial estimate of NMF using PCA
[u,s,v]=svd(Z/nm^0.5);
% Z=A*P therefore using first 3 pc can be 
%written as A=scacle*u1*s1 and P=v1'
A_pca=(nm^0.5)*u(:,1:3)*s(1:3,1:3);
P_usingpca=v(:,1:3);
A_pca=abs(A_pca);
P_usingpca=abs(P_usingpca);
% in above matlab code estimated non negative matrices are calculated using
% pca, now calculating on the mixture of all types that is 26
Z_new=Z(1:5:130,:);
[u1,s1,v1]=svd(Z_new/26^0.5);
A_pca_new=(26^0.5)*u1(:,1:3)*s1(1:3,1:3);
P_usingpca_new=v1(:,1:3)';
A_pca_new=abs(A_pca_new);       
P_usingpca_new=abs(P_usingpca_new);
disp("Mixing matrix A using PCA");
disp(A_pca_new)

%Now calculating the NMF using the function nmf
[A_nmf,P_nmf]=nmf(Z_new,A_pca_new,P_usingpca_new,0.01,3,10000);

% calculating correlation and determining which components are extracted
% well
B=[PureCo;PureCr;PureNi];
correlation1=[];
% because do not know which correspond to which so have to iteratee over
% all the possible combinations
for i=1:3
    for j=1:3
        c1=cov(B(i,:)');
        c2=cov(P_nmf(j,:)');
        c3=cov(B(i,:)',P_nmf(j,:)');
        correlation1=[correlation1;c3(1,2)/sqrt(c1(1,1)*c2(1,1))];
    end
end
disp("........................................................")
disp("Correlation matrix using NMF without averaging");
correlation1=reshape(correlation1,[3,3]);
disp(correlation1)
% from correlation matrix can identify which corresponds to Ni,Cr,Co
figure
plot(WAV,P_nmf(1,:))
hold on
plot(WAV,P_nmf(2,:))
plot(WAV,P_nmf(3,:))
xlabel("Wavelengths");
ylabel("Absorption")
title("Extracted pure components spectra using NMF without averaging");
figure
plot(WAV,PureCr)
hold on
plot(WAV,PureCo)
plot(WAV,PureNi)
legend('Cr','Co','Ni')
xlabel("Wavelengths");
ylabel("Absorption")
title("Pure components absorption spectra");
hold off

% now using average and doing the calculation
for i =0:25
   Z_new(i+1,:)=(Z(5*i+1,:)+Z(5*i+2,:)+Z(5*i+3,:)+Z(5*i+4,:)+Z(5*i+5,:))/5;
end
[u1,s1,v1]=svd(Z_new/26^0.5);
A_pca_new=u1(:,1:3)*s1(1:3,1:3);
P_usingpca_new=(26^0.5)*v1(:,1:3)';
A_pca_new=abs(A_pca_new);       
P_usingpca_new=abs(P_usingpca_new);
[A_nmf,P_nmf]=nmf(Z_new,A_pca_new,P_usingpca_new,0.0001,5000,10000);
B=[PureCo;PureCr;PureNi];
correlation2=[];
for i=1:3
    for j=1:3
        c1=cov(B(i,:)');
        c2=cov(P_nmf(j,:)');
        c3=cov(B(i,:)',P_nmf(j,:)');
        correlation2=[correlation2;c3(1,2)/sqrt(c1*c2)];
    end
end
 correlation2=reshape(correlation2,[3,3]);
 disp("........................................................")
 disp("Correlation matrix using NMF after averaging");
 disp(correlation2)
figure
plot(WAV,P_nmf(1,:))
hold on
plot(WAV,P_nmf(2,:))
plot(WAV,P_nmf(3,:))
xlabel("Wavelengths");
ylabel("Absorption")
title("Extracted pure components spectra using NMF and averaging");
figure
plot(WAV,PureCo)
hold on
plot(WAV,P_nmf(2,:))
legend('Pure Co','Estimated Co');
xlabel("Wavelengths");
ylabel("Absorption")
hold off
figure
plot(WAV,PureCr)
hold on
plot(WAV,P_nmf(3,:))
legend('Pure Cr','Estimated Cr');
xlabel("Wavelengths");
ylabel("Absorption")
hold off
figure
plot(WAV,PureNi)
hold on
plot(WAV,P_nmf(1,:))
legend('Pure Ni','Estimated Ni');
xlabel("Wavelengths");
ylabel("Absorption")
hold off


