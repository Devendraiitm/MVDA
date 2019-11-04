clc
clear variables
load ncadata
B=pureabs;
C=measabs;
D=C/(176)^0.5;
%above statement calculates the scaled data because svd we nee covariance
%of yy'/N , so introduce N^0.5 here only
%[u,s,v]=svds(D);
[u,s,v]=svd(D,0);
% denoising using pca
% caluculating the percentage of variation captured
% as the number of independent components is 3 therefore 3 principal
% directions only
var_capt=(((s(1,1))^2)+((s(2,2))^2)+((s(3,3))^2))/trace(s*s')
% in this statement calculating denoised estimate using first three
% principal directions
Denoised_data=(176^0.5)*u(:,1:3)*s(1:3,:)*v';
disp("Denoised data");
disp(Denoised_data)
%no of independent componensts=3 given
Astruct=[1 1 0;1 0 1;0 1 1;1 0 1; 1 1 0;1 0 1;0 1 1];
% now u can use denoised_data=u1*s1*v1' and group it therefore u1*s1=Apca
% and v1 is S matrix lets say this as matrix p.
% calculating Apca
%scaling is important do not miss it
Apca=(176^0.5)*u(:,1:3)*s(1:3,1:3);
% now we have got some value of Apca and we have to adjust it so that Apca
% and Astruct have same structure for this. Z=AS therefore Z=(A1*M)*(M^-1*P)
% now we have to calculate that M for which A1*M and Astruct have the same
% structure.
% here M is multiplication matrix
M=[];
m=3;
for i =1:m
    Ind_Zero_Ast=find(~Astruct(:,i));
    %Apca(Ind_Zeros_Ast,:)*M(:,i)=0 solve this to get M(:,i)
    A_1=Apca(Ind_Zero_Ast,:);
    %can use svd to solve the equation Ax=0 as x is the null space of A
    %therefore eigenvectors corresponding zeros singular value of A are the
    %nullspace od A
    [u1,s1,v1]=svd(A_1);
    % do not do economic svd because array size will be different
    M=[M;v1(:,m)];
    A_1*v1(:,m) ;       
end
% reshaping M matrix
M=reshape(M,[3,3]);
disp("Rotation matrix");
disp(M);
% using this M calculating M^-1*P
P=v(:,1:3)';
% getting A_new whose structure is similar to Astruct and corresponding
% true component matrix.
P_new=pinv(M)*P;
A_new=Apca*M;


%calculating the correlation between true component and calculated S ie P
% here we have calculated the overall correlation
correlation1=[];
for i=1:3
c1=cov(B(i,:)');
c2=cov(P_new(i,:)');
c3=cov(B(i,:)',P_new(i,:)');
correlation1=[correlation1;c3(1,2)/sqrt(c1(1,1)*c2(1,1))];
end
disp("Correlation matrix using only PCA");
disp(correlation1);
i=1:176;
figure
plot(i,P_new(1,:))
hold on
plot(i,P_new(2,:))
plot(i,P_new(3,:))
xlabel("x");
ylabel("Absorption")
title("Extracted pure components spectra using PCA ");
figure
plot(i,pureabs(1,:))
hold on
plot(i,pureabs(2,:))
plot(i,pureabs(3,:))
xlabel("x");
ylabel("Absorption")
title("Pure components spectra");
% doing part1c
% using NCA toolbox
[A,P,iter,ss]= gnca_fast(D,A_new,P_new);
correlation2=[];
for i=1:3
c1=cov(B(i,:)');
c2=cov(P(i,:)');
c3=cov(B(i,:)',P(i,:)');
correlation2=[correlation2;c3(1,2)/sqrt(c1(1,1)*c2(1,1))];
end
disp("Correlation matrix after applying NCA using the toolbox")
disp(correlation2);
i=1:176;
figure
plot(i,P(1,:))
hold on
plot(i,P(2,:))
plot(i,P(3,:))
xlabel("x");
ylabel("Absorption")
title("Extracted pure components spectra using NCA");



