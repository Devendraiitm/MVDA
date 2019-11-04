clc
clear variables
load yeastdata.mat
% you have to load NCA_Yeast_dataset using import data
Z=microarraydata;
[u1,s1,v1]=svd(Z,'econ');
A_est=u1(:,1:33)*s1(1:33,1:33);
P_est=(v1(:,1:33))';
% making the structure of A_est and Astruct same
M=[];
m=33;
for i =1:m
    Ind_Zero_Ast=find(~Astruct(:,i));
    %Apca(Ind_Zeros_Ast,:)*M(:,i)=0 solve this to get M(:,i)
    A_1=A_est(Ind_Zero_Ast,:);
    %can use svd to solve the equation Ax=0 as x is the null space of A
    %therefore eigenvectors corresponding zeros singular value of A are the
    %nullspace od A
    [u1,s1,v1]=svd(A_1);
    % do not do economic svd because array size will be different
    M=[M;v1(:,m)];
    A_1*v1(:,m) ;       
end
% reshaping M matrix
M=reshape(M,[33,33]);
% using this M calculating M^-1*P
P=(M^-1)*P_est;
% calculating A
A=A_est*M;
% calculating variance of all tfs one by one
VAR=zeros(33,1);
indices=zeros(11,1);
for i=1:33
    VAR(i)=var(P(i,:));
end
Max=maxk(VAR,11);
%finding indices of those elements
for i =1:11
    for j =1:33
        if VAR(j)==Max(i)
            indices(i)=j;
        end
    end
end

% using NCA toolbox solving this
[Anca,Pnca,iter,ss]= gnca_fast(Z,A,P);

VAR=zeros(33,1);
indices1=zeros(11,1);
for i=1:33
    VAR(i)=var(Pnca(i,:));
end
Max=maxk(VAR,11);
%finding indices of those elements
for i =1:11
    for j =1:33
        if VAR(j)==Max(i)
            indices1(i)=j;
        end
    end
end
disp("P matrix using NCA");
disp(Pnca);
disp("Indices of the Tfs identified after applying NCA using the toolbox");
disp(indices1);

