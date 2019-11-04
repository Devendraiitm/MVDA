function [A,P] = fastNCA(Z,Astruct,p)
%fastNCA code
% Step1 :finding W* first lets name it just W.
[u,s,v]= svd(Z,'econ');
W=u(:,1:p);
[n,~]=size(Z);
% find S that is a projection matrix, but first has to calculate Wr
% doing here for just 1 column that is k=1 but have to do it for all three
% columns of A
a=zeros(n,p);
for k=1:p
    [Wc,Wr]=rearrange(W, Astruct, k);
    [u1,s1,v1]=svd(Wr);
    S=v1(:,p);
    [u2,s2,v2]=svd(Wc*S);
    % take j representing j no. of nonzero entries have to assign those
    % only as if try to assign other then it will throw error because size
    % of u2 is different
    j=size(u2(:,1));
    for l=1:j
        a(l,k)=u2(l,1);
    end
end
[A]=reconstitute(a,Astruct);
P=pinv(A)*Z;
end

