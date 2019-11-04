function [Zc Zr] = rearrange(Z, Astruct, k)
% This function extracts the rows of matrix Z according to the non-zero
% and zero elements of k'th column of Astruct.  Zc contains the rows of Z
% corresponding to the non-zero elements of column k of Astruct and Zr
% contains remaining rows of Z 
nzind = find(Astruct(:,k));
zind = find(~Astruct(:,k));
Zc = Z(nzind,:);
Zr = Z(zind,:);