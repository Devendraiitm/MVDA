function [A] = reconstitute(Amix, Astruct)
% This functions rearranges elements of Amix such that it matches the
% structure of Astruct.  Each column of Amix contains the non-zero elements
% followed by zero elements.  The number of non-zero elements of each
% column of Amix should correspond to the number of non-zeros in
% corresponding column of Astruct
[nsamples, nvar] = size(Amix);
A = [];
for k = 1:nvar
    nzest = length(find(Amix(:,k)));
    nzind = find(Astruct(:,k));
    nztrue = length(nzind);
    if ( nzest ~= nztrue )
        disp('Number of non zeros in Amix and Astruct do not match for column ',k)
        return
    else
      temp = zeros(nsamples,1);
      count = 0;
      for i = 1:nztrue
          count = count + 1;
          temp(nzind(i)) = Amix(count,k);
      end
      A = [A temp];
    end
end