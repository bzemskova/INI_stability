function [D2,del] = ddz2(z,endvalues)
%%  function D2 = ddz2(z,endvalues)
%   Generate a second derivative matrix for independent variable z using
%   second order centered differences.  z is assumed to be evenly spaced
%   for j=1,N.
%
%   If endvalues = 0, the function is assumed to vanish at j=0 and j=N+1.
%   If endvalues ~=0, the function is not assumed to vanish, and backward
%      and forward differences are calculated at the ends.
%
%

%  check for even spacing
if abs(std(diff(z))/mean(diff(z))) > 1e-6
    disp(['ddz:  values of z not evenly spaced!'])
    D2 = NaN ;
    return
end

del = (max(z)-min(z))/length(z) ;%mean(diff(z)) ;
N = length(z) ;
D2 = zeros(N,N) ;
for ii = 2:N-1
    D2(ii,ii-1) = 1 ;
    D2(ii,ii) = -2 ;
    D2(ii,ii+1) = 1 ;
end
if endvalues == 0
    D2(1,1) = -2 ;
    D2(1,2) = 1 ;
    D2(N,N) = -2 ;
    D2(N,N-1) = 1 ;
else
    D2(1,1) = 2 ;
    D2(1,2) = -5 ;
    D2(1,3) = 4 ;
    D2(1,4) = -1 ;
    D2(N,N-3) = -1 ;
    D2(N,N-2) = 4 ;
    D2(N,N-1) = -5 ;
    D2(N,N) = 2 ;
end
D2 = D2/(del.^2) ;
    

return

