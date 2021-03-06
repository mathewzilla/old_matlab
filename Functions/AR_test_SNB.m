% AR_test_SNB    Naive Bayes classification, modified for publication in Evans et al 2014
% 'State of the art in whisker-based object localisation'.
%
%   LOGP = test_NB(DATA,D,LOGL) for DATA an array.
%   Outputs log-posterior LOGP, for stationary LOGL binned over D
%   Assumes columns of DATA are independent data streams.

function logp = AR_test_SNB(x,d,logl)

% dimensions
[nx nn] = size(x); rx = 1:nx; rn = 1:nn;

% check arrays of correct size
if nn~=size(logl,2); disp('incompatible dimensions'); return; end
%x = diff(x);

% if d a vector, then repeat into array
[nd nnd] = size(d); 
if nnd==1; d = repmat(d,[1,nn]); end

% find values of distn 
for n = rn; [ig,bx(:,n)] = histc(x(:,n),d(:,n)); end
bx(bx==0) = 1;

% calculate posterior
for n = rn
    cp_n(n) = sum( logl(bx(:,n),n) ); 
end
logp = sum(cp_n); % naive over rn

end
