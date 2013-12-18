% train_SNB    Calculate log-likelihood: stationary naive Bayes
%
%   LOGL = train_SNB(DATA,D) for DATA an array.
%   Outputs stationary log-likelihood LOGL binned over range D.
%   Assumes columns of DATA are independent input streams.
%
%   LOGL = train_SNB(DATA,D,NS) includes Gaussian smoothing of width NS.

function logl = train_SNB(x,d,ns)

% defaults if not set
if ~exist('ns','var'); ns = 0; end
%keyboard
% dimensions
[nx nn] = size(x); rn = 1:nn;

%x = diff(x);

% if d a vector, then repeat into array
[nd nnd] = size(d); 
if nnd==1; d = repmat(d,[1,nn]); end
%keyboard
% calculate likelihood
for n = rn
    h(:,n) = histc(x(:,n),d(:,n)); 
    l(:,n) = h(:,n)/sum(h(:,n));  
end
%keyboard
% smooth likelihood
g = @(n,si) exp(-(-n:n).^2/(2*si^2))/sum( exp((-n:n).^2/(2*si^2)) );   
if ns~=0
    for n = rn
        ln = filter(g(ns,1),1,l(:,n)); 
        ln(1:ns) = []; ln(nd) = 0; 
        l(:,n) = ln/sum(ln); 
    end
end

% take logs
logl = log(l+eps);
%keyboard
end